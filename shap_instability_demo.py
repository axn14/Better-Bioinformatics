"""
SHAP ranking instability on a co-abundant guild (compositional, defensible version)
-----------------------------------------------------------------------------------
Pipeline mirrors a real microbiome workflow so the co-abundance is a *measured* quantity,
not a stipulated one:

  1. Absolute abundances from a log-normal community. A 6-member guild shares a latent
     activity factor L (they are genuinely co-abundant); all other taxa are independent.
  2. Close to relative abundances and draw COUNTS via multinomial(depth) -> real zeros.
  3. Prevalence filter (>=10% of samples), then CLR transform (natural log, pseudocount).
  4. Co-abundance is computed as SPEARMAN on CLR values (the standard definition), so the
     guild's correlation is reported from data, not assumed.
  5. A metabolite target tracks the latent guild activity L. A random forest predicts the
     metabolite from CLR taxa; SHAP is computed; the model is refit on 25 bootstraps and we
     track how each guild member's SHAP rank moves.

Data are SYNTHETIC; species names are illustrative. seed=42. Author: Anirudh D. Nair. MIT.
"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
import matplotlib.gridspec as gridspec
from sklearn.ensemble import RandomForestRegressor
from scipy.stats import spearmanr
import shap

plt.rcParams.update({
    "figure.dpi":130,"savefig.dpi":240,"font.family":"DejaVu Sans","font.size":11,
    "axes.titlesize":13,"axes.titleweight":"bold","axes.labelsize":11.5,
    "axes.edgecolor":"#444444","axes.linewidth":0.9,"xtick.color":"#333333","ytick.color":"#333333",
    "axes.grid":True,"grid.color":"#e6e6e6","grid.linewidth":0.8})
SHAP_CMAP = LinearSegmentedColormap.from_list("shap_rb",["#1E88E5","#B23CE0","#FF0D57"])
DIV_CMAP  = LinearSegmentedColormap.from_list("div",["#2166AC","#f7f7f7","#B2182B"])

rng = np.random.default_rng(42)
N, P, K = 300, 220, 6
DEPTH = 20000                     # reads/sample
GUILD_LOADING = 1.25              # strength of shared latent in guild log-abundance
GUILD_NOISE   = 0.55              # per-member independent log noise
B = 25

DRIVERS = ["Faecalibacterium prausnitzii","Roseburia intestinalis","Agathobacter rectalis",
           "Anaerostipes hadrus","Blautia_A obeum","Ruminococcus_E bromii"]
DRIVERS_ABBR = ["F. prausnitzii","R. intestinalis","A. rectalis","A. hadrus","B. obeum","R. bromii"]
GENERA = ["Bacteroides","Phocaeicola","Alistipes","Parabacteroides","Akkermansia","Bifidobacterium",
          "Prevotella","Dialister","Collinsella","CAG-81","UMGS1071","Dorea","Fusicatenibacter",
          "Gemmiger","Coprococcus","Lachnospira","Oscillibacter","Ruthenibacterium","Flavonifractor",
          "CAG-170","Anaerobutyricum","Mediterraneibacter","UBA738","Butyrivibrio","Blautia","Enterocloster"]
def noise_names(n):
    seen=set(); out=[]
    while len(out)<n:
        nm=f"{GENERA[rng.integers(len(GENERA))]} sp{rng.integers(900000000,916000000)}"
        if nm not in seen: seen.add(nm); out.append(nm)
    return out

# ---- 1. absolute abundances ----
L = rng.normal(size=N)                                   # latent guild activity
log_abund = np.zeros((N, P))
# guild: abundant, share L
guild_base = rng.normal(2.2, 0.3, size=K)                # high baseline -> survives filtering
for k in range(K):
    log_abund[:, k] = guild_base[k] + GUILD_LOADING*L + rng.normal(0, GUILD_NOISE, size=N)
# background: independent, spread of baselines (some rare -> zeros after multinomial)
bg_base = rng.normal(-1.4, 2.2, size=P-K)   # many rare taxa -> realistic zeros
for j in range(P-K):
    log_abund[:, K+j] = bg_base[j] + rng.normal(0, 1.0, size=N)
absA = np.exp(log_abund)
rel = absA/absA.sum(1, keepdims=True)

# ---- 2. counts via multinomial -> real zeros ----
counts = np.vstack([rng.multinomial(DEPTH, rel[i]) for i in range(N)]).astype(float)

# ---- 3. prevalence filter + CLR ----
prev = (counts>0).mean(0)
keep = prev >= 0.10
counts_f = counts[:, keep]
names = np.array(DRIVERS + noise_names(P-K))[keep]
guild_mask = np.array([nm in DRIVERS for nm in names])
comp = counts_f + 0.5                                    # pseudocount
clr = np.log(comp) - np.log(comp).mean(1, keepdims=True)
Xc = pd.DataFrame(clr, columns=names)
Pf = Xc.shape[1]

# ---- 4. measured co-abundance = Spearman on CLR ----
guild_names = [n for n in names if n in DRIVERS]
rho_clr = Xc[guild_names].corr(method="spearman").values
iu = np.triu_indices(len(guild_names),1)
guild_rho = rho_clr[iu]
# background reference: random 400 background pairs
bg_names = [n for n in names if n not in DRIVERS]
bgi = rng.choice(len(bg_names), size=(400,2))
bg_rho = [spearmanr(Xc[bg_names[a]], Xc[bg_names[b]]).statistic for a,b in bgi if a!=b]
sparsity = (counts_f==0).mean()

# ---- 5. metabolite target tracks latent guild activity ----
y = pd.Series(1.5*L + rng.normal(0,0.6,size=N), name="metabolite")

def fit_shap(Xtr,ytr,Xall,seed):
    m=RandomForestRegressor(n_estimators=300,min_samples_leaf=3,max_features="sqrt",
                            random_state=seed,n_jobs=-1).fit(Xtr,ytr)
    return shap.TreeExplainer(m).shap_values(Xall,check_additivity=False)

sv_ref=fit_shap(Xc,y,Xc,42)
mean_abs=np.abs(sv_ref).mean(0)
top_idx=np.argsort(mean_abs)[::-1][:14]
guild_idx=[i for i,n in enumerate(names) if n in DRIVERS]
guild_share_ref=mean_abs[guild_idx].sum()/mean_abs.sum()

ranks=np.zeros((B,Pf)); imp=np.zeros((B,Pf)); top1=[]
for b in range(B):
    idx=rng.integers(0,N,N)
    im=np.abs(fit_shap(Xc.iloc[idx],y.iloc[idx],Xc,int(rng.integers(1e6)))).mean(0)
    imp[b]=im; ranks[b]=pd.Series(im,index=names).rank(ascending=False).values
    top1.append(names[int(np.argmax(im))])
ranks_df=pd.DataFrame(ranks,columns=names)
gshare=imp[:,guild_idx].sum(1)/imp.sum(1)
top1_switch=int((pd.Series(top1).shift()!=pd.Series(top1)).iloc[1:].sum())
n_top1=pd.Series(top1).nunique()
gcorr=np.median([spearmanr(imp[i],imp[j]).statistic for i in range(B) for j in range(i+1,B)])

print(f"retained taxa after >=10% prevalence filter: {Pf} / {P}")
print(f"overall zero fraction (sparsity): {sparsity*100:.1f}%")
print(f"MEASURED guild co-abundance (CLR Spearman): mean {guild_rho.mean():.2f} (range {guild_rho.min():.2f}-{guild_rho.max():.2f})")
print(f"background CLR Spearman: mean {np.mean(bg_rho):.2f} (sd {np.std(bg_rho):.2f})")
print(f"guild SHAP mass (ref): {guild_share_ref*100:.1f}% | bootstrap mean {gshare.mean()*100:.1f}%")
print(f"#1 taxon switches {top1_switch}/{B-1} | distinct #1 {n_top1} | between-run rankcorr {gcorr:.2f}")
print("top-1 taxa:", dict(pd.Series(top1).value_counts()))

# ================= figures =================
def beeswarm(ax,sv,Xdf,feat_idx):
    nm=list(Xdf.columns)
    for row,fi in enumerate(feat_idx[::-1]):
        s=sv[:,fi]; fr=pd.Series(Xdf.iloc[:,fi].values).rank(pct=True).values
        yj=row+(rng.random(len(s))-0.5)*0.7
        ax.scatter(s,yj,c=fr,cmap=SHAP_CMAP,s=13,alpha=0.75,linewidths=0,rasterized=True)
    ax.axvline(0,color="#888",lw=1.0,zorder=0)
    ax.set_yticks(range(len(feat_idx))); ax.set_yticklabels([nm[i] for i in feat_idx[::-1]],fontsize=10,fontstyle="italic")
    ax.set_ylim(-0.7,len(feat_idx)-0.3); ax.set_xlabel("SHAP value (impact on predicted metabolite)")
    ax.grid(axis="y",visible=False)
    for sp in ("top","right"): ax.spines[sp].set_visible(False)

palette=plt.cm.tab10(np.linspace(0,1,K))
txt=(f"guild = {gshare.mean()*100:.0f}% of SHAP mass (stable)\n"
     f"measured co-abundance (CLR Spearman) = {guild_rho.mean():.2f}\n"
     f"#1 taxon changed {top1_switch}/{B-1} resamples\n"
     f"between-run rank corr: median {gcorr:.2f}")

# combined main (A beeswarm, B bump)
fig=plt.figure(figsize=(15,6.4))
gs=gridspec.GridSpec(1,2,width_ratios=[1.12,1.0],wspace=0.42)
axA=fig.add_subplot(gs[0]); beeswarm(axA,sv_ref,Xc,top_idx)
axA.set_title("A  SHAP summary: a co-abundant guild dominates",loc="left")
cbA=fig.colorbar(ScalarMappable(Normalize(0,1),SHAP_CMAP),ax=axA,pad=0.01,fraction=0.045)
cbA.set_ticks([0,1]); cbA.set_ticklabels(["Low","High"]); cbA.set_label("CLR abundance",fontsize=9)
axB=fig.add_subplot(gs[1])
gseries={n:ranks_df[n] for n in guild_names}
for i,(n,ab) in enumerate(zip(DRIVERS,DRIVERS_ABBR)):
    if n in ranks_df: axB.plot(range(1,B+1),ranks_df[n],marker="o",ms=4,lw=1.5,color=palette[i],alpha=0.9,label=ab)
axB.set_ylim(K+0.5,0.5); axB.set_yticks(range(1,K+1))
axB.set_xlabel("Bootstrap resample"); axB.set_ylabel("SHAP rank (1 = top)")
axB.set_title("B  ...but the #1 taxon is a coin-flip",loc="left")
axB.legend(ncol=3,fontsize=8,loc="upper center",bbox_to_anchor=(0.5,-0.14),frameon=False)
axB.text(0.985,0.05,txt,transform=axB.transAxes,ha="right",va="bottom",fontsize=8.5,
         bbox=dict(boxstyle="round,pad=0.4",fc="#f5f5f5",ec="#cccccc"))
for sp in ("top","right"): axB.spines[sp].set_visible(False)
fig.suptitle("SHAP ranking is unstable under collinearity: the guild is real, the ordering is noise",
             fontsize=14,fontweight="bold",y=1.03)
fig.savefig("shap_instability_main.png",bbox_inches="tight"); plt.close(fig)
__import__("shutil").copy("shap_instability_main.png","shap_instability.png")

# standalone beeswarm + bump
f1,a1=plt.subplots(figsize=(8.8,6.4)); beeswarm(a1,sv_ref,Xc,top_idx)
a1.set_title("SHAP summary - synthetic gut-omics model\ntop 14 taxa predicting a metabolite")
cb=f1.colorbar(ScalarMappable(Normalize(0,1),SHAP_CMAP),ax=a1,pad=0.01,fraction=0.045)
cb.set_ticks([0,1]); cb.set_ticklabels(["Low","High"]); cb.set_label("CLR abundance",fontsize=10)
f1.savefig("shap_beeswarm.png",bbox_inches="tight"); plt.close(f1)

f2,a2=plt.subplots(figsize=(8.6,5.2))
for i,(n,ab) in enumerate(zip(DRIVERS,DRIVERS_ABBR)):
    if n in ranks_df: a2.plot(range(1,B+1),ranks_df[n],marker="o",ms=4.5,lw=1.6,color=palette[i],alpha=0.9,label=ab)
a2.set_ylim(K+0.5,0.5); a2.set_yticks(range(1,K+1)); a2.set_xlabel("Bootstrap resample"); a2.set_ylabel("SHAP rank (1 = top)")
a2.set_title("Which guild member ranks #1 flips every resample")
a2.legend(ncol=3,fontsize=8.5,loc="upper center",bbox_to_anchor=(0.5,-0.14),frameon=False)
a2.text(0.985,0.05,txt,transform=a2.transAxes,ha="right",va="bottom",fontsize=9,
        bbox=dict(boxstyle="round,pad=0.4",fc="#f5f5f5",ec="#cccccc"))
for sp in ("top","right"): a2.spines[sp].set_visible(False)
f2.savefig("shap_rank_instability.png",bbox_inches="tight"); plt.close(f2)

# co-abundance verification heatmap (top-14 taxa, CLR Spearman)
shown=[names[i] for i in top_idx]
M=Xc[shown].corr(method="spearman")
f3,a3=plt.subplots(figsize=(7.6,6.6))
im3=a3.imshow(M.values,cmap=DIV_CMAP,vmin=-1,vmax=1)
a3.set_xticks(range(len(shown))); a3.set_xticklabels(shown,rotation=90,fontsize=7.5,fontstyle="italic")
a3.set_yticks(range(len(shown))); a3.set_yticklabels(shown,fontsize=7.5,fontstyle="italic")
a3.set_title("Measured co-abundance (CLR Spearman)\nguild block is genuinely correlated")
f3.colorbar(im3,ax=a3,fraction=0.046,pad=0.02,label="Spearman rho")
f3.savefig("coabundance_check.png",bbox_inches="tight"); plt.close(f3)
print("saved: shap_instability_main.png, shap_beeswarm.png, shap_rank_instability.png, coabundance_check.png")
