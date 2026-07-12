# Making defensible inferences from SHAP (an omics-focused checklist)

SHAP is a description of **how a model used a feature**, not proof that the feature is causal or
even individually meaningful. This guide is the workflow I now follow so that a SHAP ranking
survives a reviewer, a rerun, and my own second-guessing. It is written with microbiome /
multi-omics data in mind, but most of it generalises.

Companion case study: [01 - SHAP ranking instability](https://github.com/axn14/Better-Bioinformatics/tree/main/mistakes/01-shap-ranking-instability) shows the core failure mode this guide defends against.

---

## Part 0 — Precursors: what must be true before you trust SHAP at all

SHAP explains whatever model you give it, including a model that learned nothing useful or
learned an artifact. Clear these gates first.

1. **The model actually predicts.** Report cross-validated R² / AUC against a permutation null.
   If the model is not meaningfully better than shuffled-label chance, its SHAP values explain
   noise. (In the CRC work, mean regression R² was ~0.05 — that alone caps how hard any single
   SHAP claim can be pushed.)
2. **No leakage.** Feature selection, scaling, and imputation happen **inside** each CV fold, not
   before splitting. A feature that looks all-important is often a leak or a batch label.
3. **Confounders are handled.** Age, sex, BMI, sequencing batch, DNA-extraction kit, study/cohort.
   If a "top taxon" co-varies with batch, SHAP will happily rank the batch signal.
4. **The right explainer for the model.** `TreeExplainer` for tree ensembles, `LinearExplainer`
   for linear models, `KernelExplainer`/`PermutationExplainer` for black boxes. Do not mix.
5. **Feature representation is consistent and appropriate.** For compositional data, explain the
   model in the same space you trained it (e.g. CLR), not raw relative abundance.
6. **You know your perturbation semantics.** Tree SHAP has interventional vs tree-path-dependent
   modes; they answer subtly different questions under correlation. Pick one deliberately and
   state it.

If any gate fails, fix it before reading a single SHAP bar.

---

## Part 1 — The inference workflow (step by step)

1. **Establish predictive signal first.** Lock the CV metric and the permutation null. This is
   the licence to interpret at all.
2. **Compute SHAP with the correct explainer**, on held-out or out-of-fold samples where feasible,
   with a sensible background/reference set.
3. **Cluster correlated features and report group importance.** Compute the feature-feature
   correlation (CLR-Spearman for microbiome), cut it into blocks (e.g. hierarchical clustering at
   |ρ| ≥ 0.6), and sum |SHAP| within each block. The block is your unit of inference.
4. **Quantify stability, don't assume it.** Bootstrap or repeat CV, recompute SHAP, and report:
   between-run rank correlation, and top-k overlap (Jaccard) for the features you plan to name.
   If the #1 changes every resample, you have a tie, not a winner.
5. **Check direction, not just magnitude.** On the beeswarm, high feature value should push the
   prediction consistently one way. Smeared colour means interactions or nonlinearity — inspect
   dependence plots before claiming a direction.
6. **Cross-check against an orthogonal line of evidence.** Differential abundance, simple
   correlation, known biology/enzymatic capacity, or a mechanistic model. Convergence across
   independent methods is what turns a candidate into a claim.
7. **Name candidates only at the supported resolution.** "This guild is associated with the
   metabolite" is often defensible; "*Species X* is the producer" usually is not, unless step 4
   shows X is stably #1 and step 6 backs it mechanistically.

---

## Part 2 — Troubleshooting (symptom → likely cause → fix)

| Symptom | Likely cause | Fix |
|---|---|---|
| Top features reshuffle every run | Collinear features sharing Shapley credit | Group correlated features; report block importance + rank stability |
| One feature dominates implausibly | Leakage or a confound (batch, depth) | Audit the pipeline; add the confound as a covariate or regress it out |
| SHAP ranking disagrees with model performance | Wrong explainer, or additivity broken | Match explainer to model; check `check_additivity`; verify expected value |
| All |SHAP| values tiny / flat | Weak or no predictive signal | Do not interpret; revisit features, sample size, target definition |
| Beeswarm colour not monotonic | Interaction / nonlinearity | Use dependence plots; consider interaction SHAP |
| Global bar vs beeswarm tell different stories | Sign cancellation in the mean | Report mean(|SHAP|), and look at signed distributions |
| TreeExplainer slow or memory-heavy | Large N × features, deep forest | Subsample rows for explanation; cap trees/depth; use path-dependent mode |
| Results flip between interventional and path-dependent | Correlation + perturbation semantics | Choose the mode that matches your question; state it explicitly |
| "Important" taxon is a murine/reference artifact | Reference-database contamination | Sanity-check taxonomy; drop non-credible hits before interpreting |

---

## Part 3 — Follow-ups: what to do after SHAP

1. **Validate candidates orthogonally.** Independent cohort, held-out study, or wet-lab
   (co-culture, isotope tracing). SHAP prioritises hypotheses; it does not confirm them.
2. **Report stability as a first-class result.** Put the rank-correlation / top-k Jaccard in the
   paper next to the importance plot. Reviewers trust honesty about resolution.
3. **Prefer group-level claims; state the resolution limit** explicitly when the data cannot
   separate members of a correlated block.
4. **Run sensitivity analyses.** Different seeds, a second model family, an alternate background
   set. Claims that survive all three are the ones worth defending.
5. **Pre-register or freeze the interpretation plan** where possible, so the "top hit" is not
   chosen after seeing the ranking.

---

## Further reading
- Lundberg & Lee (2017), *A Unified Approach to Interpreting Model Predictions* (SHAP).
- Lundberg et al. (2020), *From local explanations to global understanding with explainable AI for trees* (Tree SHAP).
- Aas, Jullum & Løland (2021), *Explaining individual predictions when features are dependent* (correlated-feature SHAP).
- Molnar, *Interpretable Machine Learning* — chapters on Shapley values and correlated features.
- Gloor et al. (2017), *Microbiome datasets are compositional* (why CLR).

*Author: Anirudh D. Nair. Shared under MIT alongside the demo.*
