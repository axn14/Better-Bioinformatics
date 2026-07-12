# Publishing these pages to the GitHub Wiki

GitHub wikis are stored in a **separate git repository** from the code: `*.wiki.git`. The
markdown files in this `wiki/` folder are the source; publish them one of two ways.

## Option A — paste in the web UI (quickest for a first page)
1. On GitHub, open the repo → **Wiki** tab → **Create the first page** (this initialises the wiki repo).
2. Create a page named exactly **Home**, paste the contents of `Home.md`, save.
3. Create a page named **SHAP-Inference-Guide**, paste `SHAP-Inference-Guide.md`, save.
4. The sidebar comes from a special page named **_Sidebar**; create it and paste `_Sidebar.md`.

Page name = wiki link target. "SHAP-Inference-Guide.md" becomes the page `SHAP-Inference-Guide`,
linked as `[text](SHAP-Inference-Guide)`.

## Option B — push the wiki as a git repo (scales better)
After the wiki has been initialised once (create any page in the UI so `*.wiki.git` exists):
```bash
git clone https://github.com/axn14/Better-Bioinformatics.wiki.git
cd Better-Bioinformatics.wiki
cp ../Better-Bioinformatics/wiki/Home.md .
cp ../Better-Bioinformatics/wiki/SHAP-Inference-Guide.md .
cp ../Better-Bioinformatics/wiki/_Sidebar.md .
git add .
git commit -m "Publish wiki: Home, SHAP Inference Guide, sidebar"
git push
```

> Note: this `wiki/` folder is kept **inside the code repo too**, so the guides are versioned and
> reviewable alongside the code. The GitHub Wiki is just the rendered, browsable copy.
