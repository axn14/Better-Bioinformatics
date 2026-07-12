# Better Bioinformatics

A growing, hands-on series on **common mistakes in bioinformatics and how to avoid them**.
Each entry is small, reproducible, and honest: a real pitfall, a minimal demo you can run in a
couple of minutes, the figures that show it, and a short guide on how to not fall for it.

The goal is a reference people can actually use — not "best practices" in the abstract, but
"here is the specific way this goes wrong, here is code that reproduces it, here is the fix."

## How this repo is organised

- **`mistakes/`** — one folder per pitfall. Each is self-contained: a `README.md` write-up, a
  runnable script, and its `figures/`. This is the *case study* layer (code + images + the
  specific inference for that example).
- **Wiki** — the *conceptual* layer: reusable checklists and how-to guidance that generalise
  beyond any single example (e.g. how to make defensible inferences from SHAP). See the
  [Wiki](https://github.com/axn14/Better-Bioinformatics/wiki).

Rule of thumb: if it is "how do I reason about X in general", it goes in the wiki; if it is
"here is a concrete example with data and figures", it goes in `mistakes/`.

## Entries

| # | Mistake | What it covers | Code | Concept (wiki) |
|---|---------|----------------|------|----------------|
| 01 | **SHAP ranking instability under collinearity** | Why the #1 SHAP feature is often a coin-flip when features are correlated, shown on a compositional microbiome pipeline | [`mistakes/01-shap-ranking-instability`](mistakes/01-shap-ranking-instability) | [SHAP Inference Guide](https://github.com/axn14/Better-Bioinformatics/wiki/SHAP-Inference-Guide) |

*More entries coming. Candidate topics: batch effects mistaken for biology, compositional data
analysed without CLR, p-value hunting across contrasts, leakage from pre-split feature
selection, over-reading fold-change without effect size, reference-database contamination,
alpha-diversity compared without rarefaction/normalisation. Suggestions welcome via Issues.*

## Using an entry

```bash
cd mistakes/01-shap-ranking-instability
pip install -r requirements.txt
python shap_instability_demo.py     # regenerates the figures/
```

## Contributing / adding a mistake

See [`CONTRIBUTING.md`](CONTRIBUTING.md). In short: copy `mistakes/TEMPLATE/`, drop in a runnable
script and figures, write the README, and (if there is a generalizable lesson) add a wiki page.

## Ethos

Every demo is reproducible and every claim traces to something the code actually produced.
Synthetic data is labelled as synthetic. Failures are reported loudly — they are the point.

## License
MIT — see [`LICENSE`](LICENSE). Author: Anirudh D. Nair.
