# Contributing a mistake

Each entry is a small, self-contained case study. To add one:

1. **Copy the template.** `cp -r mistakes/TEMPLATE mistakes/NN-short-slug`
   (use the next number and a short kebab-case slug, e.g. `02-batch-effect-as-biology`).
2. **Add a runnable script.** Prefer synthetic or public data so anyone can run it in a couple
   of minutes. Save figures into the entry's `figures/` folder.
3. **Write the entry `README.md`.** Follow the template sections: the problem, why it happens,
   the demo step by step, results (traced to the script output), how to read the figures, how to
   run, and honest caveats.
4. **If there is a generalizable lesson,** add a wiki page under `wiki/` (concept layer) and
   cross-link it with the entry. Keep "how to reason about X in general" in the wiki and "here is
   a concrete example" in `mistakes/`.
5. **Update the entries table** in the root `README.md`.

## Standards
- Reproducible: fix seeds, pin nothing exotic, list deps in the entry `requirements.txt`.
- Honest: every number in the write-up must come from the script; label synthetic data as such;
  report failures plainly.
- Small: one pitfall per entry. If it needs two figures and 150 lines, that is about right.
- Neutral, useful tone. Teach the fix, do not dunk on anyone.

## Naming
- Entry folders: `NN-short-slug` (`01-shap-ranking-instability`).
- Wiki pages: `Title-With-Hyphens.md` (becomes the wiki link target).
