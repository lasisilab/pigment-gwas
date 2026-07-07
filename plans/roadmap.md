# pigment-gwas — roadmap / next steps

Working to-do list for the compendium. The README has the overview; this is the actionable
version. (Origin: lifted from `SAPPHIRE_genetics`, 2026-07-06.)

## Done

- **Reproducible build** — project `.venv`, `jupyter: python3` pinned, deps capped, `.python-version`.
- **Deployment** — GitHub Actions builds + deploys to Pages on every push to `main`; `docs/`
  gitignored; repo public. Live at https://lasisilab.github.io/pigment-gwas/.
- **Scope broadened to full pigmentation** — eye + hair color roots added (catalog 440 → 1072),
  MC1R red-hair anchor.
- **Visualizations** — an interactive Plotly view per workbook (catalog Manhattan, Crawford
  effect-size lollipop, compendium composition, HIrisPlex-S coverage) via `scripts/vizhelpers.py`.
- **License + citation** — MIT (`LICENSE`) + `CITATION.cff`.
- **Widened catalog pull** — `tidy()` recovers effect units (OR vs β, from the 95% CI), a
  CI-derived standard error, and discovery ancestry + N. Scoring-ready loci: **56 → 484 / 1124**;
  ancestry 56 → 894; SE 0 → 610. (OR/β is inferred from the CI convention — spot-check for
  high-stakes scoring.)
- **Tests + CI gate** — `pytest` suite in `tests/`, run as a required job before deploy.
- **HIrisPlex-S reference panel** — the 41-SNP forensic eye/hair/skin predictor added as a source
  (`analysis/sources/hirisplex_s.qmd`); **26/41 already in the atlas**, 15 added as reference rows.

## 1. Full summary statistics for the remaining loci  ← now the biggest lever

~640 / 1139 loci still carry `needs_sumstats = True` — catalog rows without a usable CI or
recorded ancestry. For the key studies, pull full summary statistics (GWAS Catalog sumstats FTP /
paper supplements) for proper per-SNP effect sizes + SE, and fold them in as sources. Also
spot-check the heuristic OR-vs-β classification from the widened pull against reported units.

## 2. More non-catalog sources

One workbook per source under `analysis/sources/`, each writing `output/sources/<name>_harmonized.csv`;
add the tag to `PRECEDENCE` in `analysis/02_combine.qmd` (the combine globs them in automatically).
Candidates: other African / underrepresented-ancestry pigmentation GWAS; the HIrisPlex-S model
weights (to go from panel membership to actual prediction).

## 3. Crawford raw table

`data/sources/crawford2017_table1.tsv` was **reconstructed** from SAPPHIRE's committed union. If the
original is recovered, swap it in and re-run `analysis/sources/crawford2017.qmd`. Also revisit the
sign/direction calibration (calibrate against `rs1426654`-A = light).

## 4. Smaller / open

- **Melanoma-comorbid rows** — 45 loci under `"cutaneous melanoma, hair color"` are kept as
  hair-color signal; revisit whether to exclude.
- **Freshness** — a scheduled CI `--refresh` drift-check that opens an issue if the catalog grew.
- **ORCID** in `CITATION.cff` (currently a placeholder).
