# pigment-gwas — roadmap / next steps

Working to-do list for the compendium. The README has the overview; this is the actionable
version. (Origin: lifted from `SAPPHIRE_genetics`, 2026-07-06.)

## Done

- **Reproducible build** — project `.venv`, `jupyter: python3` pinned in `_quarto.yml`, deps
  capped in `requirements.txt`, `.python-version`.
- **Deployment** — GitHub Actions builds + deploys the site to Pages on every push to `main`
  (`.github/workflows/publish.yml`); `docs/` is gitignored. Repo is public. Live at
  https://lasisilab.github.io/pigment-gwas/.
- **Scope broadened to full pigmentation** — added eye color + hair color roots (440 → 1072
  catalog loci; Output B 492 → 1124), MC1R red-hair anchor. (HIrisPlex-S scope: skin+eye+hair.)
- **Visualizations** — each workbook ends in an interactive Plotly view (catalog Manhattan,
  Crawford effect-size lollipop, compendium composition) via `scripts/vizhelpers.py`.
- **License + citation** — MIT (`LICENSE`) + `CITATION.cff`.

## 1. Widen the catalog pull  ← highest priority for PRS / HIrisPlex-S

Still the biggest lever: **1068 / 1124 loci are flagged `needs_sumstats`** (not scoring-ready),
because the pull drops fields a polygenic score needs. Today `standard_error` and `sample_size`
are 0/1124, `ancestry` is populated only for the 56 Crawford rows, and 252 catalog rows have an
unknown effect allele (`risk_allele = "?"`). In `scripts/gwas_catalog.py`, extend `tidy()` to keep
columns it currently discards from the GWAS Catalog association download:

- **`INITIAL SAMPLE SIZE` / `REPLICATION SAMPLE SIZE`** → parse ancestry label(s) + N →
  populate `ancestry` and `sample_size`.
- **`95% CI (TEXT)`** → keep the raw text (today read only for increase/decrease direction, then
  thrown away) → recover SE and the effect unit.
- **`OR or BETA` + the CI unit text** → set `effect_type` to `OR` vs `beta` instead of
  `"unknown"`, so odds ratios and linear betas are never mixed in a score.
- Improve the `STRONGEST SNP-RISK ALLELE` parse so fewer effect alleles come through as `?`.

Then: add the new fields to `TIDY_COLS`, update the catalog mapping in `analysis/02_combine.qmd`,
re-pull with `--refresh`, and drop `needs_sumstats` for catalog rows that now have trustworthy
effect units + ancestry. Update the composition viz "scoring-ready" count as it improves.

## 2. Tests + CI gate

No tests exist yet. Add a `pytest` suite for the easy-to-break bits — `vizhelpers.neglog10p`
(string underflow: HERC2 `2E-9237`), `trait_category`, `harmonize.combine` (precedence +
provenance join + dedup), `to_harmonized` (column contract), `assert_anchors` — and run it as a
job that must pass before the Actions deploy.

## 3. HIrisPlex-S comparison + more sources

Add sources under `analysis/sources/`, each writing `output/sources/<name>_harmonized.csv`; add the
tag to `PRECEDENCE` in `analysis/02_combine.qmd` (the combine globs them in automatically).
Priority candidate: the **HIrisPlex-S SNP model** (Chaitanya et al. 2018) as a reference set, then
show overlap with the atlas. Also: other African / underrepresented-ancestry pigmentation GWAS.

## 4. Full summary statistics (where the catalog isn't enough)

For key studies, pull full sumstats (GWAS Catalog sumstats FTP / paper supplements) for proper
per-SNP effect sizes + SE. Fold in as a source (see §3).

## 5. Crawford raw table

`data/sources/crawford2017_table1.tsv` was **reconstructed** from SAPPHIRE's committed union (the
original transcribed TSV was gitignored and not on disk). If the original is recovered, swap it in
and re-run `analysis/sources/crawford2017.qmd`. Also revisit the sign/direction calibration
(Crawford codes one allele/beta convention; calibrate against `rs1426654`-A = light).

## 6. Smaller / open

- **Melanoma-comorbid rows** — 45 loci come in under `"cutaneous melanoma, hair color"` and are
  currently kept as hair-color signal; revisit whether to exclude.
- **Freshness** — a scheduled CI `--refresh` drift-check that opens an issue if the catalog grew,
  so the atlas doesn't silently drift out of date.
