# pigment-gwas — roadmap / next steps

Working to-do list for the compendium. The README has the overview; this is the actionable
version. (Origin: lifted from `SAPPHIRE_genetics`, 2026-07-06.)

## 1. Widen the catalog pull  ← highest priority for PRS / HIrisPlex

Every catalog row is currently flagged `needs_sumstats = True` because the v1 pull drops
fields a polygenic score needs. In `scripts/gwas_catalog.py`, extend `tidy()` (the
classic-download → tidy-row mapper) to keep columns it currently discards from the GWAS
Catalog association download:

- **`INITIAL SAMPLE SIZE` / `REPLICATION SAMPLE SIZE`** → parse the ancestry label(s) + N →
  populate `ancestry` and `sample_size` in the harmonized schema.
- **`95% CI (TEXT)`** → keep the raw text (today it's read only to derive increase/decrease
  direction, then thrown away) → recover SE and the effect unit.
- **`OR or BETA` + the CI unit text** → set `effect_type` to `OR` vs `beta` instead of
  `"unknown"`, so odds ratios and linear betas are never mixed in a score.

Then: add the new fields to `TIDY_COLS`, update the catalog mapping in `analysis/02_combine.qmd`,
re-pull with `--refresh` (archives a new dated version under `output/catalog/versions/`), and
drop `needs_sumstats` for catalog rows that now have trustworthy effect units + ancestry.

## 2. Full summary statistics (where the catalog isn't enough)

For key studies, pull full sumstats (GWAS Catalog sumstats FTP / paper supplements) for proper
per-SNP effect sizes + SE. Fold in as a source (see §3).

## 3. Add more non-catalog sources

One workbook per source under `analysis/sources/`, each writing
`output/sources/<name>_harmonized.csv`; add its tag to `PRECEDENCE` in `analysis/02_combine.qmd`.
The combine step picks it up automatically. Candidates: other African / underrepresented-
ancestry pigmentation GWAS; the HIrisPlex(-S) SNP set (for direct comparison).

## 4. Crawford raw table

`data/sources/crawford2017_table1.tsv` was **reconstructed** from SAPPHIRE's committed union
(the original transcribed TSV was gitignored and not on disk). If the original is recovered,
swap it in and re-run `analysis/sources/crawford2017.qmd`. Also revisit the sign/direction
calibration (Crawford codes one allele/beta convention; calibrate against `rs1426654`-A = light).

## 5. Repo housekeeping

- ~~Add a git remote and push.~~ Done — `origin` = github.com/lasisilab/pigment-gwas.
- ~~Decide docs deployment.~~ Done — deploy via **GitHub Actions**
  (`.github/workflows/publish.yml`); `docs/` is gitignored and rebuilt in CI on every push.
