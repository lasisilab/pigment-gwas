# pigment-gwas

A small, versioned, **harmonized compendium of skin-pigmentation GWAS loci** — built to be
reused across projects (an improved polygenic score, HIrisPlex-style comparisons, locus
lookups). It pulls the NHGRI-EBI GWAS Catalog, folds in key GWAS papers the catalog doesn't
hold, and standardizes everything onto one effect-size- and ancestry-aware schema.

Lifted from the GWAS-pulling code originally developed in the `SAPPHIRE_genetics` project.

## The model: two workbooks, two outputs

```
scripts/gwas_catalog.py   ── engine: config-driven GWAS Catalog pull (child-trait aware, fail-loud, provenance-stamped)
scripts/traits_pigmentation.json ── which trait roots define "pigmentation" (the audit point)
scripts/harmonize.py      ── the ONE shared schema + mapping helpers (imported by every source + the combine)

analysis/01_catalog_pull.qmd     ── run the pull            → output/catalog/pigmentation_gwas_catalog.csv   (Output A, versioned)
analysis/sources/<name>.qmd      ── map one non-catalog src → output/sources/<name>_harmonized.csv           (one per source)
analysis/02_combine.qmd          ── union + dedup           → output/pigmentation_gwas_harmonized.csv         (Output B)
```

- **Output A — the catalog pull, versioned.** Every pull writes a stable
  `output/catalog/pigmentation_gwas_catalog.csv` (+ a `.meta.json` provenance sidecar) and
  archives a dated copy under `output/catalog/versions/`. The `queried_utc` stamp ties each
  row set to the moment it hit the catalog, so you can always reproduce or diff a past pull.
- **Output B — the harmonized compendium.** The catalog plus every `output/sources/*` file,
  unioned and deduplicated by rsID onto the shared schema. For a shared rsID, the
  higher-precedence source wins and `source` records all contributors (e.g.
  `crawford2017+gwas_catalog`).

## Quickstart

```bash
pip install -r requirements.txt

# refresh Output A from the live catalog (skip to use the committed pull)
python scripts/gwas_catalog.py --config scripts/traits_pigmentation.json --refresh

# render the workbooks (also (re)writes the per-source files and Output B)
quarto render
```

Both workbooks default to the committed data (no network), so a fresh checkout renders
identically. `--refresh` / `REFRESH = True` pulls fresh and archives a new version.

## The harmonized schema (`scripts/harmonize.py`)

`rsid, chr, pos, build, effect_allele, other_allele, eaf, effect_type, effect_size,
standard_error, ci_text, pvalue, trait, ancestry, sample_size, mapped_gene, pubmed,
study_accession, source, source_version, needs_sumstats`

`needs_sumstats` is an honest per-row flag: `True` where the source table's effect
units / SE / ancestry aren't trustworthy on their own and would need full summary statistics
before scoring. **Catalog rows are currently flagged** — the catalog's `or_beta` mixes odds
ratios and betas and the v1 pull didn't capture the sample ancestry/CI (see *Roadmap*).

## Add a source

1. Put the source table in `data/sources/`.
2. Copy `analysis/sources/crawford2017.qmd`, point it at the new table, adjust the `mapping`
   passed to `harmonize.to_harmonized(...)`, and write `output/sources/<name>_harmonized.csv`.
3. Add the source tag to `PRECEDENCE` in `analysis/02_combine.qmd` and re-render.

No existing workbook changes — the combine globs `output/sources/*_harmonized.csv`.

## Sources currently included

| source | tag | build | ancestry | notes |
|---|---|---|---|---|
| NHGRI-EBI GWAS Catalog | `gwas_catalog` | hg38 | mixed (per study) | trait-root pull with child traits; `needs_sumstats` flagged |
| Crawford et al. 2017 (PMID 29025994) | `crawford2017` | hg19 | African | Table 1; effect on the derived allele; not in the catalog |

## Roadmap (known limitations, made explicit)

- **Widen the catalog pull** to keep the association file's `INITIAL/REPLICATION SAMPLE`
  (ancestry + N) and the numeric 95% CI, and to disambiguate OR vs beta — then catalog rows
  can drop the `needs_sumstats` flag. Requires a small change to `tidy()` in
  `scripts/gwas_catalog.py` plus a fresh pull.
- **Full summary statistics** for the key studies (catalog sumstats FTP / paper supplements)
  for PRS-grade effect sizes and standard errors.
- More non-catalog sources as needed (each is one workbook).
