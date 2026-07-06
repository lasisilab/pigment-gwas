#!/usr/bin/env python3
"""
harmonize.py — the ONE shared definition of the combined pigmentation-GWAS schema,
plus the two helpers that keep every source mapping onto it identically.

Why this file exists
--------------------
Each GWAS source (the catalog pull, Crawford 2017, and whatever gets added next) is
mapped onto a single common column set in its own workbook. If each workbook defined
its own columns, they would drift and the final combine would be a mess of renames.
Instead every source workbook and the combine workbook import `HARMONIZED_COLS` and
`to_harmonized()` from here, so the schema is defined exactly once.

The schema is effect-size- and ancestry-aware on purpose: this compendium is meant to
feed an improved polygenic score and HIrisPlex-style comparisons later, so we keep the
fields those need (effect allele, effect size + type, EAF, ancestry, sample size) even
when a given source can't fill them all. `needs_sumstats` flags rows whose effect
units / SE / ancestry are not trustworthy from the source table alone and would need
full summary statistics before scoring — an honest, per-row TODO rather than a silent gap.
"""

import pandas as pd

# ── The shared schema (order matters for output stability) ───────────────────
HARMONIZED_COLS = [
    "rsid",             # rsID (join key across sources)
    "chr",              # chromosome
    "pos",              # position (in `build`)
    "build",            # genome build for pos, e.g. hg38 / hg19
    "effect_allele",    # allele the effect_size is expressed for
    "other_allele",     # the other allele where known ('' if not)
    "eaf",              # effect-allele frequency where reported
    "effect_type",      # 'beta' | 'OR' | 'unknown'
    "effect_size",      # numeric effect (interpret per effect_type)
    "standard_error",   # SE where reported
    "ci_text",          # raw 95% CI text where reported (carries units)
    "pvalue",
    "trait",            # human-readable trait for the association
    "ancestry",         # discovery-sample ancestry/population where known
    "sample_size",      # discovery N where known
    "mapped_gene",      # nearest / mapped gene
    "pubmed",
    "study_accession",
    "source",           # provenance tag, e.g. 'gwas_catalog', 'crawford2017'
    "source_version",   # e.g. the catalog queried_utc, or a paper's year/table
    "needs_sumstats",   # bool: effect units/SE/ancestry need full sumstats to trust
]


def const(value):
    """Wrap a scalar so a mapping can set a whole column to a constant."""
    return lambda df: pd.Series([value] * len(df), index=df.index)


def to_harmonized(df, mapping, source, source_version, needs_sumstats=False):
    """Map an arbitrary source DataFrame onto HARMONIZED_COLS.

    `mapping` is {harmonized_col: spec} where spec is one of:
      * a str        -> take that column from df (missing col -> empty)
      * a callable   -> spec(df) returning a Series (use const(x) for a constant)
    `source` / `source_version` are set on every row. `needs_sumstats` may be a bool
    (applied to all rows) or a callable(df)->Series. Any harmonized column not named
    in `mapping` is filled with '' (or the passed defaults for source/version/flag).
    Returns a DataFrame with exactly HARMONIZED_COLS, in order.
    """
    n = len(df)
    out = pd.DataFrame(index=range(n))
    for col in HARMONIZED_COLS:
        if col == "source":
            out[col] = source
            continue
        if col == "source_version":
            out[col] = source_version
            continue
        if col == "needs_sumstats":
            out[col] = needs_sumstats(df).values if callable(needs_sumstats) else needs_sumstats
            continue
        spec = mapping.get(col)
        if spec is None:
            out[col] = ""
        elif callable(spec):
            out[col] = pd.Series(spec(df)).values
        elif isinstance(spec, str) and spec in df.columns:
            out[col] = df[spec].values
        else:
            out[col] = ""  # named a column that isn't in df -> empty, not a crash
    return out[HARMONIZED_COLS]


def combine(frames, precedence):
    """Union harmonized frames and dedup by rsID, preferring higher-precedence sources.

    `precedence` is a list of source tags, highest priority first. For an rsID present
    in several sources, the highest-precedence row is kept and its `source` is rewritten
    to the '+'-joined set of contributing sources (precedence order), e.g.
    'crawford2017+gwas_catalog' — so multi-source provenance survives dedup, mirroring
    the SAPPHIRE convention. Rows without a real rsID are dropped.
    """
    all_rows = pd.concat(frames, ignore_index=True)
    all_rows = all_rows[all_rows["rsid"].astype(str).str.match(r"^rs\d+$", na=False)].copy()

    rank = {s: i for i, s in enumerate(precedence)}
    big = len(precedence) + 1
    all_rows["_rank"] = all_rows["source"].map(lambda s: rank.get(s, big))

    # provenance: the full set of sources touching each rsID, in precedence order
    def _src_set(sources):
        uniq = sorted(set(sources), key=lambda s: rank.get(s, big))
        return "+".join(uniq)
    src_by_rsid = all_rows.groupby("rsid")["source"].apply(_src_set)

    kept = (all_rows.sort_values("_rank", kind="stable")   # stable: deterministic winner on ties
                    .drop_duplicates("rsid", keep="first")
                    .drop(columns="_rank"))
    kept["source"] = kept["rsid"].map(src_by_rsid)
    kept = kept.sort_values(["chr", "pos"], na_position="last").reset_index(drop=True)
    return kept[HARMONIZED_COLS]
