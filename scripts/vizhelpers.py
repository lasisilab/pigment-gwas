#!/usr/bin/env python3
"""
vizhelpers.py — shared plotting helpers for the workbook visualizations.

Kept in one place (imported by the qmd workbooks) so the two non-obvious, easy-to-get-wrong
pieces of GWAS-plot data prep are defined exactly once:

  * neglog10p()   — compute -log10(p) from a p-value STRING. The catalog reports p-values like
                    "2E-9237" (HERC2); pd.to_numeric underflows anything below ~1e-308 to 0.0,
                    and -log10(0) = inf, which silently destroys the strongest loci in the atlas.
                    Parsing the exponent out of the string keeps them (HERC2 -> 9236.7).
  * trait_category() — collapse the messy free-text mapped_trait / trait labels (13+ variants,
                    incl. compound strings like "cutaneous melanoma, hair color") into 4 stable
                    pigmentation categories for coloring: Eye / Hair / Skin/UV / Multi.

Plus a genome_layout() helper that lays chromosomes end-to-end for a Manhattan x-axis, and a
shared colorblind-safe palette so the three workbooks look like one atlas.
"""

import re
import numpy as np
import pandas as pd

# ── Shared palettes (colorblind-safe) ────────────────────────────────────────
TRAIT_COLORS = {
    "Eye":     "#4C72B0",   # blue
    "Hair":    "#DD8452",   # orange
    "Skin/UV": "#55A868",   # green
    "Multi":   "#8C8C8C",   # gray
}
# effect-direction palette (Crawford betas; sign only, convention under review)
DIR_POS = "#C44E52"        # + beta (red)
DIR_NEG = "#4C72B0"        # - beta (blue)

# autosomes only ever appear in this pull; X/Y/MT coded for robustness but absent
CHROM_ORDER = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]

_SCI = re.compile(r"^\s*([0-9]*\.?[0-9]+)\s*[eE]\s*([+-]?[0-9]+)\s*$")


def neglog10p(value):
    """-log10(p) from a p-value string/number, underflow-safe.

    "2E-9237" -> 9236.70 (NOT inf); "5e-8" -> 7.30; "0.01" -> 2.0; ""/bad -> NaN.
    Parses the exponent directly so values below the float64 floor (~1e-308) survive
    instead of being crushed to 0.0 by float conversion.
    """
    s = str(value).strip()
    if not s:
        return np.nan
    m = _SCI.match(s)
    if m:
        mant, exp = float(m.group(1)), int(m.group(2))
        return -(np.log10(mant) + exp) if mant > 0 else np.nan
    try:
        v = float(s)
    except ValueError:
        return np.nan
    return -np.log10(v) if v > 0 else np.nan


def trait_category(trait):
    """Collapse a raw mapped_trait / trait label into one of 4 pigmentation categories.

    Rules (a label is scored across three pigmentation *domains*):
      eye   : contains "eye"
      hair  : contains "hair"   (incl. "cutaneous melanoma, hair color" — a hair-color assoc.)
      skin  : skin pigmentation / facial pigmentation / suntan / sunburn / freckle /
              sun sensitivity / melanin
    Two or more domains present -> "Multi" (genuinely multi-trait strings). Otherwise the single
    domain wins. Nothing matched -> "Multi" (rare; keeps the plot honest rather than dropping it).
    """
    t = (trait or "").lower()
    eye = "eye" in t
    hair = "hair" in t
    skin = any(k in t for k in (
        "skin", "facial pigment", "sun", "freckle", "melanin",
    ))  # "skin" + "sun" cover skin pigmentation, suntan, sunburn, skin sensitivity to sun
    domains = eye + hair + skin
    if domains >= 2:
        return "Multi"
    if eye:
        return "Eye"
    if hair:
        return "Hair"
    if skin:
        return "Skin/UV"
    return "Multi"


def genome_layout(df, chr_col="chr", pos_col="pos", gap=20_000_000):
    """Lay chromosomes end-to-end for a Manhattan x-axis.

    Returns (plot_df, ticks, ticklabels) where plot_df is df restricted to rows with a real
    chromosome (1..22, X/Y/MT) and a numeric position, plus an added `_x` cumulative-genome
    coordinate and `_chr` (ordered categorical). `ticks`/`ticklabels` place each chromosome's
    number at the midpoint of its span. Rows with blank/foreign chr or missing pos are dropped
    (report the count separately — never silently). `gap` is the blank space between chromosomes.
    """
    d = df.copy()
    d["_chr"] = d[chr_col].astype(str).str.strip()
    d = d[d["_chr"].isin(CHROM_ORDER)]
    d["_pos"] = pd.to_numeric(d[pos_col], errors="coerce")
    d = d.dropna(subset=["_pos"])
    present = [c for c in CHROM_ORDER if c in set(d["_chr"])]

    offset, cum, ticks, ticklabels = {}, 0, [], []
    for c in present:
        offset[c] = cum
        cmax = d.loc[d["_chr"] == c, "_pos"].max()
        ticks.append(cum + cmax / 2)
        ticklabels.append(c)
        cum += cmax + gap
    d["_x"] = d["_pos"] + d["_chr"].map(offset)
    d["_chr"] = pd.Categorical(d["_chr"], categories=present, ordered=True)
    return d, ticks, ticklabels
