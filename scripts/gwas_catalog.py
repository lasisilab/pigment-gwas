#!/usr/bin/env python3
"""
gwas_catalog.py — Portable, config-driven extraction of trait-associated loci from
the live NHGRI-EBI GWAS Catalog.

The pull engine behind this repo's catalog workbook (analysis/01_catalog_pull.qmd).
There are no hardcoded gene or trait lists here — the traits to query are described
by a small JSON config next to this file (scripts/traits_pigmentation.json). Lifted
verbatim (bar output paths) from the SAPPHIRE_genetics project, where it was the
engine behind the 02a/02b query notebooks; kept axis-agnostic so the same code can
pull any trait axis from just a different config.

WHY QUERY BY TRAIT, NOT BY A STATIC GENE LIST
---------------------------------------------
A hardcoded gene list goes stale and hides its provenance. Instead we name a few
ontology ROOT traits (EFO/OBA/MONDO short forms) and let the catalog return every
locus filed under each root AND its child traits (includeChildTraits=true) — so a
study added after we authored the list is still captured. Design = "discover the
roots once, freeze them, pull fresh": completeness at authoring time, auditability
at score time. The frozen roots + anchors live in the config, which is the audit
point (scripts/traits_pigmentation.json).

HOW THE CATALOG IS QUERIED (proven mechanics, do not change)
------------------------------------------------------------
Uses the catalog DOWNLOAD endpoint with child-trait expansion — the mechanism
behind the website's "Download Associations (with child traits)" button:

    GET https://www.ebi.ac.uk/gwas/api/search/downloads
        ?q=shortForm:"<EFO_ID>"&includeChildTraits=true&efo=true&facet=association

Traps that silently return 0 rows (verified 2026-06-25): the REST endpoint
/efoTraits/{id}/associations (does NOT expand children; empty for OBA grouping
nodes) and the query field efoTrait:"..." (use shortForm: instead).

CACHING & REPRODUCIBILITY (the replication feature)
---------------------------------------------------
Every pull is stamped with the UTC time it hit the catalog. Output = a stable CSV
(cfg.out_csv) carrying a `queried_utc` column, a sidecar `<out_csv>.meta.json`
(machine-readable provenance), and a timestamped archive copy under a `versions/`
folder next to the output. Runs DEFAULT TO THE CACHE (no network) so a frozen
analysis re-renders identically; pass --refresh to pull fresh from the API.

FAIL LOUD (so a broken pull can never masquerade as a good one)
---------------------------------------------------------------
  * per-root row counts printed;
  * every root MUST return > 0 rows (else PullError) — a real empty trait is
    indistinguishable from a dropped HTTP response, so both raise;
  * ANCHOR ASSERTIONS: config anchors (e.g. skin rs1426654 SLC24A5, rs12913832
    HERC2) must appear in the output or the run exits non-zero;
  * DRIFT REPORT: mapped-trait values not in the per-axis reviewed snapshot are
    flagged (config.reviewed_traits_file); --strict-drift makes it a hard fail.

Examples
--------
  # cache-first (no network if a cache exists)
  python scripts/gwas_catalog.py --config scripts/traits_pigmentation.json

  # force a fresh, timestamped pull (archived under output/catalog/versions/)
  python scripts/gwas_catalog.py --config scripts/traits_pigmentation.json --refresh

  # hard-fail if the catalog grew a new unreviewed trait under a root
  python scripts/gwas_catalog.py --config scripts/traits_pigmentation.json --refresh --strict-drift

Docs: https://www.ebi.ac.uk/gwas/rest/docs/api  (the downloads endpoint is undocumented)
"""

import argparse
import csv
import datetime
import io
import json
import math
import os
import re
import sys
import time
from dataclasses import dataclass, field

import pandas as pd
# NOTE: `requests` is imported lazily (only the live pull needs it) so that reading the
# committed atlas — the common case for consumers and for cache-first workbook renders —
# needs nothing but pandas.

__version__ = "1.0"

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

DOWNLOAD = "https://www.ebi.ac.uk/gwas/api/search/downloads"
TIMEOUT = 180
RETRIES = 4
SLEEP = 0.3

_SESSION = None


def _session():
    """Lazily create the requests session (imported here, not at module load)."""
    global _SESSION
    if _SESSION is None:
        import requests
        _SESSION = requests.Session()
    return _SESSION

# tidy schema (order matters for output stability)
TIDY_COLS = [
    "axis", "source_trait", "source_short_form", "mapped_trait", "reported_trait",
    "rsid", "chr", "pos_hg38", "risk_allele", "direction_raw", "risk_freq",
    "or_beta", "pvalue", "mapped_gene", "pubmed", "study_accession",
    # widened pull (roadmap #1): effect units, SE recovered from the 95% CI, discovery
    # ancestry + N — so scoring-ready catalog rows can drop the needs_sumstats flag.
    "effect_type", "standard_error", "ci_text", "ancestry", "sample_size",
]

# Broad ancestry descriptors as they appear in the catalog's free-text sample columns.
# Ordered most-specific-first so a compound match (e.g. "African American") is consumed
# before its components ("African") can re-match. Unmatched text -> ancestry stays "".
_ANCESTRIES = [
    "African American", "Afro-Caribbean", "Sub-Saharan African", "African",
    "East Asian", "South East Asian", "Southeast Asian", "South Asian", "Central Asian",
    "Han Chinese", "Chinese", "Japanese", "Korean", "Asian",
    "Hispanic", "Latin American", "Latino", "Native American", "Amerindian",
    "European", "British", "Finnish", "Icelandic", "Sardinian", "Ashkenazi",
    "Middle Eastern", "Indian", "Pakistani", "Oceanian", "Melanesian",
    "Greenlandic", "Inuit", "admixed", "multi-ethnic", "multiethnic",
]

_CI_RANGE = re.compile(r"([0-9]*\.?[0-9]+)\s*[-–]\s*([0-9]*\.?[0-9]+)")
_Z95 = 2 * 1.959964  # 95% CI half-width in SEs (~3.92)


class PullError(RuntimeError):
    """Raised on any fail-loud condition so the run exits non-zero."""


# ── Config ───────────────────────────────────────────────────────────────────
@dataclass
class TraitConfig:
    axis: str                       # short axis id used in output rows, e.g. "skin"
    label: str                      # human title, e.g. "Skin pigmentation"
    roots: dict                     # {short_form: {"label": str, "note": str}}
    anchors: dict                   # {rsid: gene} — must survive the pull
    out_csv: str                    # absolute path (resolved from REPO)
    reviewed_traits_file: str       # absolute path (resolved from REPO)
    description: str = ""
    hp_export: dict = field(default_factory=dict)
    path: str = ""                  # the config file this came from


def _resolve(p):
    return p if os.path.isabs(p) else os.path.join(REPO, p)


def load_config(path):
    """Load and validate a trait-axis config (stdlib json; no pyyaml dependency).

    `roots` may be given as {short_form: {"label":..., "note":...}} or the compact
    {short_form: [label, note]}; both normalize to the dict form. Paths in the
    JSON are relative to the repo root so the same file works from the repo root
    (CLI) or from analysis/ (qmd). Raises PullError with a clear message if a
    required key is missing.
    """
    with open(path) as f:
        raw = json.load(f)
    for key in ("axis", "roots", "anchors"):
        if key not in raw:
            raise PullError(f"config {path}: missing required key '{key}'")

    roots = {}
    for sf, v in raw["roots"].items():
        if isinstance(v, dict):
            roots[sf] = {"label": v.get("label", sf), "note": v.get("note", "")}
        elif isinstance(v, (list, tuple)):
            roots[sf] = {"label": v[0], "note": (v[1] if len(v) > 1 else "")}
        else:
            roots[sf] = {"label": str(v), "note": ""}

    axis = raw["axis"]
    out_csv = raw.get("out_csv", f"output/catalog/gwas_api_{axis}.csv")
    reviewed = raw.get("reviewed_traits_file",
                       f"output/catalog/gwas_api_{axis}_reviewed_traits.txt")
    return TraitConfig(
        axis=axis, label=raw.get("label", axis), roots=roots, anchors=raw["anchors"],
        out_csv=_resolve(out_csv), reviewed_traits_file=_resolve(reviewed),
        description=raw.get("description", ""), hp_export=raw.get("hp_export", {}),
        path=os.path.abspath(path),
    )


# ── Download one root trait (with child traits) ──────────────────────────────
def download_trait(short_form):
    """Fetch all associations for a trait + its children. Returns a DataFrame.

    Distinguishes three outcomes:
      - HTTP/parse failure after retries          -> raise PullError (NOT silent)
      - HTTP 200 but no TSV body                  -> raise PullError (unexpected)
      - HTTP 200 with a valid (possibly 0-row) table -> return the DataFrame
    """
    params = {
        "q": f'shortForm:"{short_form}"',
        "pvalfilter": "", "orfilter": "", "betafilter": "", "datefilter": "",
        "genomicfilter": "", "traitfilter[]": "", "dateaddedfilter": "",
        "facet": "association", "efo": "true", "includeChildTraits": "true",
    }
    import requests  # lazy: only the live pull needs the network library
    last = None
    for attempt in range(RETRIES):
        try:
            r = _session().get(DOWNLOAD, params=params, timeout=TIMEOUT)
            r.raise_for_status()
            if "\t" not in r.text[:2000]:
                raise PullError(f"{short_form}: 200 OK but no TSV body returned "
                                f"(first 200 chars: {r.text[:200]!r})")
            time.sleep(SLEEP)
            return pd.read_csv(io.StringIO(r.text), sep="\t", dtype=str,
                               quoting=csv.QUOTE_NONE, na_filter=False,
                               on_bad_lines="skip")
        except (requests.RequestException, pd.errors.ParserError) as e:
            last = e
            if attempt < RETRIES - 1:
                time.sleep(2 ** attempt)
    raise PullError(f"{short_form}: download failed after {RETRIES} tries ({last})")


# ── Flatten the classic-download rows to one tidy row per association ─────────
def parse_dir(ci_text):
    t = (ci_text or "").lower()
    if "increase" in t:
        return "increase"
    if "decrease" in t:
        return "decrease"
    return ""


def parse_sample(*sample_texts):
    """Parse ancestry label(s) + total N from the catalog's free-text sample columns.

    e.g. "7,148 Japanese ancestry female cases, 4,034 Japanese ancestry female controls"
    -> ("Japanese", "11182"). Sums every comma-formatted count across the INITIAL and
    REPLICATION columns and collects known ancestry descriptors (compound-first so
    "African American" isn't split into "African"). Unknown ancestry -> "" (honest).
    """
    total, anc = 0, []
    for t in sample_texts:
        t = (t or "").strip()
        if not t or t.upper() in ("NA", "NR", "N/A"):
            continue
        for m in re.findall(r"\d[\d,]*", t):
            try:
                total += int(m.replace(",", ""))
            except ValueError:
                pass
        tl = t.lower()
        for a in _ANCESTRIES:
            if a.lower() in tl:
                if a not in anc:
                    anc.append(a)
                tl = tl.replace(a.lower(), " ")   # consume so components don't re-match
    return ", ".join(anc), (str(total) if total else "")


def parse_effect(or_beta, ci_text):
    """Classify OR vs beta and recover a standard error from the 95% CI text.

    The catalog does not flag OR vs beta directly, but its convention disambiguates: a beta
    carries a direction/unit word ("unit increase/decrease") while an odds ratio is a bare
    ratio range "[lo-hi]". Returns (effect_type, standard_error) with effect_type
    "beta" | "OR" | "unknown" and SE on the reported scale (beta: (hi-lo)/3.92; OR:
    (ln hi - ln lo)/3.92, i.e. SE of log-OR). Missing/ambiguous -> ("unknown", "").
    """
    ci = (ci_text or "").strip()
    try:
        val = float((or_beta or "").strip())
    except (TypeError, ValueError):
        val = None
    lo = hi = None
    m = _CI_RANGE.search(ci)
    if m:
        lo, hi = float(m.group(1)), float(m.group(2))
        if hi < lo:
            lo, hi = hi, lo
    cl = ci.lower()
    if any(w in cl for w in ("increase", "decrease", "unit")):
        se = (hi - lo) / _Z95 if lo is not None else None
        return "beta", ("" if se is None else f"{se:.5g}")
    if val is not None and lo is not None and lo > 0:
        return "OR", f"{(math.log(hi) - math.log(lo)) / _Z95:.5g}"
    return "unknown", ""


def tidy(df, axis, short_form, label):
    """Classic download -> tidy rows. Keeps effect units, a CI-derived SE, and discovery
    ancestry/N (the widened pull) alongside the original fields."""
    rows = []
    for _, r in df.iterrows():
        snp = (r.get("SNP_ID_CURRENT") or "").strip()
        snps = (r.get("SNPS") or "").strip()
        rsid = ("rs" + snp) if snp and snp.isdigit() else snps
        sra = r.get("STRONGEST SNP-RISK ALLELE") or ""
        risk = sra.split("-")[-1].strip() if "-" in sra else ""
        risk = risk if re.fullmatch(r"[ACGT]+", risk) else "?"   # accept SNV + indel alleles
        ci = r.get("95% CI (TEXT)") or ""
        or_beta = r.get("OR or BETA") or ""
        effect_type, se = parse_effect(or_beta, ci)
        ancestry, n = parse_sample(r.get("INITIAL SAMPLE SIZE"), r.get("REPLICATION SAMPLE SIZE"))
        rows.append(dict(
            axis=axis, source_trait=label, source_short_form=short_form,
            mapped_trait=r.get("MAPPED_TRAIT") or "",
            reported_trait=r.get("DISEASE/TRAIT") or "",
            rsid=rsid, chr=r.get("CHR_ID") or "", pos_hg38=r.get("CHR_POS") or "",
            risk_allele=risk, direction_raw=parse_dir(ci),
            risk_freq=r.get("RISK ALLELE FREQUENCY") or "",
            or_beta=or_beta, pvalue=r.get("P-VALUE") or "",
            mapped_gene=r.get("MAPPED_GENE") or "",
            pubmed=r.get("PUBMEDID") or "",
            study_accession=r.get("STUDY ACCESSION") or "",
            effect_type=effect_type, standard_error=se, ci_text=ci,
            ancestry=ancestry, sample_size=n,
        ))
    return rows


# ── The reusable core the qmds import ────────────────────────────────────────
def query_traits(roots, axis, progress=print):
    """Query every root (with child traits), tidy, and dedup to one lead row per
    rsID for this axis (lowest p-value wins). Pure: no file writes, no exit — safe
    to call from a read-only notebook render. Raises PullError if a root returns 0.
    """
    rows = []
    for short_form, meta in roots.items():
        label = meta.get("label", short_form)
        df = download_trait(short_form)
        n = len(df)
        progress(f"  {label:34s} {short_form:16s} rows={n}")
        if n == 0:
            raise PullError(
                f"{axis}: root '{label}' ({short_form}) returned 0 rows. "
                "A reviewed root should always carry data — investigate before "
                "trusting this pull (catalog change, retired trait, or API issue).")
        rows.extend(tidy(df, axis, short_form, label))

    out = pd.DataFrame(rows, columns=TIDY_COLS)
    out = out[out.rsid.str.match(r"^rs\d+$", na=False)].copy()
    out["pval_num"] = pd.to_numeric(out["pvalue"], errors="coerce")
    out = (out.sort_values("pval_num", na_position="last")
              .drop_duplicates(["axis", "rsid"], keep="first")
              .drop(columns="pval_num")
              .sort_values(["chr", "pos_hg38"], na_position="last")
              .reset_index(drop=True))
    return out


# ── Fail-loud checks ─────────────────────────────────────────────────────────
def assert_anchors(axis, df, anchors, progress=print):
    present = set(df.rsid)
    missing = {rs: g for rs, g in anchors.items() if rs not in present}
    if missing:
        raise PullError(
            f"{axis}: ANCHOR(S) MISSING from the pull — "
            + ", ".join(f"{rs} ({g})" for rs, g in missing.items())
            + ". The query is not capturing the loci this axis is built on; "
              "do NOT feed this into the PRS until resolved.")
    progress(f"  ✓ {axis} anchors present: "
             + ", ".join(f"{g} {rs}" for rs, g in anchors.items()))


def drift_report(df, reviewed_file, strict, progress=print):
    """Flag mapped traits not seen in the per-axis reviewed snapshot; refresh it.

    First run writes the baseline (review it once); later runs diff against it.
    Tolerates an unwritable directory (e.g. a read-only render host): warns and
    returns [] rather than crashing.
    """
    observed = sorted({t for t in df.mapped_trait if t})
    if not os.path.exists(reviewed_file):
        try:
            os.makedirs(os.path.dirname(reviewed_file), exist_ok=True)
            with open(reviewed_file, "w") as f:
                f.write("# Reviewed mapped-trait snapshot. One per line.\n"
                        "# New traits appearing in a later pull are flagged for review.\n")
                f.write("\n".join(observed) + "\n")
            progress(f"  (baseline reviewed-trait snapshot written → "
                     f"{reviewed_file}; review it once.)")
        except OSError as e:
            progress(f"  (could not write reviewed-trait snapshot: {e})")
        return []
    with open(reviewed_file) as f:
        reviewed = {ln.strip() for ln in f if ln.strip() and not ln.startswith("#")}
    new = [t for t in observed if t not in reviewed]
    if new:
        progress("  ⚠ NEW mapped traits since last review (classify before scoring):\n"
                 + "\n".join(f"      - {t}" for t in new))
        if strict:
            raise PullError(f"{len(new)} unreviewed mapped trait(s); classify them "
                            f"or update {os.path.basename(reviewed_file)}, then rerun.")
    else:
        progress("  ✓ no new mapped traits since last review")
    return new


# ── Cache (timestamped, replication-grade) ───────────────────────────────────
def _utc_now():
    return datetime.datetime.now(datetime.timezone.utc)


def _cache_dir(out_csv):
    # versioned pull history lives in a `versions/` folder next to the stable output
    # (repo-local — derived from out_csv, not anchored to any host project's layout)
    return os.path.join(os.path.dirname(os.path.abspath(out_csv)), "versions")


def load_cache(cfg, progress=print):
    """Load a previously-pulled cache CSV (+ sidecar meta) with no network call."""
    df = pd.read_csv(cfg.out_csv, dtype=str, keep_default_na=False)
    meta = {}
    meta_path = cfg.out_csv + ".meta.json"
    if os.path.exists(meta_path):
        with open(meta_path) as f:
            meta = json.load(f)
    stamp = meta.get("queried_utc") or (df["queried_utc"].iloc[0]
                                        if "queried_utc" in df and len(df) else "unknown")
    progress(f"  using cached pull from {stamp} ({len(df)} loci) — "
             f"{os.path.relpath(cfg.out_csv, REPO)}")
    return df, meta


def write_cache(df, meta, cfg, progress=print):
    """Write stable CSV + sidecar meta.json + a timestamped archive copy.

    The archive (a dated copy of every pull) lands in `versions/` next to out_csv and
    is named from the output's basename — so the pull history is repo-local and reads
    like `<output-stem>_<UTC>.csv`, letting you always go back to a previous pull.
    """
    cdir = _cache_dir(cfg.out_csv)
    os.makedirs(os.path.dirname(cfg.out_csv), exist_ok=True)
    os.makedirs(cdir, exist_ok=True)
    df.to_csv(cfg.out_csv, index=False)
    with open(cfg.out_csv + ".meta.json", "w") as f:
        json.dump(meta, f, indent=2)
    stem = os.path.splitext(os.path.basename(cfg.out_csv))[0]
    stamp = _utc_now().strftime("%Y%m%dT%H%M%SZ")
    archive = os.path.join(cdir, f"{stem}_{stamp}.csv")
    n = 2  # never clobber an existing archive if two pulls land in the same second
    while os.path.exists(archive):
        archive = os.path.join(cdir, f"{stem}_{stamp}_{n}.csv")
        n += 1
    df.to_csv(archive, index=False)
    progress(f"  wrote {os.path.relpath(cfg.out_csv, REPO)} "
             f"(+ .meta.json, archive versions/{os.path.basename(archive)})")


def run_axis(cfg, refresh=False, strict_drift=False, progress=print):
    """End-to-end for one axis. Returns (loci_df, meta).

    Cache-first: if not refresh and cfg.out_csv exists, load it (no network).
    Otherwise pull fresh from the API, run all fail-loud checks, stamp the UTC
    query time onto every row, and return. Does NOT write to disk — the caller
    (CLI, or a notebook that wants to persist) calls write_cache(). This keeps
    run_axis safe to call from a read-only render.
    """
    if not refresh and os.path.exists(cfg.out_csv):
        return load_cache(cfg, progress=progress)

    progress(f"\n=== {cfg.axis} axis ({cfg.label}) — live pull ===")
    df = query_traits(cfg.roots, cfg.axis, progress=progress)
    assert_anchors(cfg.axis, df, cfg.anchors, progress=progress)
    new_traits = drift_report(df, cfg.reviewed_traits_file, strict_drift, progress=progress)

    queried_utc = _utc_now().strftime("%Y-%m-%dT%H:%M:%SZ")
    df = df.copy()
    df["queried_utc"] = queried_utc
    meta = {
        "queried_utc": queried_utc,
        "axis": cfg.axis,
        "label": cfg.label,
        "config_file": os.path.relpath(cfg.path, REPO) if cfg.path else "",
        "roots": list(cfg.roots.keys()),
        "anchors": cfg.anchors,
        "n_loci": int(len(df)),
        "new_unreviewed_traits": new_traits,
        "endpoint": DOWNLOAD,
        "query": 'shortForm:"<root>" includeChildTraits=true efo=true facet=association',
        "source": "NHGRI-EBI GWAS Catalog download endpoint (with child traits)",
        "script": f"scripts/gwas_catalog.py v{__version__}",
    }
    return df, meta


# ── Main ─────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True, action="append", metavar="PATH",
                    help="trait-axis JSON config (repeatable, e.g. scripts/traits_pigmentation.json)")
    ap.add_argument("--out", help="output CSV path; overrides config.out_csv for both "
                                  "the cache-read and the write (only valid with a "
                                  "single --config)")
    ap.add_argument("--refresh", action="store_true",
                    help="pull fresh from the API (default: use the timestamped cache if present)")
    ap.add_argument("--strict-drift", action="store_true",
                    help="hard-fail if any new mapped trait appears since last review")
    args = ap.parse_args()

    if args.out and len(args.config) != 1:
        ap.error("--out is only valid with a single --config")

    summary = ["GWAS-catalog API pull — diagnostic", ""]
    try:
        for cfg_path in args.config:
            cfg = load_config(cfg_path)
            if args.out:                      # redirect cache-read + output to --out
                cfg.out_csv = _resolve(args.out)
            df, meta = run_axis(cfg, refresh=args.refresh,
                                strict_drift=args.strict_drift)
            if args.refresh or not os.path.exists(cfg.out_csv):
                write_cache(df, meta, cfg)
            line = (f"{cfg.axis}: {len(df)} lead SNPs "
                    f"[queried {meta.get('queried_utc', 'cache')}] → "
                    f"{os.path.relpath(cfg.out_csv, REPO)}")
            print(line)
            summary.append(line)
    except PullError as e:
        print(f"\n✗ FAIL-LOUD: {e}", file=sys.stderr)
        sys.exit(2)

    try:
        os.makedirs(os.path.join(REPO, "output"), exist_ok=True)
        with open(os.path.join(REPO, "output", "pull_summary.txt"), "w") as f:
            f.write("\n".join(summary) + "\n")
    except OSError:
        pass


if __name__ == "__main__":
    main()
