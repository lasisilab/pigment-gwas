"""Guards on the committed atlas itself — catches a bad pull/combine before it deploys."""
import json
import os

import pandas as pd

from harmonize import HARMONIZED_COLS

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _read(rel):
    return pd.read_csv(os.path.join(REPO, rel), dtype=str, keep_default_na=False)


def test_output_b_schema_exact():
    b = _read("output/pigmentation_gwas_harmonized.csv")
    assert list(b.columns) == HARMONIZED_COLS


def test_output_b_rsids_unique_and_valid():
    b = _read("output/pigmentation_gwas_harmonized.csv")
    assert b["rsid"].is_unique
    assert b["rsid"].str.match(r"^rs\d+$").all()


def test_output_b_identity_columns_populated():
    b = _read("output/pigmentation_gwas_harmonized.csv")
    for col in ("rsid", "source", "source_version"):
        assert (b[col] != "").all(), f"{col} has blanks"


def test_config_anchors_present_in_catalog_and_atlas():
    cfg = json.load(open(os.path.join(REPO, "scripts/traits_pigmentation.json")))
    catalog = set(_read("output/catalog/pigmentation_gwas_catalog.csv")["rsid"])
    atlas = set(_read("output/pigmentation_gwas_harmonized.csv")["rsid"])
    for rs in cfg["anchors"]:
        assert rs in catalog, f"anchor {rs} missing from catalog pull"
        assert rs in atlas, f"anchor {rs} missing from harmonized atlas"
