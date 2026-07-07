import json

import pandas as pd
import pytest

from gwas_catalog import (assert_anchors, tidy, parse_dir, parse_sample, parse_effect,
                          load_config, PullError)


def test_assert_anchors_ok():
    df = pd.DataFrame({"rsid": ["rs1426654", "rs12913832", "rsX"]})
    assert_anchors("skin", df, {"rs1426654": "SLC24A5", "rs12913832": "HERC2"},
                   progress=lambda *_: None)


def test_assert_anchors_missing_raises():
    df = pd.DataFrame({"rsid": ["rs1"]})
    with pytest.raises(PullError):
        assert_anchors("skin", df, {"rs1426654": "SLC24A5"}, progress=lambda *_: None)


def test_parse_dir():
    assert parse_dir("[0.1-0.2] unit increase") == "increase"
    assert parse_dir("blah decrease blah") == "decrease"
    assert parse_dir("") == ""
    assert parse_dir(None) == ""


def test_tidy_row_mapping():
    df = pd.DataFrame([{
        "SNP_ID_CURRENT": "1426654", "SNPS": "rs1426654",
        "STRONGEST SNP-RISK ALLELE": "rs1426654-A", "CHR_ID": "15", "CHR_POS": "48134287",
        "MAPPED_TRAIT": "skin pigmentation", "DISEASE/TRAIT": "Skin pigmentation",
        "RISK ALLELE FREQUENCY": "0.5", "OR or BETA": "1.2", "P-VALUE": "4E-150",
        "MAPPED_GENE": "SLC24A5", "PUBMEDID": "123", "STUDY ACCESSION": "GCST1",
        "95% CI (TEXT)": "[x] increase",
    }])
    r = tidy(df, "skin", "EFO_x", "skin pigmentation")[0]
    assert r["rsid"] == "rs1426654"        # "rs" + SNP_ID_CURRENT when numeric
    assert r["risk_allele"] == "A"         # from STRONGEST SNP-RISK ALLELE
    assert r["direction_raw"] == "increase"
    assert r["pos_hg38"] == "48134287"
    assert r["mapped_gene"] == "SLC24A5"


def test_tidy_unknown_allele_becomes_question():
    df = pd.DataFrame([{"SNP_ID_CURRENT": "", "SNPS": "rs1",
                        "STRONGEST SNP-RISK ALLELE": "rs1-?", "P-VALUE": "1e-9"}])
    assert tidy(df, "a", "b", "c")[0]["risk_allele"] == "?"


def test_load_config_missing_key_raises(tmp_path):
    p = tmp_path / "bad.json"
    p.write_text(json.dumps({"axis": "x"}))   # no roots / anchors
    with pytest.raises(PullError):
        load_config(str(p))


def test_parse_sample():
    anc, n = parse_sample(
        "7,148 Japanese ancestry female cases, 4,034 Japanese ancestry female controls", "NA")
    assert anc == "Japanese"
    assert n == "11182"                                   # 7148 + 4034, commas stripped
    assert parse_sample("1,200 African American individuals", "")[0] == "African American"
    assert parse_sample("NR", "NA") == ("", "")           # unknown ancestry, no N


def test_parse_effect_or_vs_beta():
    t, se = parse_effect("1.45", "[1.33-1.56]")           # bare ratio CI -> OR
    assert t == "OR" and float(se) > 0
    t, se = parse_effect("0.12", "[0.08-0.16] unit increase")   # direction word -> beta
    assert t == "beta" and abs(float(se) - (0.16 - 0.08) / (2 * 1.959964)) < 1e-4
    assert parse_effect("1.2", "")[0] == "unknown"        # value but no CI -> unknown
    assert parse_effect("", "") == ("unknown", "")


def test_load_config_normalizes_roots(tmp_path):
    p = tmp_path / "cfg.json"
    p.write_text(json.dumps({
        "axis": "x", "anchors": {"rs1": "G"},
        "roots": {"EFO_1": ["label one", "note one"], "EFO_2": {"label": "two"}},
    }))
    cfg = load_config(str(p))
    assert cfg.roots["EFO_1"] == {"label": "label one", "note": "note one"}
    assert cfg.roots["EFO_2"]["label"] == "two"
