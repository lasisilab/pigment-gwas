import pandas as pd

from harmonize import to_harmonized, combine, const, HARMONIZED_COLS


def _row(rsid, source, chr="1", pos="100"):
    d = {c: "" for c in HARMONIZED_COLS}
    d.update(rsid=rsid, source=source, chr=chr, pos=pos)
    return d


def test_to_harmonized_contract():
    df = pd.DataFrame({"a": ["rs1", "rs2"], "b": ["1", "2"]})
    mapping = {"rsid": "a", "chr": "b", "build": const("hg38"), "pos": "missing_col"}
    out = to_harmonized(df, mapping, source="src", source_version="v1", needs_sumstats=True)

    assert list(out.columns) == HARMONIZED_COLS            # exact schema, in order
    assert out["rsid"].tolist() == ["rs1", "rs2"]
    assert out["build"].tolist() == ["hg38", "hg38"]        # const()
    assert out["pos"].tolist() == ["", ""]                  # mapped-but-missing col -> '' (no crash)
    assert out["eaf"].tolist() == ["", ""]                  # unmapped col -> ''
    assert out["source"].tolist() == ["src", "src"]
    assert out["needs_sumstats"].tolist() == [True, True]


def test_combine_precedence_and_provenance_join():
    cat = pd.DataFrame([_row("rs1", "gwas_catalog"), _row("rs2", "gwas_catalog", "2", "200")])
    craw = pd.DataFrame([_row("rs1", "crawford2017")])      # rs1 is shared
    out = combine([cat, craw], precedence=["crawford2017", "gwas_catalog"])

    assert set(out["rsid"]) == {"rs1", "rs2"}
    r1 = out[out.rsid == "rs1"].iloc[0]
    assert r1["source"] == "crawford2017+gwas_catalog"     # higher-precedence wins, both recorded
    r2 = out[out.rsid == "rs2"].iloc[0]
    assert r2["source"] == "gwas_catalog"


def test_combine_drops_non_rsid_rows():
    df = pd.DataFrame([_row("rs1", "gwas_catalog"), _row("chr1:123", "gwas_catalog"), _row("", "gwas_catalog")])
    out = combine([df], precedence=["gwas_catalog"])
    assert out["rsid"].tolist() == ["rs1"]


def test_const_length_matches():
    df = pd.DataFrame({"x": [1, 2, 3]})
    s = const("hg19")(df)
    assert s.tolist() == ["hg19", "hg19", "hg19"]
