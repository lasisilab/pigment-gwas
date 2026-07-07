import math

import pandas as pd

from vizhelpers import neglog10p, trait_category, genome_layout


def test_neglog10p_survives_float_underflow():
    # these underflow float64 to 0.0 -> -log10 would give inf; string parse must keep them
    assert abs(neglog10p("2E-9237") - 9236.699) < 0.01   # HERC2
    assert abs(neglog10p("1E-320") - 320.0) < 1e-6
    assert abs(neglog10p("5e-8") - 7.301) < 0.01         # genome-wide threshold
    assert neglog10p("0.01") == 2.0
    assert neglog10p("1") == 0.0


def test_neglog10p_numeric_input():
    assert abs(neglog10p(5e-8) - 7.301) < 0.01


def test_neglog10p_bad_values_are_nan():
    for bad in ["", "  ", "NR", "0", "-1", "abc"]:
        assert math.isnan(neglog10p(bad)), bad


def test_trait_category_single_domain():
    assert trait_category("eye color") == "Eye"
    assert trait_category("eye colour measurement") == "Eye"
    assert trait_category("hair color") == "Hair"
    assert trait_category("skin pigmentation") == "Skin/UV"
    assert trait_category("suntan") == "Skin/UV"
    assert trait_category("sunburn") == "Skin/UV"
    assert trait_category("freckles") == "Skin/UV"
    assert trait_category("facial pigmentation") == "Skin/UV"
    assert trait_category("melanin measurement") == "Skin/UV"


def test_trait_category_regressions():
    # "skin sensitivity to sun" once fell through to Multi (keyword bug)
    assert trait_category("skin sensitivity to sun") == "Skin/UV"
    # hair-color association reported in a melanoma cohort is hair signal, not skin/multi
    assert trait_category("cutaneous melanoma, hair color") == "Hair"


def test_trait_category_multi_and_empty():
    assert trait_category("eye colour measurement, strand of hair color, skin pigmentation") == "Multi"
    assert trait_category("") == "Multi"
    assert trait_category(None) == "Multi"


def test_genome_layout_orders_numerically_and_drops_unplaceable():
    df = pd.DataFrame({
        "chr": ["2", "10", "1", "", "99"],
        "pos": ["100", "200", "300", "400", "500"],
        "rsid": ["rs1", "rs2", "rs3", "rs4", "rs5"],
    })
    g, ticks, labels = genome_layout(df, "chr", "pos")
    assert set(g["rsid"]) == {"rs1", "rs2", "rs3"}      # blank + foreign chr dropped
    assert labels == ["1", "2", "10"]                    # numeric order, not lexical
    assert ticks == sorted(ticks)                        # cumulative offsets increase
