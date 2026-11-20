"""Legend names of fits and the j(z) figure."""
import matplotlib
matplotlib.use("Agg")

from arcwind.plotting._style import law_label


def test_law_label():
    assert law_label("broken-power-law") == "broken power law"
    assert law_label("broken-power-law-jz") == "broken power law"
    assert law_label("power-law-same-traces") == "power law"
    assert law_label("power-law-first") == "power law"
    assert "2.00" in law_label("broken-power-law_ymax2.00")
    assert law_label("my-fit") == "my fit"
    assert law_label("broken-power-law-constj") == "broken power law, constant j"
