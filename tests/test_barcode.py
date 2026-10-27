import time
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from tinbarcode import export, loader, plots
from tinbarcode.barcode import Settings, find_barcodes
from tinbarcode.isotopes import NIST_SN, expected_pattern, isotope_set
from tinbarcode.loader import hill_formula

EXAMPLE = Path(__file__).resolve().parent.parent / "Rory Trial.xlsx"


def barcode_peaks(anchor_mz, anchor_intensity=1000.0, shift=0.0, drop=(), scale=None):
    """Peaks for an ideal Sn barcode (NIST values) plus unrelated peaks."""
    pattern = expected_pattern(isotope_set({}, include_optional=False))
    rows = [(anchor_mz + shift, anchor_intensity)]
    for label, (offset, ratio) in pattern.items():
        if label not in drop:
            factor = (scale or {}).get(label, 1.0)
            rows.append((anchor_mz + offset + shift, anchor_intensity * ratio / 100 * factor))
    rows += [(50.0, 300.0), (60.5, 20.0), (200.0, 500.0)]
    return pd.DataFrame(rows, columns=["mz", "intensity"]).sort_values("mz", ignore_index=True)


def test_nist_spacings_match_brief():
    m = {k: v.mass for k, v in NIST_SN.items()}
    assert m["120Sn"] - m["118Sn"] == pytest.approx(2.00060, abs=1e-5)
    assert m["118Sn"] - m["116Sn"] == pytest.approx(1.99986, abs=1e-5)
    assert m["120Sn"] - m["116Sn"] == pytest.approx(4.00046, abs=1e-5)
    assert m["120Sn"] - m["119Sn"] == pytest.approx(0.99889, abs=1e-5)
    assert m["120Sn"] - m["117Sn"] == pytest.approx(2.99925, abs=1e-5)


def test_expected_ratios():
    pattern = expected_pattern(isotope_set({}, include_optional=True))
    assert pattern["118Sn"][1] == pytest.approx(74.3, abs=0.1)
    assert pattern["116Sn"][1] == pytest.approx(44.6, abs=0.1)


def test_ideal_barcode_is_strong():
    hits = find_barcodes(barcode_peaks(176.9), Settings())
    assert list(hits["Class"]) == ["Strong"]
    assert hits.iloc[0]["Anchor m/z (120Sn)"] == pytest.approx(176.9)
    assert hits.iloc[0]["Peaks found"] == "3 of 3"


def test_calibration_shift_does_not_matter():
    hits = find_barcodes(barcode_peaks(176.9, shift=-0.045), Settings(mass_tol=0.002))
    assert list(hits["Class"]) == ["Strong"]


def test_missing_peak_gives_possible():
    hits = find_barcodes(barcode_peaks(176.9, drop=("116Sn",)), Settings())
    assert list(hits["Class"]) == ["Possible"]
    assert "116Sn: no peak" in hits.iloc[0]["Notes"]


def test_wrong_ratio_is_rejected_and_explained():
    peaks = barcode_peaks(176.9, scale={"116Sn": 3.0})
    hits = find_barcodes(peaks, Settings())
    assert list(hits["Class"]) == ["Possible"]
    assert "ratio" in hits.iloc[0]["Notes"]
    assert find_barcodes(peaks, Settings(possible_missing=0)).empty


def test_ppm_tolerance():
    peaks = barcode_peaks(176.9)
    peaks.loc[peaks["mz"].sub(172.9).abs() < 0.01, "mz"] += 0.003  # 116Sn off by ~17 ppm
    assert find_barcodes(peaks, Settings(mass_tol=10, mass_tol_unit="ppm"))["Class"].tolist() == ["Possible"]
    assert find_barcodes(peaks, Settings(mass_tol=25, mass_tol_unit="ppm"))["Class"].tolist() == ["Strong"]


def test_min_intensity_ignores_small_peaks():
    hits = find_barcodes(barcode_peaks(176.9, anchor_intensity=100.0), Settings(min_intensity=50.0))
    assert list(hits["Class"]) == ["Possible"]
    assert "threshold" in hits.iloc[0]["Notes"]


def test_no_hits_returns_empty_table_with_columns():
    peaks = pd.DataFrame({"mz": [10.0, 20.0, 30.0], "intensity": [1.0, 2.0, 3.0]})
    hits = find_barcodes(peaks, Settings())
    assert hits.empty and "Class" in hits.columns


def test_mfp_cross_check():
    peaks = barcode_peaks(176.9)
    mfp = pd.DataFrame({"mz": [176.9, 176.9], "ppm": [-5.0, 1.0], "formula": ["C2H3NOSn", "C9H4O"],
                        "has_sn": [True, False], "Sn": [1, 0]})
    hit = find_barcodes(peaks, Settings(), mfp=mfp).iloc[0]
    assert hit["MFP check"] == "Sn formula"
    assert hit["MFP candidates"].startswith("C2H3NOSn")
    assert "C9H4O" in hit["MFP candidates"]  # every candidate is kept
    assert find_barcodes(peaks, Settings(), mfp=mfp.iloc[1:]).iloc[0]["MFP check"] == "No Sn formula"
    assert find_barcodes(peaks, Settings(), mfp=mfp.iloc[0:0]).iloc[0]["MFP check"] == "No MFP candidates"


def test_hill_formula():
    assert hill_formula({"H": 3, "C": 2, "N": 1, "O": 1, "Sn": 1, "F": 0}) == "C2H3NOSn"
    assert hill_formula({"H": 1, "Na": 1, "F": 1}) == "FHNa"


def test_fifteen_thousand_peaks_under_a_second():
    rng = np.random.default_rng(0)
    peaks = pd.DataFrame({"mz": np.sort(rng.uniform(10, 1500, 15000)), "intensity": rng.lognormal(5, 2, 15000)})
    start = time.perf_counter()
    find_barcodes(peaks, Settings())
    assert time.perf_counter() - start < 1.0


def test_csv_export_header():
    hits = find_barcodes(barcode_peaks(176.9), Settings())
    text = export.hits_csv(hits, Settings(), ["a.xlsx"], "Raw", isotope_set({}, False))
    assert text.startswith("# Sn barcode finder")
    assert "# mass_tol = 0.008" in text and "a.xlsx" in text
    assert "Anchor m/z (120Sn)" in text


def test_plots_build_and_export_points():
    peaks = barcode_peaks(176.9)
    hits = find_barcodes(peaks, Settings())
    fig = plots.full_spectrum([("x", peaks)], "Raw", "t", hits, ["118Sn", "116Sn"])
    assert len(plots.plotted_points(fig)) >= len(peaks)
    zoom = plots.hit_zoom(peaks, hits.iloc[0], "Raw", "t", {}, False)
    assert list(zoom.layout.xaxis.range) == pytest.approx([172.4, 178.4])


def test_rescale_hits_between_datasets():
    peaks = barcode_peaks(176.9)
    hits = find_barcodes(peaks, Settings())
    normalised = peaks.assign(intensity=peaks["intensity"] / 1000)
    scaled = plots.rescale_hits(hits, normalised)
    assert scaled.iloc[0]["120Sn intensity"] == pytest.approx(1.0)
    assert scaled.iloc[0]["118Sn ratio (%)"] == hits.iloc[0]["118Sn ratio (%)"]


def test_plain_csv_peak_list(tmp_path):
    path = tmp_path / "sample_(-).csv"
    barcode_peaks(176.9).rename(columns={"mz": "m/z"}).to_csv(path, index=False)
    ds = loader.load(path)
    assert ds.polarity == "-"
    assert list(ds.peaks) == ["Peaks"] and len(ds.peaks["Peaks"]) == 6


@pytest.fixture(scope="module")
def ds():
    return loader.load(EXAMPLE)


@pytest.mark.skipif(not EXAMPLE.exists(), reason="example workbook not present")
class TestExampleWorkbook:
    def test_loads_every_table(self, ds):
        assert ds.polarity == "+"
        assert ds.sample == "NAG+SnF2_1-100_LMIG_ToF_(+)1"
        assert {k: len(v) for k, v in ds.peaks.items()} == {"Normalised": 302, "Raw": 302}
        assert len(ds.mfp) == 7889 and ds.mfp["has_sn"].sum() == 913
        assert ds.isotopes["120Sn"].abundance == 32.4
        assert ds.isotopes["118Sn"].mass == pytest.approx(117.901607)

    def test_finds_the_sn_ion(self, ds):
        hits = find_barcodes(ds.peaks["Normalised"], Settings(), ds.isotopes, ds.mfp)
        strong = hits[hits["Class"] == "Strong"]["Anchor m/z (120Sn)"].round(3).tolist()
        assert 119.857 in strong  # Sn+
        assert 136.851 in strong  # SnOH+

    def test_raw_and_normalised_agree(self, ds):
        norm = find_barcodes(ds.peaks["Normalised"], Settings(), ds.isotopes)
        raw = find_barcodes(ds.peaks["Raw"], Settings(), ds.isotopes)
        assert norm["Class"].value_counts().get("Strong") == raw["Class"].value_counts().get("Strong")
