"""Detect the Sn isotope barcode in a peak list.

Every peak above the intensity threshold is tried as the 120Sn peak (the
anchor). For each other isotope the expected position is the anchor's
*observed* m/z plus the isotope's mass offset, so a calibration shift that
moves all peaks together does not matter. A partner peak counts as found
when it lies within the mass tolerance of that position and its intensity
relative to the anchor is within the ratio tolerance of the expected value.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np
import pandas as pd

from .isotopes import ANCHOR, Isotope, expected_pattern, isotope_set

STRONG, POSSIBLE = "Strong", "Possible"


@dataclass
class Settings:
    """All matching rules. Defaults are proposals for the scientist to confirm."""

    mass_tol: float = 0.008  # allowed error for each expected isotope peak
    mass_tol_unit: str = "Da"  # "Da" or "ppm"
    ratio_tol: float = 25.0  # allowed drift, % of the expected relative intensity
    min_intensity: float = 0.0  # peaks below this are ignored
    include_optional: bool = False  # also test 117Sn and 119Sn
    strong_missing: int = 0  # peaks that may be missing for Strong (0 = all found)
    possible_missing: int = 1  # peaks that may be missing for Possible
    near_threshold_factor: float = 2.0  # flag peaks below this multiple of min_intensity
    mfp_match_tol: float = 0.0005  # Da, to link a peak to its MFP rows
    hide_shadows: bool = True  # drop Possible hits anchored on a peak of a Strong hit

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict) -> "Settings":
        known = cls.__dataclass_fields__
        return cls(**{k: v for k, v in data.items() if k in known})


def tolerance_da(mz: np.ndarray, settings: Settings) -> np.ndarray:
    if settings.mass_tol_unit == "ppm":
        return mz * settings.mass_tol * 1e-6
    return np.full_like(mz, settings.mass_tol, dtype=float)


def _nearest(sorted_mz: np.ndarray, targets: np.ndarray) -> np.ndarray:
    """Index of the nearest value in sorted_mz for each target."""
    idx = np.clip(np.searchsorted(sorted_mz, targets), 1, len(sorted_mz) - 1)
    left, right = sorted_mz[idx - 1], sorted_mz[idx]
    return np.where(np.abs(targets - left) <= np.abs(right - targets), idx - 1, idx)


def find_barcodes(
    peaks: pd.DataFrame,
    settings: Settings,
    provided_isotopes: dict[str, Isotope] | None = None,
    mfp: pd.DataFrame | None = None,
    sample: str = "",
    polarity: str = "?",
) -> pd.DataFrame:
    """Return one row per anchor peak classed Strong or Possible."""
    isotopes = isotope_set(provided_isotopes, settings.include_optional)
    pattern = expected_pattern(isotopes)
    labels = list(pattern)
    n_tested = 1 + len(labels)

    kept = peaks[peaks["intensity"] >= settings.min_intensity].sort_values("mz")
    mz = kept["mz"].to_numpy(float)
    inten = kept["intensity"].to_numpy(float)
    if len(mz) < 2:
        return _empty(labels)

    found = np.ones(len(mz), dtype=int)  # the anchor itself
    per_iso = {}
    for label, (offset, ratio) in pattern.items():
        target = mz + offset
        j = _nearest(mz, target)
        mass_ok = np.abs(mz[j] - target) <= tolerance_da(target, settings)
        observed = 100.0 * inten[j] / inten
        ratio_ok = np.abs(observed - ratio) <= settings.ratio_tol / 100.0 * ratio
        ok = mass_ok & ratio_ok
        found += ok
        per_iso[label] = dict(j=j, target=target, mass_ok=mass_ok, ok=ok, observed=observed, expected=ratio)

    missing = n_tested - found
    strong = missing <= settings.strong_missing
    possible = missing <= settings.possible_missing
    if settings.hide_shadows:
        # A real barcode read one isotope too low (118Sn taken as 120Sn) looks like a Possible hit.
        members = np.concatenate([info["j"][strong & info["ok"]] for info in per_iso.values()])
        possible &= ~np.isin(np.arange(len(mz)), members)
    hits = np.flatnonzero(strong | possible)
    if len(hits) == 0:
        return _empty(labels)

    rows = []
    near = settings.near_threshold_factor * settings.min_intensity
    for i in hits:
        row = {
            "Class": STRONG if strong[i] else POSSIBLE,
            "Anchor m/z (120Sn)": mz[i],
            "120Sn intensity": inten[i],
            "Peaks found": f"{found[i]} of {n_tested}",
        }
        notes = []
        low_peaks = [ANCHOR] if settings.min_intensity > 0 and inten[i] < near else []
        for label in labels:
            info = per_iso[label]
            j = info["j"][i]
            if info["ok"][i]:
                row[f"{label} m/z"] = mz[j]
                row[f"{label} error (ppm)"] = (mz[j] - info["target"][i]) / info["target"][i] * 1e6
                row[f"{label} intensity"] = inten[j]
                row[f"{label} ratio (%)"] = info["observed"][i]
                if settings.min_intensity > 0 and inten[j] < near:
                    low_peaks.append(label)
            else:
                row[f"{label} m/z"] = row[f"{label} error (ppm)"] = np.nan
                row[f"{label} intensity"] = row[f"{label} ratio (%)"] = np.nan
                expected_intensity = inten[i] * info["expected"] / 100.0
                if info["mass_ok"][i]:
                    notes.append(
                        f"{label}: peak at m/z {mz[j]:.4f} but ratio {info['observed'][i]:.0f}% "
                        f"(expected {info['expected']:.0f}%)"
                    )
                elif settings.min_intensity > 0 and expected_intensity < near:
                    notes.append(f"{label}: expected intensity is near or below the threshold")
                else:
                    notes.append(f"{label}: no peak at expected m/z {info['target'][i]:.4f}")
        if low_peaks:
            notes.append("Near intensity threshold: " + ", ".join(low_peaks))
        row.update(_mfp_check(mfp, mz[i], settings))
        row["Notes"] = "; ".join(notes)
        row["Polarity"] = polarity
        row["Sample"] = sample
        rows.append(row)
    return pd.DataFrame(rows, columns=_columns(labels))


def _mfp_check(mfp: pd.DataFrame | None, anchor_mz: float, settings: Settings) -> dict:
    if mfp is None:
        return {"MFP check": "No MFP data", "MFP candidates": ""}
    cands = mfp[np.abs(mfp["mz"] - anchor_mz) <= settings.mfp_match_tol]
    if cands.empty:
        return {"MFP check": "No MFP candidates", "MFP candidates": ""}
    cands = cands.assign(_abs=cands["ppm"].abs()).sort_values(["has_sn", "_abs"], ascending=[False, True])
    text = "; ".join(f"{f} ({p:+.0f} ppm)" for f, p in zip(cands["formula"], cands["ppm"]))
    sn = cands[cands["has_sn"]]
    if sn.empty:
        flag = "No Sn formula"
    elif (sn["Sn"] > 1).all():
        flag = "Sn formula (Sn2+ only)"
    else:
        flag = "Sn formula"
    return {"MFP check": flag, "MFP candidates": text}


def _columns(labels: list[str]) -> list[str]:
    cols = ["Class", "Anchor m/z (120Sn)", "120Sn intensity", "Peaks found"]
    for label in labels:
        cols += [f"{label} m/z", f"{label} error (ppm)", f"{label} intensity", f"{label} ratio (%)"]
    return cols + ["MFP check", "MFP candidates", "Notes", "Polarity", "Sample"]


def _empty(labels: list[str]) -> pd.DataFrame:
    return pd.DataFrame(columns=_columns(labels))
