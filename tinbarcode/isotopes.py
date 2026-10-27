"""Tin isotope reference data.

The isotope table supplied with the data is the source of truth. The NIST
values here are only a fallback for isotopes that table does not list
(for example the optional 117Sn and 119Sn) and a reference for tests.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Isotope:
    label: str  # e.g. "120Sn"
    mass: float  # exact mass in Da
    abundance: float  # natural abundance in %
    source: str = "NIST"


# NIST atomic weights and isotopic compositions.
NIST_SN = {
    "116Sn": Isotope("116Sn", 115.901744, 14.54),
    "117Sn": Isotope("117Sn", 116.902954, 7.68),
    "118Sn": Isotope("118Sn", 117.901607, 24.22),
    "119Sn": Isotope("119Sn", 118.903311, 8.59),
    "120Sn": Isotope("120Sn", 119.902199, 32.58),
}

ANCHOR = "120Sn"
CORE = ("118Sn", "116Sn")
OPTIONAL = ("119Sn", "117Sn")


def isotope_set(provided: dict[str, Isotope] | None, include_optional: bool) -> dict[str, Isotope]:
    """Return the isotopes to search for, preferring the provided table over NIST."""
    provided = provided or {}
    labels = (ANCHOR, *CORE, *(OPTIONAL if include_optional else ()))
    return {label: provided.get(label, NIST_SN[label]) for label in labels}


def expected_pattern(isotopes: dict[str, Isotope]) -> dict[str, tuple[float, float]]:
    """Offset from the anchor (Da) and intensity relative to the anchor (=100) per isotope."""
    anchor = isotopes[ANCHOR]
    return {
        label: (iso.mass - anchor.mass, 100.0 * iso.abundance / anchor.abundance)
        for label, iso in isotopes.items()
        if label != ANCHOR
    }
