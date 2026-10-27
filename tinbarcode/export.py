"""CSV export with the run parameters written in a commented header."""

from __future__ import annotations

import io
from datetime import datetime, timezone

import pandas as pd

from . import __version__
from .barcode import Settings
from .isotopes import Isotope


def hits_csv(
    hits: pd.DataFrame,
    settings: Settings,
    file_names: list[str],
    dataset: str,
    isotopes: dict[str, Isotope],
) -> str:
    lines = [
        f"# Sn barcode finder {__version__}",
        f"# Exported: {datetime.now(timezone.utc).isoformat(timespec='seconds')}",
        f"# Files: {', '.join(file_names)}",
        f"# Dataset searched: {dataset}",
    ]
    lines += [f"# {key} = {value}" for key, value in settings.to_dict().items()]
    lines += [f"# Isotope {iso.label}: mass {iso.mass} Da, abundance {iso.abundance} % ({iso.source})"
              for iso in isotopes.values()]
    buffer = io.StringIO()
    buffer.write("\n".join(lines) + "\n")
    hits.to_csv(buffer, index=False)
    return buffer.getvalue()


def points_csv(points: pd.DataFrame, title: str) -> str:
    buffer = io.StringIO()
    buffer.write(f"# Sn barcode finder {__version__}\n# Plot: {title}\n")
    points.to_csv(buffer, index=False)
    return buffer.getvalue()
