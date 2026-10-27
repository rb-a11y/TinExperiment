"""Remember the user's last settings between sessions in a small JSON file."""

from __future__ import annotations

import json
import os
from pathlib import Path

PREFS_PATH = Path(os.environ.get("TINBARCODE_PREFS", Path.home() / ".tinbarcode" / "settings.json"))


def load() -> dict:
    try:
        return json.loads(PREFS_PATH.read_text())
    except (OSError, ValueError):
        return {}


def save(prefs: dict) -> None:
    try:
        PREFS_PATH.parent.mkdir(parents=True, exist_ok=True)
        PREFS_PATH.write_text(json.dumps(prefs, indent=2, sort_keys=True))
    except OSError:
        pass  # Remembering settings is a convenience; never fail the app over it.
