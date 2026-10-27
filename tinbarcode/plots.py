"""Plotly figures: full spectrum and a zoom per barcode hit."""

from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from .isotopes import ANCHOR, Isotope, expected_pattern, isotope_set

# Okabe-Ito colours, distinguishable with common colour-vision deficiencies.
SPECTRUM = "#4d4d4d"
CLASS_COLOURS = {"Strong": "#D55E00", "Possible": "#0072B2"}
EXPECTED = "#009E73"
OVERLAY = ["#4d4d4d", "#E69F00", "#56B4E9", "#CC79A7", "#F0E442", "#009E73"]

Y_LABELS = {
    "Raw": "Intensity (counts)",
    "Normalised": "Intensity (normalised to total ion count)",
}

# Show at least four decimals on the m/z axis once the visible range is under ~1 Da.
MZ_TICK_STOPS = [
    dict(dtickrange=[None, 0.25], value=".4f"),
    dict(dtickrange=[0.25, 5], value=".1f"),
    dict(dtickrange=[5, None], value=".0f"),
]

# Plotly's camera button: PNG at 300 dpi for an 8 x 4.5 inch figure.
PLOT_CONFIG = {
    "displaylogo": False,
    "toImageButtonOptions": {"format": "png", "width": 800, "height": 450, "scale": 3},
}


def y_label(dataset: str) -> str:
    return Y_LABELS.get(dataset, "Intensity")


def _sticks(mz: np.ndarray, intensity: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """One trace of vertical lines: (m, 0) -> (m, I), separated by gaps."""
    n = len(mz)
    x = np.empty(3 * n)
    y = np.empty(3 * n)
    x[0::3], x[1::3], x[2::3] = mz, mz, np.nan
    y[0::3], y[1::3], y[2::3] = 0.0, intensity, np.nan
    return x, y


def _layout(fig: go.Figure, title: str, dataset: str, x_range=None, y_range=None) -> go.Figure:
    fig.update_layout(
        title=dict(text=title, x=0, xref="paper", xanchor="left"),
        template="plotly_white",
        hovermode="closest",
        margin=dict(l=70, r=20, t=60, b=110),
        legend=dict(orientation="h", yanchor="top", y=-0.18, xanchor="left", x=0),
        font=dict(size=14),
    )
    fig.update_xaxes(title_text="m/z", tickformatstops=MZ_TICK_STOPS, showline=True, linecolor="black",
                     ticks="outside", range=x_range)
    fig.update_yaxes(title_text=y_label(dataset), showline=True, linecolor="black", ticks="outside",
                     exponentformat="power", range=y_range, rangemode="tozero")
    return fig


def _barcode_peaks(hits: pd.DataFrame, labels: list[str]) -> pd.DataFrame:
    """Long table of every highlighted peak: m/z, intensity, isotope, class, anchor."""
    rows = []
    for _, hit in hits.iterrows():
        rows.append((hit["Anchor m/z (120Sn)"], hit["120Sn intensity"], ANCHOR, hit["Class"], hit))
        for label in labels:
            if pd.notna(hit.get(f"{label} m/z")):
                rows.append((hit[f"{label} m/z"], hit[f"{label} intensity"], label, hit["Class"], hit))
    return pd.DataFrame(rows, columns=["mz", "intensity", "isotope", "class", "hit"])


def _first_formula(hit: pd.Series) -> str:
    """The best Sn formula MFP proposed for the anchor, if any."""
    text = hit.get("MFP candidates") or ""
    first = text.split(";")[0].split(" (")[0].strip()
    return first if "Sn" in first else ""


def full_spectrum(
    series: list[tuple[str, pd.DataFrame]],
    dataset: str,
    title: str,
    hits: pd.DataFrame | None = None,
    labels: list[str] | None = None,
    highlight: bool = True,
    x_range=None,
    y_range=None,
) -> go.Figure:
    """Stick spectrum of one or more samples, with barcode peaks highlighted."""
    fig = go.Figure()
    for k, (name, df) in enumerate(series):
        x, y = _sticks(df["mz"].to_numpy(float), df["intensity"].to_numpy(float))
        fig.add_trace(go.Scattergl(
            x=x, y=y, mode="lines", name=name, line=dict(color=OVERLAY[k % len(OVERLAY)], width=1),
            hovertemplate="m/z %{x:.4f}<br>intensity %{y:.4g}<extra>" + name + "</extra>",
        ))
    if highlight and hits is not None and not hits.empty:
        _add_highlights(fig, hits, labels or [], annotate_strong=True, possible_markers=False)
    return _layout(fig, title, dataset, x_range, y_range)


def _add_highlights(fig: go.Figure, hits: pd.DataFrame, labels: list[str], annotate_strong: bool,
                    possible_markers: bool = True) -> None:
    peaks = _barcode_peaks(hits, labels)
    for cls, colour in CLASS_COLOURS.items():
        sub = peaks[peaks["class"] == cls]
        if sub.empty:
            continue
        x, y = _sticks(sub["mz"].to_numpy(float), sub["intensity"].to_numpy(float))
        fig.add_trace(go.Scatter(x=x, y=y, mode="lines", name=f"{cls} barcode peaks",
                                 line=dict(color=colour, width=2.5), hoverinfo="skip"))
        if cls == "Possible" and not possible_markers:
            continue
        fig.add_trace(go.Scatter(
            x=sub["mz"], y=sub["intensity"], mode="markers", name=f"{cls} barcode peaks", showlegend=False,
            marker=dict(color=colour, size=8, symbol="diamond" if cls == "Strong" else "circle"),
            customdata=np.stack([sub["isotope"], sub["class"]], axis=1),
            hovertemplate="%{customdata[0]} (%{customdata[1]})<br>m/z %{x:.4f}<br>intensity %{y:.4g}<extra></extra>",
        ))
    if annotate_strong:
        for _, hit in hits[hits["Class"] == "Strong"].iterrows():
            formula = _first_formula(hit)
            text = f"{formula}<br>m/z {hit['Anchor m/z (120Sn)']:.4f}" if formula else f"m/z {hit['Anchor m/z (120Sn)']:.4f}"
            fig.add_annotation(x=hit["Anchor m/z (120Sn)"], y=hit["120Sn intensity"], text=text,
                               showarrow=True, arrowhead=0, ax=0, ay=-30, font=dict(size=11))


def hit_zoom(
    peaks: pd.DataFrame,
    hit: pd.Series,
    dataset: str,
    title: str,
    provided_isotopes: dict[str, Isotope] | None,
    include_optional: bool,
    below: float = 4.5,
    above: float = 1.5,
    highlight: bool = True,
    x_range=None,
    y_range=None,
) -> go.Figure:
    """Zoom on one hit with the expected barcode as dashed lines scaled to the 120Sn peak."""
    anchor = float(hit["Anchor m/z (120Sn)"])
    height = float(hit["120Sn intensity"])
    lo, hi = anchor - below, anchor + above
    window = peaks[(peaks["mz"] >= lo) & (peaks["mz"] <= hi)]
    fig = go.Figure()
    x, y = _sticks(window["mz"].to_numpy(float), window["intensity"].to_numpy(float))
    fig.add_trace(go.Scatter(x=x, y=y, mode="lines", name="Spectrum", line=dict(color=SPECTRUM, width=1.5),
                             hovertemplate="m/z %{x:.4f}<br>intensity %{y:.4g}<extra></extra>"))

    pattern = expected_pattern(isotope_set(provided_isotopes, include_optional))
    labels = list(pattern)
    if highlight:
        exp_x = [anchor] + [anchor + off for off, _ in pattern.values()]
        exp_y = [height] + [height * ratio / 100.0 for _, ratio in pattern.values()]
        names = [ANCHOR] + labels
        ex, ey = _sticks(np.array(exp_x), np.array(exp_y))
        fig.add_trace(go.Scatter(x=ex, y=ey, mode="lines", name="Expected barcode",
                                 line=dict(color=EXPECTED, width=1.5, dash="dash"), hoverinfo="skip"))
        fig.add_trace(go.Scatter(
            x=exp_x, y=exp_y, mode="markers", name="Expected barcode", showlegend=False, text=names,
            marker=dict(color=EXPECTED, size=7, symbol="line-ew-open", line=dict(width=2)),
            hovertemplate="expected %{text}<br>m/z %{x:.4f}<br>intensity %{y:.4g}<extra></extra>",
        ))
        _add_highlights(fig, hit.to_frame().T, labels, annotate_strong=False)
        formula = _first_formula(hit)
        for label in [ANCHOR] + labels:
            mz_col = "Anchor m/z (120Sn)" if label == ANCHOR else f"{label} m/z"
            int_col = "120Sn intensity" if label == ANCHOR else f"{label} intensity"
            if pd.notna(hit.get(mz_col)):
                text = f"{label}<br>{hit[mz_col]:.4f}"
                if label == ANCHOR and formula:
                    text = f"{formula}<br>" + text
                fig.add_annotation(x=hit[mz_col], y=hit[int_col], text=text, showarrow=False, yshift=22,
                                   font=dict(size=11))
    return _layout(fig, title, dataset, x_range or [lo, hi], y_range)


def plotted_points(fig: go.Figure) -> pd.DataFrame:
    """The data points behind a figure, for CSV export (stick gaps removed)."""
    frames = []
    for trace in fig.data:
        if trace.x is None or trace.y is None or trace.hoverinfo == "skip":
            continue
        x = np.asarray(trace.x, dtype=float)
        y = np.asarray(trace.y, dtype=float)
        keep = np.isfinite(x) & np.isfinite(y)
        if trace.mode == "lines":
            keep &= y != 0  # drop the baseline end of each stick
        frames.append(pd.DataFrame({"trace": trace.name or "", "mz": x[keep], "intensity": y[keep]}))
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=["trace", "mz", "intensity"])


def to_image(fig: go.Figure, fmt: str) -> bytes:
    """PNG at 300 dpi (8 x 4.5 inch) or SVG. Needs Kaleido and a Chrome/Chromium install."""
    return fig.to_image(format=fmt, width=800, height=450, scale=3 if fmt == "png" else 1)


def rescale_hits(hits: pd.DataFrame, peaks: pd.DataFrame) -> pd.DataFrame:
    """Express hit intensities in another dataset's units (e.g. raw instead of normalised).

    Each hit is scaled by the ratio of its anchor's intensity in ``peaks`` to
    the intensity it was found with, so ratios between its peaks are unchanged.
    """
    if hits.empty or peaks.empty:
        return hits
    mz = peaks["mz"].to_numpy(float)
    inten = peaks["intensity"].to_numpy(float)
    anchors = hits["Anchor m/z (120Sn)"].to_numpy(float)
    idx = np.clip(np.searchsorted(mz, anchors), 1, len(mz) - 1)
    idx = np.where(np.abs(mz[idx - 1] - anchors) <= np.abs(mz[idx] - anchors), idx - 1, idx)
    matched = np.abs(mz[idx] - anchors) <= 1e-4
    factor = np.where(matched, inten[idx] / hits["120Sn intensity"].to_numpy(float), np.nan)
    out = hits.copy()
    for col in [c for c in hits.columns if c.endswith("intensity")]:
        out[col] = hits[col].to_numpy(float) * factor
    return out[matched].reset_index(drop=True) if not matched.all() else out
