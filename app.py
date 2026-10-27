"""Sn barcode finder: a web interface for finding tin isotope patterns in SIMS peak lists.

Run with:  streamlit run app.py
"""

from __future__ import annotations

import time
from pathlib import Path

import pandas as pd
import streamlit as st

from tinbarcode import __version__, export, loader, plots, prefs
from tinbarcode.barcode import Settings, find_barcodes
from tinbarcode.isotopes import ANCHOR, NIST_SN, expected_pattern, isotope_set

EXAMPLE_FILE = Path(__file__).with_name("Rory Trial.xlsx")
DATASET_ORDER = ["Normalised", "Raw"]
ROW_COLOURS = {"Strong": "#fbe0cf", "Possible": "#d5e8f6"}
DEFAULT_TOL = {"Da": 0.008, "ppm": 40.0}

# Every widget's default. The user's last values are restored from the prefs file.
DEFAULTS = {
    "use_example": False,
    "search_dataset": "Normalised",
    "mass_tol_unit": "Da",
    "mass_tol": DEFAULT_TOL["Da"],
    "ratio_tol": 25.0,
    "min_intensity": 0.0,
    "include_optional": False,
    "strong_missing": 0,
    "possible_missing": 1,
    "hide_shadows": True,
    "highlight": True,
    "class_filter": ["Strong", "Possible"],
    "show_plots": False,
    "plot_dataset": "Normalised",
    "overlay": False,
    "lock_axes": False,
    "x_min": None,
    "x_max": None,
    "y_min": None,
    "y_max": None,
    "zoom_below": 4.5,
    "zoom_above": 1.5,
    "zoom_include_possible": False,
}

st.set_page_config(page_title="Sn barcode finder", layout="wide")

MIN_STREAMLIT = (1, 40)
if tuple(int(part) for part in st.__version__.split(".")[:2]) < MIN_STREAMLIT:
    st.error(f"This app needs Streamlit 1.40 or newer, but this computer has {st.__version__}. "
             "Run `pip install --upgrade streamlit` in Terminal, then start the app again.")
    st.stop()


# ------------------------------------------------------------------ state


def init_state() -> None:
    if st.session_state.get("_initialised"):
        # Keep values of widgets that are hidden on this run (Streamlit drops their state otherwise).
        for key in DEFAULTS:
            if key in st.session_state:
                st.session_state[key] = st.session_state[key]
        return
    saved = prefs.load()
    for key, default in DEFAULTS.items():
        st.session_state[key] = saved.get(key, default)
    st.session_state["_initialised"] = True


def save_state() -> None:
    current = {key: st.session_state.get(key, default) for key, default in DEFAULTS.items()}
    if current != st.session_state.get("_saved"):
        prefs.save(current)
        st.session_state["_saved"] = current


def reset_settings() -> None:
    for key, default in DEFAULTS.items():
        st.session_state[key] = default


def reset_axes() -> None:
    for key in ("x_min", "x_max", "y_min", "y_max"):
        st.session_state[key] = None


def unit_changed() -> None:
    st.session_state["mass_tol"] = DEFAULT_TOL[st.session_state["mass_tol_unit"]]


@st.cache_data(show_spinner="Reading file…")
def load_file(data: bytes, name: str) -> loader.Dataset:
    return loader.load(data, name)


def pick_peaks(ds: loader.Dataset, preferred: str) -> tuple[str, pd.DataFrame]:
    name = preferred if preferred in ds.peaks else next(iter(ds.peaks))
    return name, ds.peaks[name]


# ------------------------------------------------------------------ sidebar


def sidebar() -> list[loader.Dataset]:
    sb = st.sidebar
    sb.header("1. Load data")
    uploads = sb.file_uploader(
        "Peak list or MFP workbook (.xlsx, .csv, .txt)",
        type=["xlsx", "xlsm", "csv", "txt", "tsv"],
        accept_multiple_files=True,
        help="Upload one or more files. Each file is searched separately and listed as its own sample.",
    )
    if EXAMPLE_FILE.exists():
        sb.checkbox(f"Use the example file ({EXAMPLE_FILE.name})", key="use_example")

    datasets = []
    names = set()
    if st.session_state["use_example"] and EXAMPLE_FILE.exists():
        datasets.append(load_file(EXAMPLE_FILE.read_bytes(), EXAMPLE_FILE.name))
        names.add(EXAMPLE_FILE.name)
    for up in uploads or []:
        if up.name in names:
            continue
        try:
            datasets.append(load_file(up.getvalue(), up.name))
            names.add(up.name)
        except Exception as exc:  # show the problem, keep the other files
            sb.error(f"Could not read {up.name}: {exc}")

    sb.header("2. Matching rules")
    sb.radio("Search the", DATASET_ORDER, key="search_dataset", horizontal=True,
             help="Which intensities to search. Ratios are the same in both; only the threshold units differ.")
    sb.radio("Mass tolerance unit", ["Da", "ppm"], key="mass_tol_unit", on_change=unit_changed, horizontal=True)
    sb.number_input(
        f"Mass tolerance ({st.session_state['mass_tol_unit']})",
        min_value=0.0, step=0.001 if st.session_state["mass_tol_unit"] == "Da" else 1.0,
        format="%.4f" if st.session_state["mass_tol_unit"] == "Da" else "%.1f", key="mass_tol",
        help="How far each isotope peak may sit from where the barcode predicts it, measured from the "
             "observed 120Sn peak. This ToF data needs about 0.008 Da (roughly 40 ppm at m/z 200).",
    )
    sb.number_input("Ratio tolerance (% of expected)", min_value=0.0, max_value=200.0, step=5.0, key="ratio_tol",
                    help="Allowed drift in the 116 : 118 : 120 intensity ratios. 25 means 118Sn may be 75 % ± 19 %.")
    sb.number_input("Minimum intensity", min_value=0.0, format="%g", key="min_intensity",
                    help="Peaks below this are ignored. Units follow the dataset searched.")
    sb.checkbox("Also test 117Sn and 119Sn", key="include_optional",
                help="Adds the two smaller isotopes. Their values come from NIST if the file's isotope table lacks them.")
    n_tested = 5 if st.session_state["include_optional"] else 3
    sb.number_input("Strong: peaks allowed missing", min_value=0, max_value=n_tested - 1, step=1, key="strong_missing",
                    help=f"0 means all {n_tested} peaks must be found.")
    sb.number_input("Possible: peaks allowed missing", min_value=0, max_value=n_tested - 1, step=1,
                    key="possible_missing", help=f"1 means {n_tested - 1} of {n_tested} peaks.")
    sb.checkbox("Hide echoes of Strong hits", key="hide_shadows",
                help="A real tin pattern also looks like a weaker pattern one isotope lower (118Sn read as 120Sn). "
                     "This hides Possible hits whose 120Sn peak already belongs to a Strong hit.")

    sb.header("3. Display")
    sb.toggle("Highlight barcodes", key="highlight", help="One switch for all highlighting on the plots.")
    sb.button("Reset all settings", on_click=reset_settings)
    sb.caption(f"Sn barcode finder {__version__}. Settings are remembered between sessions.")
    return datasets


def current_settings() -> Settings:
    s = st.session_state
    return Settings(
        mass_tol=float(s["mass_tol"]), mass_tol_unit=s["mass_tol_unit"], ratio_tol=float(s["ratio_tol"]),
        min_intensity=float(s["min_intensity"]), include_optional=bool(s["include_optional"]),
        strong_missing=int(s["strong_missing"]), possible_missing=int(s["possible_missing"]),
        hide_shadows=bool(s["hide_shadows"]),
    )


# ------------------------------------------------------------------ main panels


def file_notes(datasets: list[loader.Dataset]) -> None:
    polarities = {ds.polarity for ds in datasets} - {"?"}
    if len(polarities) > 1:
        st.error("Positive and negative polarity files are loaded together. They are searched separately, "
                 "but never compare or overlay them.")
    for ds in datasets:
        for warning in ds.warnings:
            st.warning(f"{ds.file_name}: {warning}")


def isotope_panel(datasets: list[loader.Dataset], settings: Settings) -> None:
    with st.expander("Isotope values used"):
        for ds in datasets:
            isotopes = isotope_set(ds.isotopes, settings.include_optional)
            pattern = expected_pattern(isotopes)
            rows = []
            for label, iso in isotopes.items():
                offset, ratio = pattern.get(label, (0.0, 100.0))
                nist = NIST_SN[label]
                rows.append({
                    "Isotope": label, "Mass (Da)": iso.mass, "Abundance (%)": iso.abundance, "Source": iso.source,
                    "Offset from 120Sn (Da)": offset, "Expected intensity (120Sn = 100)": ratio,
                    "NIST mass (Da)": nist.mass, "NIST abundance (%)": nist.abundance,
                })
            st.markdown(f"**{ds.file_name}**")
            table = pd.DataFrame(rows)
            st.dataframe(table, hide_index=True, column_config={
                "Mass (Da)": st.column_config.NumberColumn(format="%.6f"),
                "NIST mass (Da)": st.column_config.NumberColumn(format="%.6f"),
                "Offset from 120Sn (Da)": st.column_config.NumberColumn(format="%.5f"),
                "Expected intensity (120Sn = 100)": st.column_config.NumberColumn(format="%.1f"),
            })
            differs = [r["Isotope"] for r in rows if r["Source"] != "NIST"
                       and abs(r["Abundance (%)"] - r["NIST abundance (%)"]) > 0.05]
            if differs:
                st.caption("The file's abundances differ from NIST for " + ", ".join(differs)
                           + ". The file's values are used; please confirm them with the scientist.")


def run_search(datasets: list[loader.Dataset], settings: Settings) -> pd.DataFrame:
    frames = []
    summary = []
    for ds in datasets:
        name, peaks = pick_peaks(ds, st.session_state["search_dataset"])
        start = time.perf_counter()
        hits = find_barcodes(peaks, settings, ds.isotopes, ds.mfp, ds.sample, ds.polarity)
        elapsed = (time.perf_counter() - start) * 1000
        hits.insert(0, "File", ds.file_name)
        frames.append(hits)
        counts = hits["Class"].value_counts()
        if hits.empty:
            summary.append(f"**{ds.file_name}**: no tin barcode found with these settings "
                           f"({len(peaks)} {name.lower()} peaks searched). Try a wider mass or ratio tolerance, "
                           "or a lower minimum intensity.")
        else:
            summary.append(f"**{ds.file_name}**: {counts.get('Strong', 0)} Strong and "
                           f"{counts.get('Possible', 0)} Possible hits in {len(peaks)} {name.lower()} peaks "
                           f"(searched in {elapsed:.0f} ms).")
    for line in summary:
        st.markdown(line)
    frames = [f for f in frames if not f.empty]
    if not frames:
        return pd.DataFrame()
    hits = pd.concat(frames, ignore_index=True)
    order = hits["Class"].map({"Strong": 0, "Possible": 1})
    return hits.assign(_order=order).sort_values(["_order", "File", "Anchor m/z (120Sn)"]).drop(columns="_order") \
        .reset_index(drop=True)


def results_table(hits: pd.DataFrame) -> pd.DataFrame:
    """Show the filterable hit table. Returns the selected hit (0 or 1 rows)."""
    f1, f2 = st.columns(2)
    f1.multiselect("Show classes", ["Strong", "Possible"], key="class_filter")
    flags = sorted(hits["MFP check"].unique())
    st.session_state.setdefault("mfp_filter", flags)
    st.session_state["mfp_filter"] = [f for f in st.session_state["mfp_filter"] if f in flags] or flags
    f2.multiselect("Show MFP check results", flags, key="mfp_filter",
                   help="'Sn formula' means MFP proposed at least one tin-containing formula for the 120Sn peak.")
    shown = hits[hits["Class"].isin(st.session_state["class_filter"])
                 & hits["MFP check"].isin(st.session_state["mfp_filter"])].reset_index(drop=True)
    if shown.empty:
        st.info("No hits match the filters.")
        return shown

    def colour(row):
        return [f"background-color: {ROW_COLOURS.get(row['Class'], 'white')}"] * len(row)

    number_formats = {c: "{:.4f}" for c in shown.columns if c.endswith("m/z") or c.endswith("(120Sn)")}
    number_formats |= {c: "{:+.1f}" for c in shown.columns if c.endswith("(ppm)")}
    number_formats |= {c: "{:.0f}" for c in shown.columns if c.endswith("(%)")}
    number_formats |= {c: "{:.4g}" for c in shown.columns if c.endswith("intensity")}
    styled = shown.style.apply(colour, axis=1).format(number_formats, na_rep="–")
    st.caption("Select a row to zoom on that hit. Strong rows are orange, Possible rows are blue; "
               "the Class column says the same in words.")
    event = st.dataframe(styled, hide_index=True, on_select="rerun", selection_mode="single-row",
                         key="hit_table", height=min(38 + 35 * len(shown), 420))
    selected = event.selection.rows if event and event.selection else []
    return shown.iloc[selected]


def download_figure(fig, title: str, key: str) -> None:
    stem = "".join(ch if ch.isalnum() or ch in "-_." else "_" for ch in title)[:80]
    cols = st.columns(3)
    # Images take a second or two to render, so they are only made when asked for.
    signature = hash(fig.to_json())
    images = st.session_state.get(f"img_{key}")
    if images and images["signature"] == signature:
        cols[0].download_button("Download PNG (300 dpi)", images["png"], f"{stem}.png", "image/png",
                                key=f"png_{key}")
        cols[1].download_button("Download SVG", images["svg"], f"{stem}.svg", "image/svg+xml", key=f"svg_{key}")
    elif cols[0].button("Prepare PNG and SVG", key=f"make_{key}"):
        try:
            with st.spinner("Rendering images…"):
                st.session_state[f"img_{key}"] = {
                    "signature": signature, "png": plots.to_image(fig, "png"), "svg": plots.to_image(fig, "svg"),
                }
            st.rerun()
        except Exception:
            cols[1].warning("PNG/SVG export needs Google Chrome installed (or run `plotly_get_chrome` once). "
                            "The camera icon on the plot saves a PNG without it.")
    cols[2].download_button("Plotted points (CSV)", export.points_csv(plots.plotted_points(fig), title),
                            f"{stem}.csv", "text/csv", key=f"csv_{key}")


def axis_limits(for_zoom: bool):
    s = st.session_state
    if for_zoom and not s["lock_axes"]:
        return None, None
    x = [s["x_min"], s["x_max"]] if s["x_min"] is not None and s["x_max"] is not None else None
    y = [s["y_min"], s["y_max"]] if s["y_min"] is not None and s["y_max"] is not None else None
    return x, y


def zoom_figure(ds: loader.Dataset, hit: pd.Series, dataset: str, settings: Settings):
    name, peaks = pick_peaks(ds, dataset)
    scaled = plots.rescale_hits(hit.to_frame().T.reset_index(drop=True), peaks)
    if scaled.empty:
        return None, ""
    hit = scaled.iloc[0]
    title = f"120Sn at m/z {hit['Anchor m/z (120Sn)']:.4f} | {ds.file_name} | {ds.polarity} | {name}"
    x, y = axis_limits(for_zoom=True)
    plot_title = f"120Sn at m/z {hit['Anchor m/z (120Sn)']:.4f}<br><sup>{ds.file_name} | {ds.polarity} | {name}</sup>"
    fig = plots.hit_zoom(peaks, hit, name, plot_title, ds.isotopes, settings.include_optional,
                         below=st.session_state["zoom_below"], above=st.session_state["zoom_above"],
                         highlight=st.session_state["highlight"], x_range=x, y_range=y)
    return fig, title


def plot_section(datasets: list[loader.Dataset], hits: pd.DataFrame, settings: Settings) -> None:
    st.toggle("Show spectrum plots", key="show_plots",
              help="Off by default. Draws the full spectrum and a zoom for each hit.")
    if not st.session_state["show_plots"]:
        return
    s = st.session_state
    available = [d for d in DATASET_ORDER if any(d in ds.peaks for ds in datasets)] or \
        sorted({k for ds in datasets for k in ds.peaks})
    if s["plot_dataset"] not in available:
        s["plot_dataset"] = available[0]

    c1, c2, c3 = st.columns(3)
    c1.radio("Plot the", available, key="plot_dataset", horizontal=True)
    if len(datasets) > 1:
        c2.checkbox("Overlay all samples", key="overlay")
    c3.checkbox("Lock axis limits for every plot", key="lock_axes",
                help="Applies the limits below to the zoom plots too, for fair comparison.")
    a1, a2, a3, a4, a5 = st.columns([1, 1, 1, 1, 0.8], vertical_alignment="bottom")
    a1.number_input("x minimum (m/z)", key="x_min", format="%.4f", placeholder="auto")
    a2.number_input("x maximum (m/z)", key="x_max", format="%.4f", placeholder="auto")
    a3.number_input("y minimum", key="y_min", format="%g", placeholder="auto")
    a4.number_input("y maximum", key="y_max", format="%g", placeholder="auto")
    a5.button("Reset axes", on_click=reset_axes)
    st.caption("Drag on a plot to zoom; double-click it to zoom back out.")

    dataset = s["plot_dataset"]
    x, y = axis_limits(for_zoom=False)
    labels = list(expected_pattern(isotope_set({}, settings.include_optional)))
    polarities = {ds.polarity for ds in datasets} - {"?"}
    groups = [datasets] if s["overlay"] and len(datasets) > 1 and len(polarities) <= 1 else [[ds] for ds in datasets]
    if s["overlay"] and len(polarities) > 1:
        st.warning("Overlay is off because the files have different polarities.")

    st.subheader("Full spectrum")
    for k, group in enumerate(groups):
        series, group_hits = [], []
        for ds in group:
            name, peaks = pick_peaks(ds, dataset)
            series.append((ds.sample if len(group) > 1 else name, peaks))
            if not hits.empty:
                group_hits.append(plots.rescale_hits(hits[hits["File"] == ds.file_name].reset_index(drop=True), peaks))
        ds_hits = pd.concat(group_hits, ignore_index=True) if group_hits else pd.DataFrame()
        files = ", ".join(ds.file_name for ds in group)
        title = f"{files} | {group[0].polarity} | {dataset}"
        fig = plots.full_spectrum(series, dataset, title, ds_hits, labels, s["highlight"], x, y)
        st.plotly_chart(fig, config=plots.PLOT_CONFIG, theme=None, key=f"full_{k}")
        download_figure(fig, title, f"full_{k}")

    if hits.empty:
        return
    st.subheader("Zoom on each hit")
    z1, z2, z3 = st.columns(3)
    z1.number_input("Window below 120Sn (Da)", min_value=0.5, step=0.5, key="zoom_below")
    z2.number_input("Window above 120Sn (Da)", min_value=0.5, step=0.5, key="zoom_above")
    z3.checkbox("Include Possible hits", key="zoom_include_possible")
    shown = hits if s["zoom_include_possible"] else hits[hits["Class"] == "Strong"]
    if shown.empty:
        st.info("No Strong hits to zoom on. Tick 'Include Possible hits' to see the others.")
        return
    by_name = {ds.file_name: ds for ds in datasets}
    grid = st.columns(2)
    for k, (_, hit) in enumerate(shown.iterrows()):
        fig, title = zoom_figure(by_name[hit["File"]], hit, dataset, settings)
        if fig is None:
            continue
        with grid[k % 2]:
            st.plotly_chart(fig, config=plots.PLOT_CONFIG, theme=None, key=f"zoom_{k}")
            download_figure(fig, title, f"zoom_{k}")


def main() -> None:
    init_state()
    datasets = sidebar()
    settings = current_settings()

    st.title("Sn barcode finder")
    st.markdown(
        "Finds ions that carry the tin isotope pattern: peaks for **116Sn, 118Sn and 120Sn** about 2 Da apart, "
        "with heights near **45 : 75 : 100**. Every result is a putative assignment to check by eye."
    )
    if not datasets:
        st.info("Load a file in the sidebar to start" + (", or tick 'Use the example file'." if EXAMPLE_FILE.exists() else "."))
        save_state()
        return

    file_notes(datasets)
    isotope_panel(datasets, settings)

    st.header("Possible tin ions")
    hits = run_search(datasets, settings)
    if not hits.empty:
        selected = results_table(hits)
        file_names = [ds.file_name for ds in datasets]
        isotopes = isotope_set(datasets[0].isotopes, settings.include_optional)
        st.download_button(
            "Download hit table (CSV)",
            export.hits_csv(hits, settings, file_names, st.session_state["search_dataset"], isotopes),
            "sn_barcode_hits.csv", "text/csv",
        )
        if not selected.empty:
            hit = selected.iloc[0]
            by_name = {ds.file_name: ds for ds in datasets}
            dataset = st.session_state["plot_dataset"] if st.session_state["show_plots"] else st.session_state["search_dataset"]
            fig, title = zoom_figure(by_name[hit["File"]], hit, dataset, settings)
            if fig is not None:
                st.subheader("Selected hit")
                st.plotly_chart(fig, config=plots.PLOT_CONFIG, theme=None, key="selected_zoom")
                download_figure(fig, title, "selected")

    st.header("Spectra")
    plot_section(datasets, hits, settings)
    save_state()


main()
