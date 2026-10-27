# Sn barcode finder

A small app that finds tin (Sn) isotope patterns in SIMS peak lists and plots
the spectra. It is the proof of concept described in
`Sn Barcode Proof of Concept Developer Brief.docx`.

Tin comes in several isotopes, so an ion that contains tin shows a group of
peaks: 116Sn, 118Sn and 120Sn, about 2 Da apart, with heights near
45 : 75 : 100. That group is the "barcode". The app looks for it in every
peak list you load, scores each match, checks it against the MFP formulae and
highlights it in a table and on the plots.

## Run it

You need Python 3.10 or newer and Streamlit 1.40 or newer (an older Streamlit,
for example from Anaconda, shows a message asking you to upgrade).

```bash
pip install -r requirements.txt
streamlit run app.py
```

Your browser opens the app. Tick **Use the example file** in the sidebar to
load `Rory Trial.xlsx`, or upload your own files.

PNG and SVG export uses Kaleido, which needs Google Chrome or Chromium. If
neither is installed, run `plotly_get_chrome` once, or use the camera icon on
any plot, which saves a 300 dpi PNG.

## Using the app

1. **Load data.** Upload one or more `.xlsx` workbooks laid out like
   `Rory Trial.xlsx`, or plain peak lists (`.csv`/`.txt` with an m/z column
   and an intensity column). Each file is searched separately.
2. **Matching rules** (sidebar). Every rule is a setting:
   | Setting | Default | Meaning |
   |---|---|---|
   | Search the | Normalised | Which intensities to search. Ratios are the same in both. |
   | Mass tolerance | 0.008 Da | How far each isotope peak may sit from its expected position. Can be set in Da or ppm. |
   | Ratio tolerance | 25 % | Allowed drift in the 116 : 118 : 120 ratios, as a share of the expected value. |
   | Minimum intensity | 0 | Peaks below this are ignored. |
   | Also test 117Sn and 119Sn | off | Adds the two smaller isotopes. |
   | Strong / Possible: peaks allowed missing | 0 / 1 | 3 of 3 peaks for Strong, 2 of 3 for Possible. |
   | Hide echoes of Strong hits | on | A real pattern also looks like a weaker one an isotope lower. This hides those duplicates. |
3. **Results.** Strong hits come first, then Possible hits. Strong rows are
   orange and Possible rows are blue, and the Class column says the same in
   words. You can filter by class and by MFP check result. Click a row to
   zoom on that hit. **Download hit table (CSV)** writes the settings, file
   names and app version at the top of the file.
4. **Spectra** (off by default). Turn on **Show spectrum plots** to see the
   full spectrum and a zoom on each hit. You can switch between raw and
   normalised data, overlay several samples, type in axis limits, lock those
   limits for every plot, and export each plot as PNG, SVG or CSV.

The app remembers your settings between sessions in
`~/.tinbarcode/settings.json`. To keep them somewhere else, set the
`TINBARCODE_PREFS` environment variable to another path.

## How matching works

Every peak above the minimum intensity is tried as the 120Sn peak. For each
other isotope, the expected position is that peak's **observed** m/z plus the
isotope's mass offset from the isotope table. A partner peak counts as found
when:

- it lies within the mass tolerance of that expected position, and
- its height relative to the 120Sn peak is within the ratio tolerance.

Because positions are measured from the observed peak, a calibration shift
that moves the whole spectrum does not matter.

The isotope table in the input file is used when there is one. NIST values
fill in any isotope the file does not list. The **Isotope values used** panel
shows both sets side by side.

### MFP check

For each hit, the app looks up the MFP candidates at the 120Sn peak's m/z and
keeps all of them.

| MFP check | Meaning |
|---|---|
| Sn formula | At least one candidate contains tin. |
| Sn formula (Sn2+ only) | Every tin candidate has two or more Sn atoms, so the single-tin ratios may not apply. |
| No Sn formula | Candidates exist, but none contains tin. |
| No MFP candidates | MFP gave nothing for this peak. |

Every result is a putative assignment and needs checking by eye.

## Notes on the example data

- `Rory Trial.xlsx` is ToF-SIMS data (`LMIG_ToF`), not OrbiSIMS. Its spacing
  errors reach about 6 mDa (around 40 ppm), so the default tolerance is
  0.008 Da. The brief's 2 to 3 ppm would miss real tin ions.
- The spectrum is shifted about 0.045 Da low. Sn⁺ appears at m/z 119.857
  instead of 119.902. The MFP search allowed at most 300 ppm, so most real
  tin ions have **no** Sn formula in the MFP output. That is why many Strong
  hits show "No MFP candidates".
- With the defaults, the app finds 9 Strong hits. These include peaks that are likely Sn⁺ (119.857),
  SnOH⁺ (136.851), SnONa⁺ (158.825) and a Sn₂ cluster near 270.66.
- The file's abundances (32.4 / 24.3 / 14.7 %) differ slightly from NIST
  (32.58 / 24.22 / 14.54 %). The file's values are used.

## Code layout

| File | What it does |
|---|---|
| `tinbarcode/loader.py` | Reads workbooks (finding each table by its header text) and plain peak lists |
| `tinbarcode/isotopes.py` | NIST tin isotopes and the expected pattern |
| `tinbarcode/barcode.py` | `Settings` and `find_barcodes`, the matching logic |
| `tinbarcode/plots.py` | Plotly figures and image export |
| `tinbarcode/export.py` | CSV export with a parameter header |
| `tinbarcode/prefs.py` | Remembers settings between sessions |
| `app.py` | The Streamlit interface |
| `tests/` | Run with `python -m pytest` |

To support a new file format, add a reader to `loader.py` that returns a
`Dataset`. To add a new matching rule, add a field to `Settings` and use it
in `find_barcodes`.

## Open questions for the scientist

- Which instrument and mass tolerance should be the default? The example is
  ToF; the brief assumes OrbiSIMS.
- The brief refers to a "section 10" of open questions, which is missing.
- The labelled spectrum image puts Sn⁺ at m/z 119.884, but this file has it
  at 119.857. Are they from the same run?
- One example image (C₂H₃NOSn⁻) is negative polarity, while the data is
  positive.
- What minimum intensity should be used?
- The app does not yet read raw profile spectra or MATLAB output, because no
  example of either was provided. Downsampling for million-point spectra will
  be needed once it does.
