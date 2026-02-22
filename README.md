# Lysozyme Kinetics

Repository for lysozyme kinetics measurements exported from a Hitachi U-5100 spectrophotometer.

## Dependencies

This project is built with Python 3.12.3.
To run plot_selected_txt.py, create a virtual environment and install the dependencies (including matplotlib 3.10.8):

```bash
pip install -r requirements.txt
```

## Structure

- `data_samples/input/` - example of raw .TXT files
- `data_samples/processed/` - results of parse_kinetics.py and analyze_slopes.py
- `data_samples/plots/` - resulting plots, result of plot_selected_txt.py

## Parse raw TXT files

Generate per-file Data Points CSV from all `.TXT` files:

```bash
python3 src/parse_kinetics_txt.py --input-root data_samples/input --output-dir data_samples/processed
```

This command writes one CSV per input file with mirrored folder structure, e.g.:

- `datasets/processed/arbonal/10cells_970buffer_20imaz_lsz_1.csv`
- `datasets/processed/humat/10cells_970buffer_20imaz_hum_lsz_1.csv`

Each output CSV contains only:

- `time_s`
- `absorbance`

## Analyze slopes and export XLSX

Compute:

- slope of segment 1 (baseline)
- slope of segment 3 using the best 10-point linear window (max `R^2`, window start is searched within the first 2 points of segment 3)
- final slope: `abs(segment3_slope) - abs(segment1_slope)` (positive rate)
- all slopes are reported in `Abs/s`
- summary stats across replicates (`_1`, `_2`, `_3`):
  - `segment1` and `segment3`: mean only
  - `final`: mean, standard deviation, and 95% confidence interval

Command:

```bash
python3 src/analyze_slopes.py --input-root data_samples/processed --output-xlsx data_samples/processed/final.xlsx
```

Optional robustness flag for baseline fit:

```bash
--segment1-trim-points 2
```

This excludes the last N points of segment 1 before fitting its line.

Workbook sheets:

- `per_measurement` - one row per measurement file
- `summary` - aggregated statistics per condition

## Plot selected TXT files

Build plots for selected `.TXT` files:

- figure 1: full signal overview
- figure 2: close-up around segment 3 with best 10-point linear fit
- both figures are saved in a single PNG (two subplots)

Command example:

```bash
python3 src/plot_selected_txt.py \
  --input-root data_samples/input \
  --output-dir data_samples/plots \
  --txt-files 10cells_970buffer_20imaz_lsz_1.TXT 20cells_960buffer_20imaz_lsz_1.TXT 30cells_950buffer_20imaz_lsz_3.TXT
```

All files from a directory:

```bash
python3 src/plot_selected_txt.py \
  --input-root data_samples/input \
  --output-dir data_samples/plots \
  --all-files
```

Output per selected file:

- `<name>_plots.png`

## Notes

- Raw files are kept unchanged.
