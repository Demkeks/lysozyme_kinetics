# Lysozyme Kinetics

Repository for lysozyme kinetics measurements exported from a Hitachi U-5100 spectrophotometer.

## Structure

- `datasets/arbonal/` - measurements from arbonal set
- `datasets/humat/` - measurements from humat set
- `src/` - parsing scripts

## Parse raw TXT files

Generate normalized CSV tables from all `.TXT` files:

```bash
python3 src/parse_kinetics_txt.py --input-root datasets --output-dir datasets/processed
```

This command writes:

- `datasets/processed/runs.csv` - run-level metadata and kinetic summary
- `datasets/processed/peaks.csv` - peak table
- `datasets/processed/data_points.csv` - time series (`time_s`, `absorbance`)

## Notes

- Raw files are kept unchanged.
- Output folder is ignored by git (`datasets/processed/`).
