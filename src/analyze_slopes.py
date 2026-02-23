#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import re
from dataclasses import dataclass
from pathlib import Path
from statistics import mean, median, stdev
from typing import Dict, List, Optional, Sequence, Tuple

import pandas as pd


REPLICATE_RE = re.compile(r"^(?P<condition>.+)_(?P<replicate>\d+)$")
DEFAULT_WINDOW_SIZE = 10
DEFAULT_CANDIDATE_STARTS = 2
DEFAULT_SEGMENT1_TRIM_POINTS = 2

# Two-sided t critical values for alpha=0.05 (95% CI), df=1..30.
T_CRIT_95: Dict[int, float] = {
    1: 12.706,
    2: 4.303,
    3: 3.182,
    4: 2.776,
    5: 2.571,
    6: 2.447,
    7: 2.365,
    8: 2.306,
    9: 2.262,
    10: 2.228,
    11: 2.201,
    12: 2.179,
    13: 2.160,
    14: 2.145,
    15: 2.131,
    16: 2.120,
    17: 2.110,
    18: 2.101,
    19: 2.093,
    20: 2.086,
    21: 2.080,
    22: 2.074,
    23: 2.069,
    24: 2.064,
    25: 2.060,
    26: 2.056,
    27: 2.052,
    28: 2.048,
    29: 2.045,
    30: 2.042,
}


@dataclass
class LinearFit:
    slope_per_s: float
    intercept: float
    r: float
    r2: float


@dataclass
class MeasurementResult:
    cohort: str
    condition_id: str
    replicate: Optional[int]
    source_file: str
    n_points: int
    segment1_end_time_s: float
    segment3_start_time_s: float
    segment1_fit: LinearFit
    segment3_fit: LinearFit
    segment3_window_start_time_s: float
    segment3_window_end_time_s: float
    final_slope_per_s: float


def read_points(csv_path: Path) -> List[Tuple[float, float]]:
    points: List[Tuple[float, float]] = []
    with csv_path.open(encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                t = float(row["time_s"])
                a = float(row["absorbance"])
            except (KeyError, ValueError):
                continue
            points.append((t, a))
    return points


def cluster_indices(indices: Sequence[int], max_gap: int = 2) -> List[Tuple[int, int]]:
    if not indices:
        return []
    clusters: List[Tuple[int, int]] = []
    start = indices[0]
    prev = indices[0]
    for idx in indices[1:]:
        if idx - prev <= max_gap:
            prev = idx
            continue
        clusters.append((start, prev))
        start = idx
        prev = idx
    clusters.append((start, prev))
    return clusters


def linear_fit(points: Sequence[Tuple[float, float]]) -> LinearFit:
    if len(points) < 2:
        raise ValueError("Need at least 2 points for linear fit.")

    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    n = len(points)
    x_mean = sum(xs) / n
    y_mean = sum(ys) / n

    sxx = sum((x - x_mean) ** 2 for x in xs)
    sxy = sum((x - x_mean) * (y - y_mean) for x, y in points)
    syy = sum((y - y_mean) ** 2 for y in ys)

    if sxx == 0:
        raise ValueError("Degenerate X values for linear fit.")

    slope = sxy / sxx
    intercept = y_mean - slope * x_mean
    if syy == 0:
        r = 0.0
        r2 = 0.0
    else:
        r = sxy / math.sqrt(sxx * syy)
        r2 = r * r
    return LinearFit(slope_per_s=slope, intercept=intercept, r=r, r2=r2)


def detect_segment_boundaries(points: Sequence[Tuple[float, float]], window_size: int) -> Tuple[int, int]:
    if len(points) < window_size + 2:
        raise ValueError("Not enough points to detect segments.")

    ys = [p[1] for p in points]
    diffs = [abs(ys[i + 1] - ys[i]) for i in range(len(ys) - 1)]
    diff_median = median(diffs) if diffs else 0.0
    jump_threshold = max(0.4, diff_median * 25.0)

    jump_indices = [i for i, d in enumerate(diffs) if d >= jump_threshold]
    clusters = cluster_indices(jump_indices, max_gap=2)

    if len(clusters) < 2:
        top_jump_indices = sorted(range(len(diffs)), key=lambda i: diffs[i], reverse=True)[:6]
        clusters = cluster_indices(sorted(top_jump_indices), max_gap=2)

    if len(clusters) < 2:
        raise ValueError("Could not detect two jump regions.")

    first_cluster = clusters[0]
    second_cluster: Optional[Tuple[int, int]] = None
    for cluster in clusters[1:]:
        enough_gap = cluster[0] - first_cluster[1] >= 4
        enough_points_after = len(points) - (cluster[1] + 1) >= window_size
        if enough_gap and enough_points_after:
            second_cluster = cluster
            break

    if second_cluster is None:
        raise ValueError("Could not detect transition to segment 3.")

    segment1_end_idx = first_cluster[0]
    segment3_start_idx = second_cluster[1] + 1
    if segment1_end_idx < 1:
        raise ValueError("Segment 1 too short for regression.")
    if len(points) - segment3_start_idx < window_size:
        raise ValueError("Segment 3 too short for regression window.")

    return segment1_end_idx, segment3_start_idx


def pick_best_window(
    points: Sequence[Tuple[float, float]],
    segment3_start_idx: int,
    window_size: int,
    candidate_starts: int,
) -> Tuple[LinearFit, int]:
    max_offset = min(candidate_starts - 1, len(points) - segment3_start_idx - window_size)
    if max_offset < 0:
        raise ValueError("Not enough points in segment 3 for requested window.")

    best_fit: Optional[LinearFit] = None
    best_offset: Optional[int] = None

    for offset in range(max_offset + 1):
        start = segment3_start_idx + offset
        window = points[start : start + window_size]
        fit = linear_fit(window)
        if best_fit is None or fit.r2 > best_fit.r2:
            best_fit = fit
            best_offset = offset

    if best_fit is None or best_offset is None:
        raise ValueError("Could not select best window for segment 3.")
    return best_fit, best_offset


def split_condition_and_replicate(stem: str) -> Tuple[str, Optional[int]]:
    match = REPLICATE_RE.match(stem)
    if not match:
        return stem, None
    condition = match.group("condition")
    replicate = int(match.group("replicate"))
    return condition, replicate


def analyze_file(
    csv_path: Path,
    input_root: Path,
    window_size: int,
    candidate_starts: int,
    segment1_trim_points: int,
) -> MeasurementResult:
    points = read_points(csv_path)
    if len(points) < window_size + 2:
        raise ValueError(f"Not enough points in {csv_path}")

    segment1_end_idx, segment3_start_idx = detect_segment_boundaries(points, window_size=window_size)
    segment1_fit_end_idx = max(1, segment1_end_idx - segment1_trim_points)

    segment1_points = points[: segment1_fit_end_idx + 1]
    segment1_fit = linear_fit(segment1_points)

    segment3_fit, offset = pick_best_window(
        points,
        segment3_start_idx=segment3_start_idx,
        window_size=window_size,
        candidate_starts=candidate_starts,
    )
    best_start_idx = segment3_start_idx + offset
    best_end_idx = best_start_idx + window_size - 1

    final_slope = abs(segment3_fit.slope_per_s) - abs(segment1_fit.slope_per_s)

    rel = csv_path.relative_to(input_root)
    cohort = rel.parts[0] if len(rel.parts) > 1 else "unknown"
    condition_id, replicate = split_condition_and_replicate(csv_path.stem)

    return MeasurementResult(
        cohort=cohort,
        condition_id=condition_id,
        replicate=replicate,
        source_file=rel.as_posix(),
        n_points=len(points),
        segment1_end_time_s=points[segment1_fit_end_idx][0],
        segment3_start_time_s=points[segment3_start_idx][0],
        segment1_fit=segment1_fit,
        segment3_fit=segment3_fit,
        segment3_window_start_time_s=points[best_start_idx][0],
        segment3_window_end_time_s=points[best_end_idx][0],
        final_slope_per_s=final_slope,
    )


def t_critical_95(df: int) -> float:
    if df <= 0:
        return float("nan")
    if df in T_CRIT_95:
        return T_CRIT_95[df]
    return 1.96


def summarize(values: Sequence[float]) -> Tuple[float, Optional[float], Optional[float], Optional[float]]:
    m = mean(values)
    if len(values) < 2:
        return m, None, None, None
    s = stdev(values)
    half_width = t_critical_95(len(values) - 1) * s / math.sqrt(len(values))
    return m, s, m - half_width, m + half_width


def build_measurement_df(results: Sequence[MeasurementResult]) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    for r in sorted(results, key=lambda x: (x.cohort, x.condition_id, x.replicate or 0, x.source_file)):
        rows.append(
            {
                "cohort": r.cohort,
                "condition_id": r.condition_id,
                "replicate": r.replicate,
                "source_file": r.source_file,
                "n_points": r.n_points,
                "segment1_end_time_s": r.segment1_end_time_s,
                "segment3_start_time_s": r.segment3_start_time_s,
                "segment1_slope_abs_per_s": r.segment1_fit.slope_per_s,
                "segment1_r2": r.segment1_fit.r2,
                "segment3_window_start_time_s": r.segment3_window_start_time_s,
                "segment3_window_end_time_s": r.segment3_window_end_time_s,
                "segment3_slope_abs_per_s": r.segment3_fit.slope_per_s,
                "segment3_r2": r.segment3_fit.r2,
                "final_slope_abs_per_s": r.final_slope_per_s,
            }
        )
    return pd.DataFrame(rows)


def build_summary_df(results: Sequence[MeasurementResult]) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    grouped: Dict[Tuple[str, str], List[MeasurementResult]] = {}
    for r in results:
        grouped.setdefault((r.cohort, r.condition_id), []).append(r)

    for (cohort, condition_id), items in sorted(grouped.items()):
        s1 = [x.segment1_fit.slope_per_s for x in items]
        s3 = [x.segment3_fit.slope_per_s for x in items]
        sf = [x.final_slope_per_s for x in items]
        sf_stats = summarize(sf)

        rows.append(
            {
                "cohort": cohort,
                "condition_id": condition_id,
                "n_repeats": len(items),
                "segment1_mean_abs_per_s": mean(s1),
                "segment3_mean_abs_per_s": mean(s3),
                "final_mean_abs_per_s": sf_stats[0],
                "final_std_abs_per_s": sf_stats[1],
                "final_ci95_low_abs_per_s": sf_stats[2],
                "final_ci95_high_abs_per_s": sf_stats[3],
            }
        )
    return pd.DataFrame(rows)


def find_input_files(input_root: Path) -> List[Path]:
    return sorted(input_root.glob("*.csv"))


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compute segment slopes and summary statistics from per-file Data Points CSVs."
    )
    parser.add_argument("--input-root", type=Path, default=Path("datasets/processed"))
    parser.add_argument("--output-xlsx", type=Path, default=Path("datasets/processed/slopes_analysis.xlsx"))
    parser.add_argument("--window-size", type=int, default=DEFAULT_WINDOW_SIZE)
    parser.add_argument(
        "--segment1-trim-points",
        type=int,
        default=DEFAULT_SEGMENT1_TRIM_POINTS,
        help="Drop last N points of segment 1 (before the first jump) before linear regression.",
    )
    parser.add_argument(
        "--candidate-starts",
        type=int,
        default=DEFAULT_CANDIDATE_STARTS,
        help="Try window starts within first N points of segment 3 and pick the fit with max R^2.",
    )
    args = parser.parse_args()
    if args.segment1_trim_points < 0:
        raise SystemExit("--segment1-trim-points must be >= 0")

    input_files = find_input_files(args.input_root)
    if not input_files:
        raise SystemExit(f"No input CSV files found under: {args.input_root}")

    results: List[MeasurementResult] = []
    errors: List[str] = []
    for csv_path in input_files:
        try:
            result = analyze_file(
                csv_path,
                input_root=args.input_root,
                window_size=args.window_size,
                candidate_starts=args.candidate_starts,
                segment1_trim_points=args.segment1_trim_points,
            )
            results.append(result)
        except Exception as exc:
            errors.append(f"{csv_path.as_posix()}: {exc}")

    if not results:
        raise SystemExit("No files were successfully analyzed.")

    measurement_df = build_measurement_df(results)
    summary_df = build_summary_df(results)
    args.output_xlsx.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(args.output_xlsx, engine="openpyxl") as writer:
        measurement_df.to_excel(writer, sheet_name="per_measurement", index=False)
        summary_df.to_excel(writer, sheet_name="summary", index=False)

    print(f"Analyzed files: {len(results)}")
    print(f"Output: {args.output_xlsx}")
    if errors:
        print(f"Files with errors: {len(errors)}")
        for err in errors:
            print(f"- {err}")


if __name__ == "__main__":
    main()
