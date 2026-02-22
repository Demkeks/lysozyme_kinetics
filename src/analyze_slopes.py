#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import re
import zipfile
from dataclasses import dataclass
from pathlib import Path
from statistics import mean, median, stdev
from typing import Dict, Iterable, List, Optional, Sequence, Tuple
from xml.sax.saxutils import escape


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


def safe_number(value: Optional[float]) -> object:
    if value is None:
        return ""
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return ""
    return value


def build_measurement_rows(results: Sequence[MeasurementResult]) -> List[List[object]]:
    rows: List[List[object]] = [
        [
            "cohort",
            "condition_id",
            "replicate",
            "source_file",
            "n_points",
            "segment1_end_time_s",
            "segment3_start_time_s",
            "segment1_slope_abs_per_s",
            "segment1_r2",
            "segment3_window_start_time_s",
            "segment3_window_end_time_s",
            "segment3_slope_abs_per_s",
            "segment3_r2",
            "final_slope_abs_per_s",
        ]
    ]

    for r in sorted(results, key=lambda x: (x.cohort, x.condition_id, x.replicate or 0, x.source_file)):
        rows.append(
            [
                r.cohort,
                r.condition_id,
                r.replicate if r.replicate is not None else "",
                r.source_file,
                r.n_points,
                r.segment1_end_time_s,
                r.segment3_start_time_s,
                r.segment1_fit.slope_per_s,
                r.segment1_fit.r2,
                r.segment3_window_start_time_s,
                r.segment3_window_end_time_s,
                r.segment3_fit.slope_per_s,
                r.segment3_fit.r2,
                r.final_slope_per_s,
            ]
        )
    return rows


def build_summary_rows(results: Sequence[MeasurementResult]) -> List[List[object]]:
    rows: List[List[object]] = [
        [
            "cohort",
            "condition_id",
            "n_repeats",
            "segment1_mean_abs_per_s",
            "segment3_mean_abs_per_s",
            "final_mean_abs_per_s",
            "final_std_abs_per_s",
            "final_ci95_low_abs_per_s",
            "final_ci95_high_abs_per_s",
        ]
    ]

    grouped: Dict[Tuple[str, str], List[MeasurementResult]] = {}
    for r in results:
        grouped.setdefault((r.cohort, r.condition_id), []).append(r)

    for (cohort, condition_id), items in sorted(grouped.items()):
        s1 = [x.segment1_fit.slope_per_s for x in items]
        s3 = [x.segment3_fit.slope_per_s for x in items]
        sf = [x.final_slope_per_s for x in items]

        sf_stats = summarize(sf)

        rows.append(
            [
                cohort,
                condition_id,
                len(items),
                safe_number(mean(s1)),
                safe_number(mean(s3)),
                safe_number(sf_stats[0]),
                safe_number(sf_stats[1]),
                safe_number(sf_stats[2]),
                safe_number(sf_stats[3]),
            ]
        )
    return rows


def col_name(index: int) -> str:
    # 1-based column index -> Excel letters.
    name = ""
    n = index
    while n > 0:
        n, rem = divmod(n - 1, 26)
        name = chr(65 + rem) + name
    return name


def cell_ref(row_idx: int, col_idx: int) -> str:
    return f"{col_name(col_idx)}{row_idx}"


def is_number(value: object) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and not (
        isinstance(value, float) and (math.isnan(value) or math.isinf(value))
    )


def build_sheet_xml(rows: Sequence[Sequence[object]], shared_strings: Dict[str, int]) -> str:
    row_xml: List[str] = []
    max_col = max((len(r) for r in rows), default=1)
    for r_idx, row in enumerate(rows, start=1):
        cells: List[str] = []
        for c_idx, value in enumerate(row, start=1):
            ref = cell_ref(r_idx, c_idx)
            if value is None or value == "":
                continue
            if is_number(value):
                cells.append(f'<c r="{ref}"><v>{value}</v></c>')
            else:
                text = str(value)
                sid = shared_strings[text]
                cells.append(f'<c r="{ref}" t="s"><v>{sid}</v></c>')
        row_xml.append(f'<row r="{r_idx}">{"".join(cells)}</row>')

    dimension = f"A1:{col_name(max_col)}{max(1, len(rows))}"
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        f'<dimension ref="{dimension}"/>'
        '<sheetData>'
        f'{"".join(row_xml)}'
        "</sheetData>"
        "</worksheet>"
    )


def build_shared_strings(sheets: Sequence[Sequence[Sequence[object]]]) -> Tuple[Dict[str, int], str]:
    lookup: Dict[str, int] = {}
    strings: List[str] = []

    for rows in sheets:
        for row in rows:
            for value in row:
                if value is None or value == "":
                    continue
                if is_number(value):
                    continue
                text = str(value)
                if text not in lookup:
                    lookup[text] = len(strings)
                    strings.append(text)

    si = "".join(f"<si><t>{escape(s)}</t></si>" for s in strings)
    xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<sst xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        f'count="{len(strings)}" uniqueCount="{len(strings)}">'
        f"{si}</sst>"
    )
    return lookup, xml


def write_xlsx(path: Path, sheets: Sequence[Tuple[str, Sequence[Sequence[object]]]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    sheet_rows = [rows for _, rows in sheets]
    shared_lookup, shared_xml = build_shared_strings(sheet_rows)
    sheet_xml = [build_sheet_xml(rows, shared_lookup) for rows in sheet_rows]

    workbook_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
        "<sheets>"
        + "".join(
            f'<sheet name="{escape(name)}" sheetId="{idx}" r:id="rId{idx}"/>'
            for idx, (name, _) in enumerate(sheets, start=1)
        )
        + "</sheets></workbook>"
    )

    workbook_rels_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        + "".join(
            f'<Relationship Id="rId{idx}" '
            'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" '
            f'Target="worksheets/sheet{idx}.xml"/>'
            for idx in range(1, len(sheets) + 1)
        )
        + f'<Relationship Id="rId{len(sheets) + 1}" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" '
        'Target="styles.xml"/>'
        + f'<Relationship Id="rId{len(sheets) + 2}" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/sharedStrings" '
        'Target="sharedStrings.xml"/>'
        + "</Relationships>"
    )

    root_rels_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" '
        'Target="xl/workbook.xml"/>'
        "</Relationships>"
    )

    styles_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        '<fonts count="1"><font><sz val="11"/><name val="Calibri"/></font></fonts>'
        '<fills count="2"><fill><patternFill patternType="none"/></fill>'
        '<fill><patternFill patternType="gray125"/></fill></fills>'
        '<borders count="1"><border><left/><right/><top/><bottom/><diagonal/></border></borders>'
        '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>'
        '<cellXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/></cellXfs>'
        '<cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>'
        "</styleSheet>"
    )

    content_types_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
        '<Default Extension="xml" ContentType="application/xml"/>'
        '<Override PartName="/xl/workbook.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
        + "".join(
            f'<Override PartName="/xl/worksheets/sheet{idx}.xml" '
            'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
            for idx in range(1, len(sheets) + 1)
        )
        + '<Override PartName="/xl/styles.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>'
        + '<Override PartName="/xl/sharedStrings.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sharedStrings+xml"/>'
        + "</Types>"
    )

    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("[Content_Types].xml", content_types_xml)
        zf.writestr("_rels/.rels", root_rels_xml)
        zf.writestr("xl/workbook.xml", workbook_xml)
        zf.writestr("xl/_rels/workbook.xml.rels", workbook_rels_xml)
        zf.writestr("xl/styles.xml", styles_xml)
        zf.writestr("xl/sharedStrings.xml", shared_xml)
        for idx, xml in enumerate(sheet_xml, start=1):
            zf.writestr(f"xl/worksheets/sheet{idx}.xml", xml)


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

    measurement_rows = build_measurement_rows(results)
    summary_rows = build_summary_rows(results)
    write_xlsx(
        args.output_xlsx,
        sheets=[
            ("per_measurement", measurement_rows),
            ("summary", summary_rows),
        ],
    )

    print(f"Analyzed files: {len(results)}")
    print(f"Output: {args.output_xlsx}")
    if errors:
        print(f"Files with errors: {len(errors)}")
        for err in errors:
            print(f"- {err}")


if __name__ == "__main__":
    main()
