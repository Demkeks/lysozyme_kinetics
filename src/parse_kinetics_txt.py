#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple


SECTION_HEADERS = {
    "Instrument",
    "Instrument Parameters",
    "Kinetic Data",
    "Peak Detection",
    "Peaks",
    "Data Points",
}

NAME_PATTERN = re.compile(
    r"^(?P<cells>\d+)cells_(?P<buffer>\d+)buffer_(?P<imaz>\d+)imaz_(?P<sample>.+)_(?P<replicate>\d+)$"
)


@dataclass
class ParsedFile:
    run: Dict[str, object]
    peaks: List[Dict[str, object]]
    points: List[Dict[str, object]]


def to_float(value: str) -> Optional[float]:
    value = value.strip()
    if value == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def to_int(value: str) -> Optional[int]:
    num = to_float(value)
    if num is None:
        return None
    return int(num)


def parse_filename_fields(file_path: Path) -> Dict[str, object]:
    stem = file_path.stem
    match = NAME_PATTERN.match(stem)
    if not match:
        return {
            "cells": None,
            "buffer": None,
            "imaz": None,
            "sample_type": None,
            "replicate": None,
        }
    parts = match.groupdict()
    return {
        "cells": to_int(parts["cells"]),
        "buffer": to_int(parts["buffer"]),
        "imaz": to_int(parts["imaz"]),
        "sample_type": parts["sample"],
        "replicate": to_int(parts["replicate"]),
    }


def parse_kv_sections(lines: Sequence[str]) -> Dict[str, Dict[str, str]]:
    section = "header"
    sections: Dict[str, Dict[str, str]] = {section: {}}

    for raw in lines:
        line = raw.rstrip("\n")
        stripped = line.strip()
        if not stripped:
            continue
        if stripped in SECTION_HEADERS:
            section = stripped
            sections.setdefault(section, {})
            continue
        if ":" not in line:
            continue
        key, value = line.split(":", 1)
        key = key.strip()
        value = value.strip().lstrip("\t").strip("'")
        if key:
            sections.setdefault(section, {})[key] = value
    return sections


def parse_tabular_blocks(lines: Sequence[str]) -> Tuple[List[Dict[str, object]], List[Dict[str, object]], Dict[str, object]]:
    peaks: List[Dict[str, object]] = []
    points: List[Dict[str, object]] = []
    kinetic: Dict[str, object] = {}

    i = 0
    while i < len(lines):
        token = lines[i].strip()
        if token == "Kinetic Data" and i + 2 < len(lines):
            header = [h.strip() for h in lines[i + 1].split("\t")]
            row = [c.strip() for c in lines[i + 2].split("\t")]
            if len(header) == len(row):
                for h, c in zip(header, row):
                    kinetic[h] = to_float(c) if h not in {"Start (s)", "End (s)"} else to_float(c)
            i += 3
            continue

        if token == "Peaks" and i + 1 < len(lines):
            i += 1
            header = [h.strip() for h in lines[i].split("\t")]
            i += 1
            while i < len(lines):
                stripped = lines[i].strip()
                if not stripped or stripped in SECTION_HEADERS:
                    break
                cols = [c.strip() for c in lines[i].split("\t")]
                if len(cols) != len(header):
                    i += 1
                    continue
                row = {h: c for h, c in zip(header, cols)}
                peaks.append(
                    {
                        "peak_num": to_int(row.get("Peak #", "")),
                        "start_s": to_float(row.get("Start (s)", "")),
                        "apex_s": to_float(row.get("Apex (s)", "")),
                        "end_s": to_float(row.get("End (s)", "")),
                        "height_abs": to_float(row.get("Height (Abs)", "")),
                        "valley_s": to_float(row.get("Valley (s)", "")),
                        "valley_abs": to_float(row.get("Valley (Abs)", "")),
                    }
                )
                i += 1
            continue

        if token == "Data Points" and i + 1 < len(lines):
            i += 2
            while i < len(lines):
                stripped = lines[i].strip()
                if not stripped:
                    break
                cols = [c.strip() for c in lines[i].split("\t")]
                if len(cols) >= 2:
                    t_s = to_float(cols[0])
                    abs_val = to_float(cols[1])
                    if t_s is not None and abs_val is not None:
                        points.append({"time_s": t_s, "absorbance": abs_val})
                i += 1
            continue

        i += 1

    return peaks, points, kinetic


def parse_file(file_path: Path, input_root: Path) -> ParsedFile:
    lines = file_path.read_text(encoding="utf-8", errors="replace").splitlines()
    sections = parse_kv_sections(lines)
    peaks, points, kinetic = parse_tabular_blocks(lines)
    name_fields = parse_filename_fields(file_path)

    rel_path = file_path.relative_to(input_root).as_posix()
    header = sections.get("header", {})
    instr = sections.get("Instrument Parameters", {})

    run = {
        "run_id": rel_path,
        "cohort": file_path.parent.name,
        "source_file": rel_path,
        "sample": header.get("Sample"),
        "file_name_udt": header.get("File name"),
        "run_date_raw": header.get("Run Date"),
        "operator": header.get("Operator"),
        "measurement_type": instr.get("Measurement Type"),
        "data_mode": instr.get("Data Mode"),
        "scan_time_s": to_float((instr.get("Scan Time") or "").replace(" s", "")),
        "wavelength_nm": to_float((instr.get("Wavelength") or "").replace(" nm", "")),
        "sampling_interval_s": to_float((instr.get("Sampling Interval") or "").replace(" s", "")),
        "delay_s": to_float((instr.get("Delay") or "").replace(" s", "")),
        "slope_per_min": kinetic.get("Slope (/min)"),
        "activity": kinetic.get("Activity"),
        "r": kinetic.get("R"),
        "r2": kinetic.get("R2"),
        "n_peaks": len(peaks),
        "n_points": len(points),
        **name_fields,
    }

    peak_rows = []
    for row in peaks:
        peak_rows.append({"run_id": rel_path, **row})

    point_rows = []
    for row in points:
        point_rows.append({"run_id": rel_path, **row})

    return ParsedFile(run=run, peaks=peak_rows, points=point_rows)


def write_csv(path: Path, rows: Sequence[Dict[str, object]], headers: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(headers))
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def parse_all(input_root: Path, output_dir: Path) -> Tuple[int, int, int, int]:
    files = sorted(input_root.glob("*/*.TXT"))
    runs: List[Dict[str, object]] = []
    peaks: List[Dict[str, object]] = []
    points: List[Dict[str, object]] = []

    for file_path in files:
        parsed = parse_file(file_path, input_root=input_root)
        runs.append(parsed.run)
        peaks.extend(parsed.peaks)
        points.extend(parsed.points)

    run_headers = [
        "run_id",
        "cohort",
        "source_file",
        "sample",
        "file_name_udt",
        "run_date_raw",
        "operator",
        "measurement_type",
        "data_mode",
        "scan_time_s",
        "wavelength_nm",
        "sampling_interval_s",
        "delay_s",
        "slope_per_min",
        "activity",
        "r",
        "r2",
        "n_peaks",
        "n_points",
        "cells",
        "buffer",
        "imaz",
        "sample_type",
        "replicate",
    ]
    peak_headers = ["run_id", "peak_num", "start_s", "apex_s", "end_s", "height_abs", "valley_s", "valley_abs"]
    point_headers = ["run_id", "time_s", "absorbance"]

    write_csv(output_dir / "runs.csv", runs, run_headers)
    write_csv(output_dir / "peaks.csv", peaks, peak_headers)
    write_csv(output_dir / "data_points.csv", points, point_headers)

    return len(files), len(runs), len(peaks), len(points)


def main() -> None:
    parser = argparse.ArgumentParser(description="Parse U-5100 kinetics TXT exports into CSV tables.")
    parser.add_argument("--input-root", type=Path, default=Path("datasets"))
    parser.add_argument("--output-dir", type=Path, default=Path("datasets/processed"))
    args = parser.parse_args()

    files, n_runs, n_peaks, n_points = parse_all(args.input_root, args.output_dir)
    print(f"Parsed files: {files}")
    print(f"Runs: {n_runs}")
    print(f"Peaks: {n_peaks}")
    print(f"Data points: {n_points}")
    print(f"Wrote CSVs to: {args.output_dir}")


if __name__ == "__main__":
    main()
