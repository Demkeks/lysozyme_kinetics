#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import List, Optional, Tuple


def to_float(value: str) -> Optional[float]:
    value = value.strip()
    if not value:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def parse_data_points(lines: List[str]) -> List[Tuple[float, float]]:
    points: List[Tuple[float, float]] = []
    i = 0
    while i < len(lines):
        if lines[i].strip() != "Data Points":
            i += 1
            continue

        # Skip table header line ("s\tAbs")
        i += 2
        while i < len(lines):
            raw = lines[i].strip()
            if not raw:
                break
            cols = [c.strip() for c in lines[i].split("\t")]
            if len(cols) >= 2:
                time_s = to_float(cols[0])
                absorbance = to_float(cols[1])
                if time_s is not None and absorbance is not None:
                    points.append((time_s, absorbance))
            i += 1
        break
    return points


def parse_file(file_path: Path) -> List[Tuple[float, float]]:
    lines = file_path.read_text(encoding="utf-8", errors="replace").splitlines()
    return parse_data_points(lines)


def write_points_csv(path: Path, points: List[Tuple[float, float]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["time_s", "absorbance"])
        writer.writerows(points)


def parse_all(input_root: Path, output_dir: Path) -> Tuple[int, int]:
    files = sorted(input_root.glob("*.TXT"))
    total_points = 0

    for txt_file in files:
        points = parse_file(txt_file)
        total_points += len(points)

        rel = txt_file.relative_to(input_root)
        out_path = output_dir / rel.with_suffix(".csv")
        write_points_csv(out_path, points)

    return len(files), total_points


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Parse U-5100 kinetics TXT exports to per-file Data Points CSV files."
    )
    parser.add_argument("--input-root", type=Path, default=Path("datasets"))
    parser.add_argument("--output-dir", type=Path, default=Path("datasets/processed"))
    args = parser.parse_args()

    file_count, total_points = parse_all(args.input_root, args.output_dir)
    print(f"Parsed files: {file_count}")
    print(f"Data points: {total_points}")
    print(f"Wrote per-file CSVs to: {args.output_dir}")


if __name__ == "__main__":
    main()
