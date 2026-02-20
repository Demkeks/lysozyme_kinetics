#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Sequence, Tuple

from analyze_slopes import (
    DEFAULT_CANDIDATE_STARTS,
    DEFAULT_WINDOW_SIZE,
    detect_segment_boundaries,
    pick_best_window,
)
from parse_kinetics_txt import parse_file


def resolve_txt_files(input_root: Path, selected_files: Sequence[str]) -> List[Path]:
    resolved: List[Path] = []
    for item in selected_files:
        p = Path(item)
        if p.is_absolute() and p.exists():
            resolved.append(p)
            continue
        candidate = input_root / p
        if candidate.exists():
            resolved.append(candidate)
            continue
        raise FileNotFoundError(f"TXT file not found: {item}")
    return resolved


def list_txt_files(input_root: Path) -> List[Path]:
    return sorted(input_root.glob("*.TXT"))


def quantile(values: Sequence[float], q: float) -> float:
    if not values:
        raise ValueError("Cannot compute quantile of empty sequence.")
    if q <= 0:
        return min(values)
    if q >= 1:
        return max(values)
    sorted_vals = sorted(values)
    pos = (len(sorted_vals) - 1) * q
    left = int(pos)
    right = min(left + 1, len(sorted_vals) - 1)
    frac = pos - left
    return sorted_vals[left] * (1.0 - frac) + sorted_vals[right] * frac


def robust_ylim(values: Sequence[float]) -> Tuple[float, float]:
    # Ignore extreme spikes (e.g., Abs=10 during cuvette handling) for plotting limits.
    q_low = quantile(values, 0.02)
    q_high = quantile(values, 0.98)
    span = max(q_high - q_low, 1e-3)
    pad = span * 0.15
    return q_low - pad, q_high + pad


def plot_one_file(
    txt_path: Path,
    output_dir: Path,
    window_size: int,
    candidate_starts: int,
) -> Path:
    import matplotlib.pyplot as plt

    points = parse_file(txt_path)
    if len(points) < window_size + 2:
        raise ValueError(f"Not enough points in {txt_path}")

    times = [t for t, _ in points]
    values = [a for _, a in points]
    y_min, y_max = robust_ylim(values)

    segment1_end_idx, segment3_start_idx = detect_segment_boundaries(points, window_size=window_size)
    seg3_fit, offset = pick_best_window(
        points,
        segment3_start_idx=segment3_start_idx,
        window_size=window_size,
        candidate_starts=candidate_starts,
    )
    best_start_idx = segment3_start_idx + offset
    best_end_idx = best_start_idx + window_size - 1

    fit_x = times[best_start_idx : best_end_idx + 1]
    fit_y = [seg3_fit.intercept + seg3_fit.slope_per_s * x for x in fit_x]

    output_dir.mkdir(parents=True, exist_ok=True)
    stem = txt_path.stem
    combined_path = output_dir / f"{stem}_plots.png"

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 9), sharex=False)

    # Top subplot: overview.
    ax1.plot(times, values, lw=1.4, label="signal")
    ax1.axvline(times[segment1_end_idx], color="#2ca02c", ls="--", lw=1.2, label="segment1 end")
    ax1.axvline(times[segment3_start_idx], color="#ff7f0e", ls="--", lw=1.2, label="segment3 start")
    ax1.plot(fit_x, fit_y, color="#d62728", lw=2.0, label=f"segment3 fit (R²={seg3_fit.r2:.4f})")
    ax1.scatter(fit_x, values[best_start_idx : best_end_idx + 1], color="#d62728", s=18, zorder=3)
    ax1.set_title(f"{stem} - overview")
    ax1.set_xlabel("Time (s)")
    ax1.set_ylabel("Absorbance")
    ax1.set_ylim(y_min, y_max)
    ax1.legend()
    ax1.grid(alpha=0.25)

    # Close-up around segment 3 approximation window.
    left_idx = max(segment3_start_idx - 5, 0)
    right_idx = min(best_end_idx + 25, len(points) - 1)
    x_zoom = times[left_idx : right_idx + 1]
    y_zoom = values[left_idx : right_idx + 1]

    # Bottom subplot: segment 3 close-up.
    ax2.plot(x_zoom, y_zoom, lw=1.6, label="signal (zoom)")
    ax2.plot(fit_x, fit_y, color="#d62728", lw=2.2, label=f"best 10-pt fit (R²={seg3_fit.r2:.4f})")
    ax2.scatter(fit_x, values[best_start_idx : best_end_idx + 1], color="#d62728", s=20, zorder=3)
    ax2.axvline(times[segment3_start_idx], color="#ff7f0e", ls="--", lw=1.2, label="segment3 start")
    ax2.set_title(f"{stem} - segment 3 close-up")
    ax2.set_xlabel("Time (s)")
    ax2.set_ylabel("Absorbance")
    ax2.set_ylim(y_min, y_max)
    ax2.legend()
    ax2.grid(alpha=0.25)

    fig.tight_layout()
    fig.savefig(combined_path, dpi=150)
    plt.close(fig)

    return combined_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot selected TXT files: overview and segment-3 close-up with linear approximation."
    )
    parser.add_argument("--input-root", type=Path, default=Path("datasets"))
    parser.add_argument("--output-dir", type=Path, default=Path("plots"))
    parser.add_argument("--window-size", type=int, default=DEFAULT_WINDOW_SIZE)
    parser.add_argument("--candidate-starts", type=int, default=DEFAULT_CANDIDATE_STARTS)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--txt-files",
        nargs="+",
        help="Selected TXT files (relative to --input-root or absolute paths).",
    )
    group.add_argument(
        "--all-files",
        action="store_true",
        help="Plot all .TXT files found directly inside --input-root.",
    )
    args = parser.parse_args()

    try:
        import matplotlib  # noqa: F401
    except ImportError as exc:
        raise SystemExit(
            "matplotlib is required for plotting. Install it first, e.g. `pip install matplotlib`."
        ) from exc

    if args.all_files:
        files = list_txt_files(args.input_root)
        if not files:
            raise SystemExit(f"No .TXT files found under: {args.input_root}")
    else:
        files = resolve_txt_files(args.input_root, args.txt_files)

    for txt_path in files:
        output_path = plot_one_file(
            txt_path=txt_path,
            output_dir=args.output_dir,
            window_size=args.window_size,
            candidate_starts=args.candidate_starts,
        )
        print(f"{txt_path} ->")
        print(f"  {output_path}")


if __name__ == "__main__":
    main()
