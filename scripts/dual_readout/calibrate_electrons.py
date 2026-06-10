#!/usr/bin/env python3
"""Calibrate scintillation and Cherenkov channels with electron samples."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

from common import (
    find_root_files,
    finite,
    linear_slope_through_origin,
    load_events,
    mean,
    read_csv,
    require_columns,
    rms,
    stderr,
    write_csv,
    write_json,
    repo_path,
    ensure_dir,
)


def summarize_file(path: Path) -> dict:
    events = load_events(path)
    require_columns(events, ["MCtruth_energy", "Nph_Scintillation", "Nph_Cherenkov", "truthFemCrystal"], path)
    valid = [
        row for row in events
        if float(row["MCtruth_energy"]) > 0
        and float(row["Nph_Scintillation"]) > 0
        and float(row["Nph_Cherenkov"]) > 0
    ]
    if not valid:
        raise SystemExit(f"No valid electron calibration events in {path}")

    energies = [float(row["MCtruth_energy"]) for row in valid]
    s_values = [float(row["Nph_Scintillation"]) for row in valid]
    c_values = [float(row["Nph_Cherenkov"]) for row in valid]
    fem_values = [float(row["truthFemCrystal"]) for row in valid if float(row["truthFemCrystal"]) >= 0]
    energy = mean(energies)
    mean_s = mean(s_values)
    mean_c = mean(c_values)
    return {
        "file": str(path.relative_to(repo_path("."))),
        "energy_MeV": energy,
        "energy_GeV": energy / 1000.0,
        "entries": len(events),
        "valid_entries": len(valid),
        "mean_S_raw": mean_s,
        "rms_S_raw": rms(s_values),
        "stderr_S_raw": stderr(s_values),
        "mean_C_raw": mean_c,
        "rms_C_raw": rms(c_values),
        "stderr_C_raw": stderr(c_values),
        "mean_truthFemCrystal": mean(fem_values) if fem_values else math.nan,
        "kS_point": energy / mean_s,
        "kC_point": energy / mean_c,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Directory containing electron ROOT files")
    parser.add_argument("--out", required=True, help="Output table directory")
    args = parser.parse_args()

    input_dir = repo_path(args.input)
    out_dir = ensure_dir(args.out)
    files = find_root_files(input_dir, role="electron")
    if not files:
        files = find_root_files(input_dir, particle="e-")
    if not files:
        raise SystemExit(f"No electron ROOT files found in {input_dir}")

    rows = [summarize_file(path) for path in files]
    rows.sort(key=lambda row: row["energy_MeV"])
    energies = [float(row["energy_MeV"]) for row in rows]
    mean_s = [float(row["mean_S_raw"]) for row in rows]
    mean_c = [float(row["mean_C_raw"]) for row in rows]

    slope_s = linear_slope_through_origin(energies, mean_s)
    slope_c = linear_slope_through_origin(energies, mean_c)
    a_s = 1.0 / slope_s if finite(slope_s) and slope_s > 0 else math.nan
    a_c = 1.0 / slope_c if finite(slope_c) and slope_c > 0 else math.nan

    for row in rows:
        row["S_cal_over_E_mean"] = a_s * float(row["mean_S_raw"]) / float(row["energy_MeV"])
        row["C_cal_over_E_mean"] = a_c * float(row["mean_C_raw"]) / float(row["energy_MeV"])

    calibration = {
        "signal_definitions": {
            "S_raw": "Nph_Scintillation",
            "C_raw": "Nph_Cherenkov",
        },
        "calibration": {
            "aS_MeV_per_count": a_s,
            "aC_MeV_per_count": a_c,
            "slope_S_counts_per_MeV": slope_s,
            "slope_C_counts_per_MeV": slope_c,
        },
        "points": rows,
    }

    write_csv(out_dir / "electron_calibration.csv", rows)
    write_json(out_dir / "electron_calibration.json", calibration)
    print(f"Wrote {out_dir / 'electron_calibration.csv'}")
    print(f"Wrote {out_dir / 'electron_calibration.json'}")
    print(f"aS = {a_s:.8g} MeV/count, aC = {a_c:.8g} MeV/count")


if __name__ == "__main__":
    main()
