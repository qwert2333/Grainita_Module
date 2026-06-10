#!/usr/bin/env python3
"""Analyze hadron samples with electron calibration constants."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

from common import (
    find_root_files,
    finite,
    fit_response_vs_fem,
    load_events,
    mean,
    read_json,
    require_columns,
    rms,
    write_csv,
    write_json,
    repo_path,
    ensure_dir,
)


def reco_events(path: Path, a_s: float, a_c: float) -> tuple[list[dict], dict]:
    events = load_events(path)
    require_columns(
        events,
        ["MCtruth_energy", "Nph_Scintillation", "Nph_Cherenkov", "truthFemCrystal", "leakageEnergy"],
        path,
    )
    rows = []
    rejected = {
        "truthFem_lt_0": 0,
        "truthFem_eq_0": 0,
        "S_raw_le_0": 0,
        "C_raw_le_0": 0,
        "energy_le_0": 0,
    }
    for row in events:
        energy = float(row["MCtruth_energy"])
        fem = float(row["truthFemCrystal"])
        s_raw = float(row["Nph_Scintillation"])
        c_raw = float(row["Nph_Cherenkov"])
        if energy <= 0:
            rejected["energy_le_0"] += 1
            continue
        if fem < 0:
            rejected["truthFem_lt_0"] += 1
            continue
        if fem == 0:
            rejected["truthFem_eq_0"] += 1
            continue
        if s_raw <= 0:
            rejected["S_raw_le_0"] += 1
            continue
        if c_raw <= 0:
            rejected["C_raw_le_0"] += 1
            continue

        s_cal = a_s * s_raw
        c_cal = a_c * c_raw
        leakage = float(row["leakageEnergy"])
        rows.append({
            "source_file": str(path.relative_to(repo_path("."))),
            "eventID": int(row["eventID"]),
            "particle": row["particle"],
            "energy_MeV": energy,
            "energy_GeV": energy / 1000.0,
            "S_raw": s_raw,
            "C_raw": c_raw,
            "S_cal": s_cal,
            "C_cal": c_cal,
            "S_cal_over_E": s_cal / energy,
            "C_cal_over_E": c_cal / energy,
            "truthFemCrystal": fem,
            "leakageEnergy": leakage,
            "leakageFraction": leakage / energy,
            "EdepCrystal": float(row["EdepCrystal"]),
            "truthEdepCrystalTotal": float(row["truthEdepCrystalTotal"]),
        })
    return rows, rejected


def summarize_subset(rows: list[dict], label: str, h_s: float, h_c: float) -> dict:
    if not rows:
        return {"subset": label, "n_events": 0}
    chi = (1.0 - h_s) / (1.0 - h_c) if abs(1.0 - h_c) > 1e-12 else math.nan
    for row in rows:
        energy = float(row["energy_MeV"])
        s_cal = float(row["S_cal"])
        c_cal = float(row["C_cal"])
        if finite(chi) and abs(1.0 - chi) > 1e-12:
            e_dr = (s_cal - chi * c_cal) / (1.0 - chi)
        else:
            e_dr = math.nan
        row["chi"] = chi
        row["E_DR"] = e_dr
        row["E_DR_over_E"] = e_dr / energy if finite(e_dr) and energy > 0 else math.nan
        if finite(e_dr) and abs(1.0 - h_s) > 1e-12 and abs(e_dr) > 1e-12:
            row["fem_DR"] = (s_cal / e_dr - h_s) / (1.0 - h_s)
        else:
            row["fem_DR"] = math.nan

    e_s = [float(row["S_cal_over_E"]) for row in rows]
    e_c = [float(row["C_cal_over_E"]) for row in rows]
    e_dr = [float(row["E_DR_over_E"]) for row in rows if finite(float(row["E_DR_over_E"]))]
    fem_delta = [
        float(row["fem_DR"]) - float(row["truthFemCrystal"])
        for row in rows
        if finite(float(row["fem_DR"]))
    ]
    return {
        "subset": label,
        "n_events": len(rows),
        "h_over_e_S": h_s,
        "h_over_e_C": h_c,
        "chi": chi,
        "mean_S_over_E": mean(e_s),
        "rms_S_over_E": rms(e_s),
        "mean_C_over_E": mean(e_c),
        "rms_C_over_E": rms(e_c),
        "mean_DR_over_E": mean(e_dr),
        "rms_DR_over_E": rms(e_dr),
        "resolution_DR": rms(e_dr) / mean(e_dr) if e_dr and mean(e_dr) != 0 else math.nan,
        "mean_fem_DR_minus_truth": mean(fem_delta) if fem_delta else math.nan,
        "rms_fem_DR_minus_truth": rms(fem_delta) if fem_delta else math.nan,
    }


def fit_and_summarize(rows: list[dict], label: str) -> tuple[dict, list[dict]]:
    fem = [float(row["truthFemCrystal"]) for row in rows]
    s_response = [float(row["S_cal_over_E"]) for row in rows]
    c_response = [float(row["C_cal_over_E"]) for row in rows]
    fit_s = fit_response_vs_fem(fem, s_response)
    fit_c = fit_response_vs_fem(fem, c_response)
    summary = summarize_subset(rows, label, fit_s["h_over_e"], fit_c["h_over_e"])
    summary.update({
        "h_over_e_S_err": fit_s["h_over_e_err"],
        "h_over_e_C_err": fit_c["h_over_e_err"],
        "fit_rms_S": fit_s["fit_rms"],
        "fit_rms_C": fit_c["fit_rms"],
    })
    response_rows = [
        {"subset": label, "channel": "S", **fit_s},
        {"subset": label, "channel": "C", **fit_c},
    ]
    return summary, response_rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Directory containing hadron ROOT files")
    parser.add_argument("--calib", required=True, help="electron_calibration.json")
    parser.add_argument("--out", required=True, help="Output table directory")
    parser.add_argument("--contained-leakage-max", type=float, default=0.05)
    args = parser.parse_args()

    calibration = read_json(args.calib)["calibration"]
    a_s = float(calibration["aS_MeV_per_count"])
    a_c = float(calibration["aC_MeV_per_count"])
    if not finite(a_s) or not finite(a_c):
        raise SystemExit("Calibration constants are not finite")

    input_dir = repo_path(args.input)
    out_dir = ensure_dir(args.out)
    files = [path for path in find_root_files(input_dir) if "electron_" not in path.name]
    if not files:
        raise SystemExit(f"No hadron ROOT files found in {input_dir}")

    all_rows = []
    rejected_by_file = {}
    for path in files:
        rows, rejected = reco_events(path, a_s, a_c)
        all_rows.extend(rows)
        rejected_by_file[str(path.relative_to(repo_path(".")))] = rejected
    if not all_rows:
        raise SystemExit("No valid hadron events after quality cuts")

    contained = [row for row in all_rows if float(row["leakageFraction"]) < args.contained_leakage_max]
    summary_all, response_all = fit_and_summarize(all_rows, "all")
    if contained:
        summary_contained, response_contained = fit_and_summarize(contained, "contained")
    else:
        summary_contained, response_contained = ({"subset": "contained", "n_events": 0}, [])

    response_rows = response_all + response_contained
    summary = {
        "calibration_file": str(repo_path(args.calib).relative_to(repo_path("."))),
        "contained_leakage_max": args.contained_leakage_max,
        "rejected_by_file": rejected_by_file,
        "subsets": [summary_all, summary_contained],
    }

    write_csv(out_dir / "hadron_response.csv", response_rows)
    write_csv(out_dir / "event_reco.csv", all_rows)
    write_json(out_dir / "dual_readout_summary.json", summary)
    print(f"Wrote {out_dir / 'hadron_response.csv'}")
    print(f"Wrote {out_dir / 'event_reco.csv'}")
    print(f"Wrote {out_dir / 'dual_readout_summary.json'}")
    print(f"all-events h/e S={summary_all['h_over_e_S']:.5g}, C={summary_all['h_over_e_C']:.5g}")
    if contained:
        print(f"contained h/e S={summary_contained['h_over_e_S']:.5g}, C={summary_contained['h_over_e_C']:.5g}")


if __name__ == "__main__":
    main()
