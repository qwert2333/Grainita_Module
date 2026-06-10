#!/usr/bin/env python3
"""Create standard plots for the dual-readout study."""

from __future__ import annotations

import argparse
import os
import tempfile

cache_dir = tempfile.mkdtemp(prefix="dual-readout-plot-cache-")
if not os.environ.get("MPLCONFIGDIR"):
    os.environ["MPLCONFIGDIR"] = cache_dir
if not os.environ.get("XDG_CACHE_HOME"):
    os.environ["XDG_CACHE_HOME"] = cache_dir

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from common import ensure_dir, read_csv, read_json, repo_path, rms


def savefig(path):
    plt.tight_layout()
    plt.savefig(path, dpi=160)
    plt.close()
    print(f"wrote {path}")


def plot_electron_linearity(tables, out):
    rows = read_csv(tables / "electron_calibration.csv")
    energies = [float(row["energy_MeV"]) / 1000.0 for row in rows]
    s = [float(row["mean_S_raw"]) for row in rows]
    c = [float(row["mean_C_raw"]) for row in rows]
    plt.figure()
    plt.plot(energies, s, "o-", label="Scintillation")
    plt.plot(energies, c, "s-", label="Cherenkov")
    plt.xlabel("Beam energy [GeV]")
    plt.ylabel("Mean raw signal [counts]")
    plt.legend()
    plt.grid(True, alpha=0.3)
    savefig(out / "electron_raw_linearity.png")

    plt.figure()
    plt.plot(energies, [float(row["S_cal_over_E_mean"]) for row in rows], "o-", label="S calibrated")
    plt.plot(energies, [float(row["C_cal_over_E_mean"]) for row in rows], "s-", label="C calibrated")
    plt.axhline(1.0, color="black", linewidth=1)
    plt.xlabel("Beam energy [GeV]")
    plt.ylabel("Mean calibrated response / E")
    plt.legend()
    plt.grid(True, alpha=0.3)
    savefig(out / "electron_calibrated_closure.png")


def plot_response_vs_fem(tables, out):
    rows = read_csv(tables / "event_reco.csv")
    response = read_csv(tables / "hadron_response.csv")
    fit = {(row["subset"], row["channel"]): row for row in response}
    fem = [float(row["truthFemCrystal"]) for row in rows]
    s = [float(row["S_cal_over_E"]) for row in rows]
    c = [float(row["C_cal_over_E"]) for row in rows]
    xline = [0.0, 1.0]

    plt.figure()
    plt.scatter(fem, s, s=12, alpha=0.65, label="S events")
    plt.scatter(fem, c, s=12, alpha=0.65, label="C events")
    if ("all", "S") in fit:
        h = float(fit[("all", "S")]["h_over_e"])
        plt.plot(xline, [x + h * (1 - x) for x in xline], label=f"S fit h/e={h:.3g}")
    if ("all", "C") in fit:
        h = float(fit[("all", "C")]["h_over_e"])
        plt.plot(xline, [x + h * (1 - x) for x in xline], label=f"C fit h/e={h:.3g}")
    plt.xlabel("truth f_em in crystal")
    plt.ylabel("Calibrated response / E")
    plt.legend()
    plt.grid(True, alpha=0.3)
    savefig(out / "hadron_response_vs_truth_fem.png")


def plot_h_over_e_vs_energy(tables, out):
    path = tables / "h_over_e_vs_energy.csv"
    if not path.exists():
        print(f"skip {path}: run fit_h_over_e_by_energy.py first")
        return
    rows = read_csv(path)
    energies = [float(row["energy_GeV"]) for row in rows]
    s = [float(row["h_over_e_S"]) for row in rows]
    s_err = [float(row["h_over_e_S_err"]) for row in rows]
    c = [float(row["h_over_e_C"]) for row in rows]
    c_err = [float(row["h_over_e_C_err"]) for row in rows]

    plt.figure(figsize=(7.2, 5.2))
    plt.errorbar(energies, s, yerr=s_err, fmt="o-", capsize=4, label="(h/e)_S")
    plt.errorbar(energies, c, yerr=c_err, fmt="s-", capsize=4, label="(h/e)_C")
    plt.xlabel("Beam energy [GeV]")
    plt.ylabel("h/e")
    plt.xticks(energies)
    plt.ylim(0.0, 1.05)
    plt.legend()
    plt.grid(True, alpha=0.3)
    savefig(out / "h_over_e_vs_energy.png")


def plot_energy_distributions(tables, out):
    rows = read_csv(tables / "event_reco.csv")
    values = {
        "S": [float(row["S_cal_over_E"]) for row in rows],
        "C": [float(row["C_cal_over_E"]) for row in rows],
        "DR": [float(row["E_DR_over_E"]) for row in rows],
    }
    plt.figure()
    for label, data in values.items():
        data = [x for x in data if x == x]
        plt.hist(data, bins=30, histtype="step", linewidth=1.8, label=label)
    plt.xlabel("Reconstructed energy / beam energy")
    plt.ylabel("Events")
    plt.legend()
    plt.grid(True, alpha=0.3)
    savefig(out / "energy_response_distributions.png")


def plot_fem_reco(tables, out):
    rows = read_csv(tables / "event_reco.csv")
    truth = [float(row["truthFemCrystal"]) for row in rows if float(row["fem_DR"]) == float(row["fem_DR"])]
    reco = [float(row["fem_DR"]) for row in rows if float(row["fem_DR"]) == float(row["fem_DR"])]
    plt.figure()
    plt.scatter(truth, reco, s=12, alpha=0.7)
    plt.plot([0, 1], [0, 1], color="black", linewidth=1)
    plt.xlabel("truth f_em in crystal")
    plt.ylabel("dual-readout reconstructed f_em")
    plt.grid(True, alpha=0.3)
    savefig(out / "reco_fem_vs_truth_fem.png")


def plot_fem_distributions_by_sample(tables, out):
    rows = read_csv(tables / "event_reco.csv")
    energies = sorted({float(row["energy_GeV"]) for row in rows})
    bins = [(-0.5 + i * 0.025) for i in range(81)]

    for energy in energies:
        sample = [
            row for row in rows
            if float(row["energy_GeV"]) == energy
            and float(row["fem_DR"]) == float(row["fem_DR"])
        ]
        truth = [float(row["truthFemCrystal"]) for row in sample]
        reco = [float(row["fem_DR"]) for row in sample]

        plt.figure(figsize=(7.2, 5.2))
        plt.hist(
            truth,
            bins=bins,
            histtype="step",
            linewidth=2.0,
            label=f"truth f_em (RMS={rms(truth):.3f})",
        )
        plt.hist(
            reco,
            bins=bins,
            histtype="step",
            linewidth=2.0,
            label=f"reco f_em (RMS={rms(reco):.3f})",
        )
        plt.axvline(0.0, color="black", linewidth=0.8, alpha=0.5)
        plt.axvline(1.0, color="black", linewidth=0.8, alpha=0.5)
        plt.xlabel("f_em")
        plt.ylabel("Events")
        plt.title(f"pi- {energy:g} GeV, truth f_em > 0 (N={len(sample)})")
        plt.xlim(-0.5, 1.5)
        plt.legend()
        plt.grid(True, alpha=0.3)
        savefig(out / f"fem_truth_reco_pim_{energy:g}GeV.png")


def plot_leakage_residual(tables, out):
    rows = read_csv(tables / "event_reco.csv")
    leakage = [float(row["leakageFraction"]) for row in rows]
    residual = [float(row["E_DR_over_E"]) - 1.0 for row in rows]
    plt.figure()
    plt.scatter(leakage, residual, s=12, alpha=0.7)
    plt.axhline(0.0, color="black", linewidth=1)
    plt.xlabel("Leakage energy / beam energy")
    plt.ylabel("E_DR / E - 1")
    plt.grid(True, alpha=0.3)
    savefig(out / "leakage_fraction_vs_dr_residual.png")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tables", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    tables = repo_path(args.tables)
    out = ensure_dir(args.out)
    plot_electron_linearity(tables, out)
    plot_response_vs_fem(tables, out)
    plot_h_over_e_vs_energy(tables, out)
    plot_energy_distributions(tables, out)
    plot_fem_reco(tables, out)
    plot_fem_distributions_by_sample(tables, out)
    plot_leakage_residual(tables, out)


if __name__ == "__main__":
    main()
