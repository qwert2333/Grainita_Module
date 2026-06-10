#!/usr/bin/env python3
"""Run generated Geant4 macros and collect ROOT files."""

from __future__ import annotations

import argparse
import shutil
import subprocess
from pathlib import Path

from common import config_dirs, iter_samples, load_config, repo_path, root_energy_label, sample_name


def find_output_root(particle: str, energy_gev: float, before: set[Path]) -> Path:
    label = root_energy_label(energy_gev)
    candidates = []
    for path in repo_path(".").glob(f"GrainitaCalo_*_{particle}_{label}.root"):
        if path not in before:
            candidates.append(path)
    if not candidates:
        for path in repo_path(".").glob(f"GrainitaCalo_*_{particle}_{label}.root"):
            candidates.append(path)
    if not candidates:
        raise SystemExit(f"Could not find ROOT output for particle={particle}, energy={label}")
    return max(candidates, key=lambda p: p.stat().st_mtime)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--executable", default="build/main")
    args = parser.parse_args()

    config = load_config(args.config)
    dirs = config_dirs(config)
    executable = repo_path(args.executable)
    if not executable.exists():
        raise SystemExit(f"Missing executable: {executable}")

    for sample in iter_samples(config):
        name = sample_name(sample["role"], sample["particle"], sample["energy_gev"])
        macro = dirs["macros"] / f"{name}.mac"
        if not macro.exists():
            raise SystemExit(f"Missing macro {macro}. Run make_macros.py first.")

        before = set(repo_path(".").glob("GrainitaCalo_*.root"))
        print(f"Running {name}: {macro}")
        subprocess.run([str(executable), str(macro)], cwd=str(repo_path(".")), check=True)
        produced = find_output_root(sample["particle"], sample["energy_gev"], before)
        destination = dirs["root"] / f"{name}.root"
        if destination.exists():
            destination.unlink()
        shutil.move(str(produced), str(destination))
        print(f"  wrote {destination.relative_to(repo_path('.'))}")


if __name__ == "__main__":
    main()
