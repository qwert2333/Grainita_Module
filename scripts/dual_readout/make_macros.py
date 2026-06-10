#!/usr/bin/env python3
"""Generate Geant4 macros for the dual-readout study."""

from __future__ import annotations

import argparse

from common import config_dirs, energy_tag, iter_samples, load_config, repo_path, sample_name


def macro_text(config, sample, sample_index: int) -> str:
    detector = config["detector"]
    seed1 = int(config.get("base_seed1", 2077)) + sample_index * 2
    seed2 = int(config.get("base_seed2", 1010)) + sample_index * 2
    return f"""/random/setSeeds {seed1} {seed2}

/run/numberOfThreads {int(config.get("threads", 1))}

/run/setCut {config.get("cut", "0.1 mm")}

/run/verbose {int(config.get("run_verbose", 1))}
/event/verbose {int(config.get("event_verbose", 0))}
/tracking/verbose {int(config.get("tracking_verbose", 0))}

# Detector construction
/detector/ModuleSize {detector.get("module_size", "1 m")}
/detector/ModuleDepth {detector.get("module_depth", "2 m")}
/detector/BoxNum {int(detector.get("box_num", 1))}
/detector/fiberNum {int(detector.get("fiber_num", 24))}
/detector/ZSegNum {int(detector.get("z_seg_num", 1))}
/detector/ResponseX0 {detector.get("response_x0", "0.856 mm")}
/detector/AttLength {detector.get("att_length", "3.4 mm")}
/detector/ResponseSlope {detector.get("response_slope", 0.93)}
/detector/ResponseIntercept {detector.get("response_intercept", 0.206)}
/detector/ReflectCoeff {detector.get("reflect_coeff", 0.9)}
/detector/ApplyLightResponse {str(detector.get("apply_light_response", False)).lower()}

/run/initialize

# Set particle gun
/gps/particle {sample["particle"]}
/gps/energy {sample["energy_gev"]:g} GeV
/gps/direction 0 0 1
/gps/pos/type Plane
/gps/pos/shape Square
/gps/pos/centre {config.get("gps_position", "0.5 0.5 -150 cm")}
/gps/pos/halfx {config.get("gps_halfx", "0. cm")}
/gps/pos/halfy {config.get("gps_halfy", "0. cm")}

# Run event number
/run/beamOn {sample["events"]}
"""


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="Study config, e.g. configs/dual_readout_smoke.yaml")
    args = parser.parse_args()

    config = load_config(args.config)
    dirs = config_dirs(config)

    rows = []
    for index, sample in enumerate(iter_samples(config)):
        name = sample_name(sample["role"], sample["particle"], sample["energy_gev"])
        path = dirs["macros"] / f"{name}.mac"
        path.write_text(macro_text(config, sample, index))
        rows.append((name, path))

    print(f"Generated {len(rows)} macros in {dirs['macros']}")
    for name, path in rows:
        print(f"  {name}: {path.relative_to(repo_path('.'))}")


if __name__ == "__main__":
    main()
