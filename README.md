# Granita Module simulation model.

## Build

```bash
conda run -n physics cmake --build build
```

The hand-edited macro `run/run.mac` is useful for quick detector and output
checks. The reusable dual-readout study below uses generated macros instead,
so production samples do not depend on the current manual state of `run/run.mac`.

## Dual Readout Study

The analysis workflow is driven by Python scripts in `scripts/dual_readout/`.
Run them through the `physics` conda environment:

```bash
conda run -n physics python scripts/dual_readout/make_macros.py --config configs/dual_readout_smoke.yaml
conda run -n physics python scripts/dual_readout/run_samples.py --config configs/dual_readout_smoke.yaml
conda run -n physics python scripts/dual_readout/calibrate_electrons.py --input outputs/dual_readout/root --out outputs/dual_readout/tables
conda run -n physics python scripts/dual_readout/analyze_hadrons.py --input outputs/dual_readout/root --calib outputs/dual_readout/tables/electron_calibration.json --out outputs/dual_readout/tables
conda run -n physics python scripts/dual_readout/fit_h_over_e_by_energy.py --input outputs/dual_readout/root --calib outputs/dual_readout/tables/electron_calibration.json --out outputs/dual_readout/tables
conda run -n physics python scripts/dual_readout/plot_dual_readout.py --tables outputs/dual_readout/tables --out outputs/dual_readout/plots
```

For full production, use the same commands with
`configs/dual_readout_full.yaml`. The default full configuration produces
`e-` and `pi-` samples at `1, 2, 5, 10, 20 GeV`; electrons are used for channel
calibration, and pions are used for hadron response and `h/e` extraction.

The config files are `.yaml` files written with JSON syntax. JSON is valid
YAML, and this avoids adding a PyYAML dependency to the current environment.

### Outputs

Generated files are written under `outputs/dual_readout/`:

- `macros/`: generated Geant4 macros, one per particle and energy point.
- `root/`: collected ROOT outputs from `build/main`.
- `tables/electron_calibration.csv`: per-electron-energy raw signal means,
  RMS values, and point calibration factors.
- `tables/electron_calibration.json`: global calibration constants
  `aS_MeV_per_count` and `aC_MeV_per_count`.
- `tables/hadron_response.csv`: fitted `(h/e)_S` and `(h/e)_C` for all events
  and leakage-contained events.
- `tables/h_over_e_vs_energy.csv`: independent `(h/e)_S` and `(h/e)_C` fits at
  every pion beam energy, including `truthFemCrystal == 0` events.
- `tables/event_reco.csv`: event-level calibrated `S`, `C`, dual-readout
  energy `E_DR`, reconstructed `f_em_DR`, truth `f_em`, and leakage fraction.
- `tables/dual_readout_summary.json`: analysis summary, selections, rejected
  event counts, response means, RMS values, and dual-readout resolution.
- `plots/electron_raw_linearity.png`: electron raw scintillation and
  Cherenkov signal linearity.
- `plots/electron_calibrated_closure.png`: electron calibrated response
  closure after applying global calibration.
- `plots/hadron_response_vs_truth_fem.png`: hadron `S/E` and `C/E` versus
  truth crystal `f_em`, with fitted `h/e` lines.
- `plots/h_over_e_vs_energy.png`: energy dependence of independently fitted
  `(h/e)_S` and `(h/e)_C`.
- `plots/energy_response_distributions.png`: `S`, `C`, and dual-readout
  reconstructed energy distributions.
- `plots/reco_fem_vs_truth_fem.png`: reconstructed dual-readout `f_em` versus
  truth crystal `f_em`.
- `plots/leakage_fraction_vs_dr_residual.png`: leakage fraction versus
  dual-readout energy residual.

### Signal Definitions

The first-pass dual-readout analysis uses:

- `S_raw = Nph_Scintillation`
- `C_raw = Nph_Cherenkov`

Electron samples set the calibration constants:

```text
S_cal = aS * S_raw
C_cal = aC * C_raw
```

Hadron samples fit:

```text
S_cal / E = f_em + (h/e)_S * (1 - f_em)
C_cal / E = f_em + (h/e)_C * (1 - f_em)
```

Then the dual-readout reconstruction uses:

```text
chi = (1 - (h/e)_S) / (1 - (h/e)_C)
E_DR = (S_cal - chi * C_cal) / (1 - chi)
f_em_DR = (S_cal / E_DR - (h/e)_S) / (1 - (h/e)_S)
```

`EdepCrystal`, fiber/frame energy deposits, and `leakageEnergy` are diagnostics
and selections. They are not used as the primary dual-readout signals.
