# Tracking Filters for a Passive Radar System

## Author
**Itumeleng Malemela**
Master's Dissertation Repository

---

## Overview

This repository contains the MATLAB implementation of tracking filters designed
for an **FM Passive Radar System**, together with the full processing chain
(FERS simulation, ECA clutter cancellation, CA-CFAR, mean-shift clustering,
multi-target tracking, and Monte Carlo evaluation).

It also hosts the source of the CIE 2026 conference paper *"Design and
Implementation of Tracking Filters for a Passive Radar System"* (paper ID
1421).

---

## Repository layout

| Folder | Contents |
| --- | --- |
| `TrackingFilter-KalmanFilter/` | Baseline Kalman filter class |
| `TrackingFilter-UKF/` | Unscented Kalman filter class |
| `TrackingFilter-ParticleFilter/` | Particle filter class |
| `TrackingFilter-RGNF/` | Recursive Gauss-Newton filter class |
| `TrackingFilter-CSKF/` `TrackingFilter-CSUKF/` `TrackingFilter-CS-ParticleFilter/` `TrackingFilter-CSRGNF/` | Adaptive covariance-scaling variants of the four filters |
| `TrackingFilter-IMM/` | Interacting Multiple Model (future work) |
| `FERS/` | FERS scenario XMLs and HDF5 loaders |
| `cfar/` | CA-CFAR implementations |
| `meanShiftCluster/` | Mean-shift clustering for CFAR centroids |
| `DPI_Suppression/` | Direct-Path-Interference / ECA cancellation |
| `multiTargetTracking/` | MTT framework (gating, GNN, M-of-N) |
| `groundTruthCalculations/` | Cubic-waypoint ground truth generation |
| `Waveform/` | FM waveform synthesis utilities |
| `mc/` | Monte Carlo evaluation drivers |
| `tuning/` | Filter parameter sweeps and optimisation |
| `cache/` | Pre-computed CFAR/cluster caches per seed |
| `diagnostics/` | One-off checks and per-seed traces |
| `runners/` | Live-plot single-run scripts |
| `evaluationScripts/` | Legacy per-metric evaluation scripts |
| `TuningFilters/` | Legacy tuning outputs |
| `Simulation_results/` | Static images from prior campaigns |
| `figures/` | Figures used in the CIE 2026 paper |
| `seeds/` | *(gitignored)* per-seed FERS thermal-noise realisations |
| `docs_notes/` | Archived paper snippets kept for comparison |
| `paper/` | CIE 2026 conference paper source (.tex, .bib, figures, compiled PDF) |
| `UCT Msc Docs/` | Dissertation-related documents |

---

## Contents by purpose

### Monte Carlo evaluation (`mc/`)

| Script | Purpose |
| --- | --- |
| `mc_evaluate_filters.m` | 100 FERS-seeded MC across three scenarios × four filters. Produces Table I of the paper. |
| `mc_kf_ukf_orbit.m` | Targeted 100-seed re-run of KF+UKF on the 360° orbit with retuned parameters. Splices into `mc_results.mat`. |
| `mtt_mc.m` | 100 FERS-seeded MTT MC across three P_fa values × four filters. Produces Table IV. |
| `mc_noisy_evaluate_filters.m` | Injected-noise variant of the filter-perf MC (early draft). |
| `mc_winsorize_summary.m` | Winsorised summary of MC results to tame outliers. |
| `multi_scenario_eval.m` | Legacy multi-scenario evaluation. |
| `per_scenario_best.m` | Reports best per-scenario tuning across `mc_results.mat`. |

### Tuning sweeps (`tuning/`)

| Script | Purpose |
| --- | --- |
| `tune_pf_quick.m` | 3-seed × orbit PF tuning sweep. |
| `tune_kf_ukf_orbit.m` | 3-seed KF+UKF orbit sweep (higher-Q candidates). |
| `tune_kf_orbit_v2.m` | KF orbit sweep: shrink Q / tighten gate. |
| `tune_kf_orbit_v3.m` | KF orbit sweep: vary α, per-axis R_s, wider gate. |
| `optimize_all_filters.m` | Full grid-search over all four filters and three scenarios. |
| `optimize_cskf.m` `optimize_cskf_akhlaghi.m` `optimize_cskf_akhlaghi_v2.m` | CSKF-specific optimisations. |
| `scheme_b_grid_sweep.m` `scheme_b2_grid_sweep.m` | Adaptive-scheme grid sweeps. |
| `run_tuned_all_scenarios.m` | Re-run the tuned filters across all scenarios. |
| `test_adaptive_q.m` `test_cskf_per_axis_q.m` | Unit-style tests for the adaptive Q update. |

### Cache generation (`cache/`)

| Script | Purpose |
| --- | --- |
| `cache_all_clusters.m` | Generate cluster centroids for the four base scenarios. |
| `cache_3targets_clusters.m` | Generate cluster centroids for the three-target MTT scenario. |
| `cache_seeded_clusters.m` | Per-seed centroid cache for single-target scenarios. |
| `cache_seeded_3targets.m` | Per-seed centroid cache for the three-target MTT scenario, sweepable over CFAR P_fa. |

### Diagnostics (`diagnostics/`)

| Script | Purpose |
| --- | --- |
| `mtt_diag_one_seed.m` | Per-seed trace showing which track IDs are born/deleted. |
| `sanity_check_covscaling.m` | Verifies R/Q covariance-scaling update mechanics. |
| `run_cskf_diagnostics.m` | CSKF-internal state dumps for debugging. |
| `check_measurement_quality.m` | Signal/detection quality checks per scan. |

### Live runners and plotting (`runners/`)

| Script | Purpose |
| --- | --- |
| `run_filter_live.m` | Live-plotting single-filter run on a chosen scenario. |
| `run_cskf_live.m` | Live-plotting CSKF run. |
| `run_imm_live.m` | Live-plotting IMM run. |
| `run_single_filter.m` `run_single_filter_ard.m` `run_single_filter_scenario.m` | One-shot single-filter runs. |
| `single_run_cspf.m` | CSPF single-run helper. |
| `runEvaluateFilters.m` `runEvaluateFiltersPmatrix.m` | Legacy filter evaluation with per-run plots. |
| `runTrackingFilterPlot.m` | Legacy live tracking plot. |
| `plot_cluster_vs_truth_ard.m` | ARD montage with mean-shift centroids overlaid. |
| `show_meanshift_live.m` | Live-updating ARD with clusters. |
| `measurementData.m` | Legacy data preparation helper. |

### Utilities (repo root)

| Script | Purpose |
| --- | --- |
| `startup.m` | Auto-runs when MATLAB opens the project — adds all subfolders to path. |
| `gated_association.m` | Mahalanobis gate + nearest-neighbour association used across MC scripts. |

---

## Implemented Tracking Filters

The following tracking filters and associated algorithms are implemented:

1. **Kalman Filter**
2. **Unscented Kalman Filter**
3. **Particle Filter**
4. **Recursive Gauss-Newton Filter**
5. **Adaptive Covariance Scaling** (Akhlaghi 2018-style innovation/residual-based R and Q updates), applied to all four filters (with a particle-cloud analogue in the PF).

---

## Other algorithms and features

- **Cell-Averaging Constant False Alarm Rate (CA-CFAR) detector**

  ![CFAR Results](/Simulation_results/CFAR/cacfar_pfa-8.svg)

- **Multi-Target Tracker** (M-of-N with GNN association)

  ![MTT Results](/Simulation_results/Picture4.png)

- **Mean-Shift clustering**

  ![CFAR Centroid Results](/Simulation_results/CFAR/cacfar_centroids.png)

- **Extensive Cancellation Algorithm (ECA)** for DPI/clutter suppression

- **Log-Likelihood as the primary statistical benchmark**, corroborated by RMSE and NIS

  ![Log Likelihood](/Simulation_results/FERS_scenarios/360_range_ll_1.svg)

- **Cubic-waypoint ground truth for FERS scenarios**

  ![GT](/Simulation_results/FERS_scenarios/3D_360.svg)
  ![GT](/Simulation_results/FERS_scenarios/rangeDoppler360.svg)

---

## How to run

The `startup.m` at the repo root adds every subfolder to the MATLAB path, so
you can invoke any script by name from the repo root:

```bash
matlab -batch "mc_evaluate_filters(100)"     # Table I MC (100 seeds)
matlab -batch "mtt_mc(100)"                  # Table IV MTT MC
matlab -batch "tune_kf_ukf_orbit"            # KF/UKF orbit tuning sweep
matlab -batch "run_filter_live"              # Live single-filter plot
```

If MATLAB is opened interactively in the project root, `startup.m` runs
automatically and no extra `addpath` is needed.

---

## Companion paper

The CIE 2026 conference paper lives in `paper/`. Its Table I, Table IV, and
computational-load tables are produced by scripts in `mc/` and `tuning/`.

Build:
```bash
cd paper
pdflatex -shell-escape IEEE_ConferencePaperCIE.tex
bibtex IEEE_ConferencePaperCIE
pdflatex -shell-escape IEEE_ConferencePaperCIE.tex
pdflatex -shell-escape IEEE_ConferencePaperCIE.tex
```
