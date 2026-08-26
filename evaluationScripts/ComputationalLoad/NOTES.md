# Computational Load — Notes

## Environment
- Precision 7550, i9-10885H, 32 GB, Gentoo, kernel 6.18, MATLAB R2025a
- FERS binary: `../../FERS/build/src/fers` (rebuilt against boost 1.90)
- CFAR `P_fa = 1e-5`, 4 guard / 6 training cells; mean-shift 10 km / 5 Hz
- MTT: confirmation=4, deletion=4, gating=20, M/N=3/5

## Scripts

| Script | Purpose | Output |
|---|---|---|
| `runComputationalLoadOriginal.m` | Original paper script | ratios only |
| `runComputationalLoadReproduction.m` | Reproduces original methodology on this box | ratios |
| `runComputationalLoadRawFilter.m` | **Table II source**: raw filter, 1000 MC, warmup=10 | `computational_load_raw_filter_results.mat` |
| `runComputationalLoadScaling.m` | **Table III source**: MTT-inclusive, N∈{1,3,5,10} | `computational_load_scaling_results.mat` |
| `summarizeScaling.m` | Computes Table III cells (mean ± σ) | stdout |

Run all scripts from the project root with:
```
matlab -batch "addpath('evaluationScripts/ComputationalLoad'); <script_name>"
```

## Correct MTT call order (per `evaluationScripts/runTrackingFilter.m`)
```
mtt = mtt.createNewTracks(centroids, i);
mtt = mtt.maintainTracks();
mtt = mtt.predictionStage();
mtt = mtt.updateStage(centroids, i);
```

## Numbers in the paper

### Table II (raw filter, μs mean±σ)
```
KF   | update  29.1±11.0  | predict  22.0±7.9
UKF  | update 121.6±44.0  | predict 146.6±57.6
PF   | update 1619 ±570   | predict 1703 ±649
RGNF | update  51.0±20.8  | predict  32.5±12.7
```

### Table III (MTT-inclusive per-scan, μs mean±σ)
```
Filter | N=1        | N=3        | N=5        | N=10
KF     |  945±72    | 1023±92    | 1566±107   |  2376±153
UKF    | 1383±57    | 1591±67    | 2527±80    |  3867±199
PF     | 9585±462   |11696±515   |19014±649   | 29858±1697
RGNF   |  808±35    |  928±66    | 1424±73    |  2207±124
```

Mean confirmed tracks: N=1 → 1.0, N=3 → 2.7, N=5 → 5.0, N=10 → 8.4.
Scaling ratios cost(N=10)/cost(N=1): KF 2.52×, UKF 2.80×, PF 3.11×, RGNF 2.73×.
Sub-linear vs 8.4× track ratio → fixed MTT overhead dominates.
