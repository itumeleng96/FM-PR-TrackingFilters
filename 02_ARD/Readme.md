# 02_ARD — Amplitude Range Doppler map

Cross-correlation of the reference and surveillance signals to produce the
Amplitude Range Doppler (ARD) map:

```
|Psi(tau, v)| = |sum_{n=0..N-1} s_E(n) * s_R*(n - tau) * exp(j 2 pi v n / N)|
```

**Status:** ARD computation currently lives inline in `Runners/` scripts
(`run_filter_live.m`, `runEvaluateFilters.m`, ...). A future refactor
should extract it into a standalone `compute_ard.m` here so the pipeline
step becomes reusable.

Inputs:
- Reference channel HDF5 (`direct.h5`)
- Surveillance channel HDF5 (`echo.h5`)

Outputs:
- ARD map matrix, indexed by range bin and Doppler bin
