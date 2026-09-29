# Design and Implementation of Tracking Filters for a Passive Radar System

**Author:** Itumeleng Malemela (UCT / CSIR)
**Supervisors:** Stephen Paine, Francois Schonken
**Venue:** 2026 CIE International Conference on Radar — Xi'an, China, 16–18 October 2026

> Single source of truth for the **PowerPoint slides** and the **conference poster**. Contents below mirror the current 14-slide deck as built in `CIE Conference Slides - Tracking Filters for a FM PR.pptx`. Slides flagged **⚠ EMPTY** are title-only placeholders that still need content.

---

## Contents (current deck — 14 slides)

| # | Title | Status |
|---|---|---|
| 1  | Title — Tracking Filters for a FM Passive Radar | ✅ done |
| 2  | Background & Motivation (Why FM PR?) | ✅ done |
| 3  | System Overview | ⚠ figure only, needs bullets |
| 4  | RADAR Processing Chain — FERS → ARD | ✅ done |
| 5  | RADAR Processing Chain — ECA → CA-CFAR | ✅ done |
| 6  | RADAR Processing Chain — (mean-shift centroiding?) | ⚠ EMPTY — placeholder |
| 7  | Tracking Filters — KF, UKF, RGNF, PF | ✅ done |
| 8  | Multi-Target Tracking | ⚠ EMPTY — placeholder |
| 9  | Adaptive Filtering | ⚠ EMPTY — placeholder |
| 10 | Performance Evaluation | ⚠ EMPTY — placeholder |
| 11 | Log-Likelihood | ✅ done |
| 12 | NEES | ⚠ EMPTY — and not in the paper; see note |
| 13 | Computation Load | ✅ done |
| 14 | Recommendations & Future Work | ✅ done |

> **Note on Slide 12 (NEES):** The camera-ready paper benchmarks filters using **Log-Likelihood only**, not NEES. Options: (a) drop the slide, (b) replace with a scenario-specific LL results slide, or (c) add NEES analysis (requires re-running Monte Carlo).
>
> **Note on ordering:** "Adaptive Filtering" (Slide 9) currently sits *after* MTT (Slide 8), but adaptive covariance scaling is a property of the single-target filters. Consider swapping so it reads Filters → Adaptive → MTT.

---

## Slide 1 — Title

- **Title:** Tracking Filters for a FM Passive Radar
- **Author:** Itumeleng Malemela¹,²
- **Supervisors:** Stephen Paine¹, Francois Schonken¹
- **Affiliations:** ¹University of Cape Town · ²Council for Scientific and Industrial Research
- **Venue:** 2026 CIE International Conference on Radar — Xi'an, China, 16–18 October 2026
- **Figure:** UCT + CSIR logos

**Speaker note:** Open in one sentence — "tracking commercial aircraft in the bistatic range–Doppler domain using FM broadcast signals as illuminators of opportunity."

---

## Slide 2 — Background & Motivation

**Framing question:** Why track commercial aircraft in the bistatic range–Doppler domain using FM PR?

Using existing FM band illuminators of opportunity:
- **Cost effective** — no dedicated transmitter.
- **No spectrum licensing** required.
- **High Doppler resolution** from long integration times.
- **Covert** — no emissions from the receiver side.

**Speaker note:** Land the motivation for ATC at small airports in developing regions. Flag the trade-off up front: long integration buys Doppler resolution but you inherit FM's limited bandwidth → poor range resolution.

---

## Slide 3 — System Overview

**Figure:** Bistatic geometry — 1 Tx, 2 Rx (reference + surveillance), targets.

- **FERS** — Flexible Extensible Radar Simulator (UCT, Brooker 2008).
- Isotropic FM illuminator at **94 MHz**, Tx–Rx baseline **12.5 km**.
- Directional surveillance antenna; omnidirectional reference antenna.
- Simulation environment used throughout this work.

> ⚠ Slide currently has diagrams but no bullets — add the four bullets above.

---

## Slide 4 — RADAR Processing Chain: FERS → ARD

**Layout:** left = FERS scene (`fers_scene_scenario1.svg` — produced by `plotFersScene_scenario1.m` in this folder); right = 3-D ARD map (`ardMapFinal.svg`); below (centred, spanning) = discrete-time cross-ambiguity equation (`ambiguity_function_eqn.svg`).

**Discrete-time cross-ambiguity function** — computed from FERS reference/surveillance I/Q:

$$
|\Psi(\tau, v)| = \left| \sum_{n=0}^{N-1} s_E(n)\, s_R^{*}(n-\tau)\, e^{j 2\pi v n / N} \right|
$$

- $s_E(n)$: surveillance echo · $s_R^{*}(n-\tau)$: conjugate reference
- Cross-correlation processed at **200 kS/s**.
- **Left figure — FERS scene:** rendered from `scenario_1_singleFile.fersxml` (PrieskaTx at (6440, 10760, 1056) m; RefRx at origin, altitude 1000 m; SurvRx colocated ~1 m offset; single target, cubic-interpolated waypoints at 3600 m altitude, 60 s duration). Shows baseline ≈ 12.5 km, target trajectory, and the two signal paths $s_R$ (Tx → target) and $s_E$ (target → SurvRx).
- **Right figure — 3-D ARD map:** dominant peak at zero Doppler = **DPI leakage**; sidelobes from the FM ambiguity function mask target returns.
- The equation is the operator that maps the *left* (I/Q from the FERS scene) to the *right* (ARD surface).

**Speaker note:** Talk through the scene left → equation → ARD. Regenerate the figure by opening `plotFersScene_scenario1.m` in MATLAB and running it — it writes both SVG and 300-dpi PNG into `Figures/`. Tweak `view(-40, 22)` at the bottom of the script for a different camera angle.

**Mesh helpers (procedural, no third-party STL required):**
- `airlinerMesh.m` — parametric twin-engine A320-class airliner (fuselage + swept wings + tail + engines). Called at t=0 s waypoint.
- `broadcastTowerMesh.m` — parametric triangular-lattice mast with optional stacked dipole bays. Used for Tx (with dipoles) and both Rx (mast only, no dipoles).

Reason for procedural over downloaded STL: avoids licence hazards for the conference deck; keeps the presentation folder small and self-contained; both meshes are ~few-hundred triangles so `patch` renders in milliseconds.

---

## Slide 5 — RADAR Processing Chain: ECA → CA-CFAR

**DPI & clutter suppression:**
- **Extensive Cancellation Algorithm (ECA)** — LMS-based, applied *before* range-Doppler processing on the surveillance channel.
- Cardinali *et al.* — near-complete clutter suppression with ECA.

**CA-CFAR detection:**

$$
T = \alpha \cdot P_n, \qquad \alpha = N \left( P_{fa}^{-1/N} - 1 \right)
$$

- **4 guard cells**, **6 training cells** each side.
- Averaging in **Doppler dimension only** (target smears across range bins due to FM bandwidth).
- Operating point: **$P_{fa} = 10^{-5}$**.
- **Figure:** CA-CFAR threshold at 15 km bistatic range for $P_{fa}$ from $10^{-3}$ to $10^{-7}$ at $t = 2$ s.

---

## Slide 6 — RADAR Processing Chain: Centroiding (Mean-Shift)

> ⚠ **SLIDE CURRENTLY EMPTY.** Suggested content below — the SVG in the deck appears to be a mean-shift centroiding figure.

- CFAR detections spread across adjacent range/Doppler cells → cluster them.
- **Mean-shift** — no prior on cluster count, returns centroids directly.
- Bandwidths: **5 Hz** Doppler, **10 km** range.
- **Figure:** CFAR output with mean-shift centroids at $t = 2$ s.

**Speaker note:** Bridge into tracking: centroids are the measurements that feed the tracking filter update step.

---

## Slide 7 — Tracking Filters

Four filters evaluated:

| Filter | Role | Notes |
|---|---|---|
| **KF**   | Linear Gaussian baseline | Efficient, robust |
| **UKF**  | Unscented transform | Handles non-linear bistatic geometry, no Jacobians |
| **RGNF** | Recursive Gauss–Newton | Iterative refinement, adaptive to unknown dynamics |
| **PF**   | Sequential Monte Carlo | Robust to non-Gaussian noise; 10 000 particles |

**Target state (Nearly Constant Velocity):**

$$
x(t) = [\, r_t,\ \dot r_t,\ f_t,\ \dot f_t \,]^{T}
$$

**Speaker note:** NCV is appropriate for non-manoeuvring commercial aircraft on the CPI timescales in this work.

---

## Slide 8 — Multi-Target Tracking

> ⚠ **SLIDE CURRENTLY EMPTY.** Suggested content below.

- **Gating** — ellipsoidal Mahalanobis region: $Z_k S_k^{-1} Z_k^{T} \le G$.
- **Global Nearest Neighbour (GNN)** association:

$$
d_{G_{ij}}^{2} = d_{ij}^{2} + \ln|S_i|
$$

- **Track management:** M-out-of-N — confirm a track if a detection is seen in ≥ M of the last N scans.
- Unassigned detections → new tentative tracks.
- **Figure (optional):** MTT output at $t = 8$ s using KF + GNN + 3-out-of-4 M/N.

---

## Slide 9 — Adaptive Filtering

> ⚠ **SLIDE CURRENTLY EMPTY.** Suggested content below. Consider moving this slide *before* Slide 8 — adaptive scaling is a per-filter property, not an MTT property.

Innovation and residual:

$$
d_k = z_k - h_k(\hat X_k^{-}), \qquad \varepsilon_k = z_k - h_k(\hat X_k^{+})
$$

Online update for $R$:

$$
R_k = \alpha R_{k-1} + (1 - \alpha)\!\left( \varepsilon_k \varepsilon_k^{T} + H_k P_k^{-} H_k^{T} \right)
$$

- Forgetting factor $0 \le \alpha \le 1$.
- Inflates measurement noise when innovations blow up.
- Prevents divergence from **CFAR outliers** and **zero-Doppler cancellation artefacts**.
- Applied to all four filters.

---

## Slide 10 — Performance Evaluation

> ⚠ **SLIDE CURRENTLY EMPTY.** Suggested content below — this slide sets up the results section.

**Three simulated flight scenarios (FERS):**
- **Straight-and-level flight** — 60 s
- **Landing** — descent + turn, 120 s
- **360° orbit** — holding pattern, 120 s

**Metrics reported:**
- **Log-likelihood** (statistical consistency, Slide 11)
- **Computational load** (relative runtime, Slide 13)
- **Data association** — true / false track confirmations across $P_{fa}$

Adaptive covariance scaling applied on all filters throughout.

---

## Slide 11 — Log-Likelihood

$$
\ell_k = -\tfrac{1}{2}\!\left[\, \ln|S_k| + e_k^{T} S_k^{-1} e_k + n \ln(2\pi) \,\right]
$$

with innovation covariance $S_k = H_k P_k^{-} H_k^{T} + R_k$ and innovation $e_k = z_k - H_k \hat X_k^{-}$.

- Captures **both** estimation accuracy **and** filter uncertainty (unlike RMSE alone).
- As filter converges, $P \downarrow$ and $\ell_k \to 0$.
- Higher (less negative) = better statistical consistency with ground truth.
- **Figure:** Doppler LL vs. time for KF/UKF/PF/RGNF on the 360° scenario.

**Ranked LL summary (Table I in paper):**

| Scenario | Domain | KF | UKF | PF | RGNF |
|---|---|---|---|---|---|
| Landing  | Range   | 2 | 3 | 4 | **1** |
| Landing  | Doppler | 3 | **1** | 4 | 2 |
| 360°     | Range   | 3 | **1** | 4 | 2 |
| 360°     | Doppler | 3 | 2 | 4 | **1** |
| Take-off | Range   | 2 | 3 | 4 | **1** |
| Take-off | Doppler | 2 | **1** | 4 | 3 |

- **RGNF** consistently top — iterative Gauss–Newton update → rapid convergence.
- **UKF** — strong Doppler LL, but outlier-sensitive during manoeuvres.
- **PF** — last despite 10 000 particles.

---

## Slide 12 — NEES

> ⚠ **SLIDE CURRENTLY EMPTY — and NEES is not analysed in the camera-ready paper.**
>
> Choose one of:
> - **(a) Drop the slide.** LL already covers consistency; NEES adds a second axis with no supporting analysis.
> - **(b) Replace with a scenario-specific LL results slide** (e.g., 360° range LL alongside Slide 11's Doppler LL).
> - **(c) Add NEES analysis.** Would require re-running the Monte Carlo trials and adding NEES computation to the pipeline.

If keeping (c), the NEES statistic for a Gaussian filter is:

$$
\epsilon_k = (x_k - \hat x_k)^{T} P_k^{-1} (x_k - \hat x_k)
$$

with the $\chi^2$ consistency bounds at the chosen confidence level.

---

## Slide 13 — Computation Load

**Table II** — Relative runtime (KF = 1.00, averaged over 1 000 Monte Carlo runs of a 60 s single-target scenario):

| Filter | Update | Prediction |
|---|---|---|
| KF   | 1.00 | 1.00 |
| UKF  | 0.96 | 1.39 |
| PF   | **6.53** | **3.56** |
| RGNF | 1.12 | 1.08 |

- All meet real-time constraints for **single-target** tracking (µs per scan; CPI ≈ 0.5–1 s).
- In **MTT**, cost scales linearly with confirmed tracks → **PF becomes prohibitive**.
- **KF and RGNF** — best accuracy/compute balance for MTT deployment.

---

## Slide 14 — Recommendations & Future Work

**Recommendations:**
- **RGNF** — best LL across scenarios at compute cost ≈ KF baseline.
- **KF** — comparable stability at the lowest cost; good MTT default.
- **UKF** — strong Doppler tracking, but sensitive to outliers during manoeuvres.
- **Adaptive covariance scaling** was essential for robustness against CFAR outliers and zero-Doppler artefacts.

**Future Work:**
- **Validation** with real-world PR system data.
- **Deep learning** for complex target dynamics.
- **Multiple Hypothesis Tracking (MHT)** for improved data association.
- **Optimisation** of tracking filter parameter tuning.
- Cognitive waveform / illuminator selection.

**Contact:** mlmitu002@uct.ac.za · imalemela@csir.co.za

> **Questions?**

---

# Poster mapping

The poster is 1.2 × 0.9 m (Radar 2026 template). Suggested layout — three columns, header banner, footer with logos + contact:

**Header banner (full width)**
- Title (from Slide 1)
- Authors + affiliations
- CIE 2026 logo (right), UCT + CSIR logos (left)

**Column 1 — Motivation & System**
- Slide 2 (Why FM PR?) — condensed to 3 bullets
- Slide 3 (System overview + bistatic geometry figure)
- Short "Contribution" callout box (not currently a slide — pull from paper introduction)

**Column 2 — Processing Chain & Filters**
- Slide 4 (ARD equation + 3-D ARD map)
- Slide 5 (ECA + CA-CFAR equation + threshold figure)
- Slide 6 (Mean-shift centroiding + figure)
- Slide 7 (four-filter table + NCV state)
- Slide 9 (Adaptive covariance scaling equation)

**Column 3 — Results & Conclusions**
- 360° hero figure — bistatic geometry + filter outputs
- Slide 11 (Doppler LL plot + Table I — ranked LL)
- Slide 13 (Table II — compute)
- MTT output figure (3-target scenario) + data association tables (Tables III–V)
- Slide 14 (Recommendations + Future Work — 3 bullets each)

**Footer**
- QR code → paper PDF
- Contact emails
- Acknowledgements line

---

# Figure map (files in `Figures/`)

Primary figures used by the slides and poster:

| Paper ref | File in `Figures/` | Slide / Poster |
|---|---|---|
| **Fig. 1** — Radar system overview | `radarOverview.drawio.png` | Slide 3 / Poster Col 1 |
| **Fig. 2** — 3-D ARD map | `ardMapFinal.svg` | Slide 4 / Poster Col 2 |
| Cross-ambiguity equation (discrete) | `ambiguity_function_eqn.svg` | Slide 4 |
| FERS scene — scenario_1 (setup diagram) | `fers_scene_scenario1.svg` (from `plotFersScene_scenario1.m`) | Slide 4 |
| **Fig. 3** — CA-CFAR threshold vs $P_{fa}$ | `ca_cfar_validate_corrections_v2.svg` | Slide 5 / Poster Col 2 |
| **Fig. 4** — CFAR output + centroids | `cacfar_pfa7_centroid_correction.svg` | Slide 6 / Poster Col 2 |
| **Fig. 5(a)** — 3-D target path (360°) | `3D_360.svg` | Slide 11 / Poster Col 3 |
| **Fig. 5(b)** — Filter outputs (360°) | `360_filter_output.svg` | Slide 11 / Poster Col 3 |
| **Fig. 6** — Doppler LL (360°) | `360_doppler_ll.svg` | Slide 11 / Poster Col 3 |
| **Fig. 7** — MTT output (3 targets) | `mtt_3_tracking_corrections_v1.svg` | Slide 8 / Poster Col 3 |

Supplementary figures available in `Figures/` (backup slides / alternative poster panels):

- `fersBistaticRadar.drawio.png` — bistatic geometry block diagram (alt Fig. 5a)
- `360_range_ll.svg` — Range LL for the 360° scenario (candidate for Slide 12 replacement)
- `landing_filter_outputs.svg`, `landingRangeLL.svg`, `dopplerLandingLL.svg` — landing scenario
- `takeoffOutputs.svg`, `ll_range_climb.svg`, `ll_doppler_climb_corrections_v1.svg` — take-off scenario
- `ardBeforeDPI_unnormalized_correction1.svg`, `ardFinalAfter_correction1.svg` — DPI cancellation before/after
- `mtt_3_cfar_corrections_v1.svg` — CFAR panel for the 3-target MTT scenario

**PowerPoint tip:** SVGs paste directly (Insert → Pictures → This Device) and stay vector. If a specific SVG doesn't render cleanly, screenshot it in your browser, or convert with `inkscape file.svg --export-type=png --export-dpi=300`.
