# Gaussian Sum Filter (GSF) — Orbital State Estimation (MATLAB)

A MATLAB toolbox for orbit determination using a Gaussian Sum Filter (GSF) with equinoctial elements, plus utilities for measurement modeling, conversions, and experiments.

---

## Table of Contents
- [Overview](#overview)
- [Repo Map](#repo-map)
- [Function Index (APIs)](#function-index-apis)
- [Script Index (Experiments / Workflows)](#script-index-experiments--workflows)
- [Typical Workflow](#typical-workflow)
- [Conventions](#conventions)
- [Reference](#reference)
- [License](#license)

---

## Overview
The toolkit initializes and propagates a GSF over equinoctial orbital elements, performs measurement updates (AER from a ground station), and compares different estimation strategies (pure GSF, GSF+UKF, Monte Carlo baselines, different Gaussian splits, and collision-probability–oriented perturbations).

---

## Repo Map


---

## Function Index (APIs)

| File | Purpose | **Inputs** | **Outputs** |
|---|---|---|---|
| `GenerateGSF.m` | Initialize & propagate a Gaussian Sum Filter over time. | `mu` (equinoctial mean), `P` (cov), `N` (num components), `numSamples`, `startTime`, `stopTime` | `ISS_mixture_mean`, `ISS_mixture_cov`, `ISS_gaussian_means`, `ISS_gaussian_covs`, `gaussian_weights`, `ISS_sigma_points_cartesian`, `ISS_all_Wm`, `ISS_all_Wc` |
| `UpdateGSF.m` | Measurement update for each Gaussian component and the mixture. | GSF state arrays (`ISS_mixture_*`, `ISS_gaussian_*`, `gaussian_weights`), `mu`, `P`, `N`, `numSamples`, `numTimeSteps`, `lat`, `lon`, `sc_mean`, `startTime`, `stopTime`, `timeVector`, `sampleTime` | `updated_mixture_mean`, `updated_mixture_cov`, `updated_gaussian_mean`, `updated_gaussian_cov`, `update_times`, `covariance_traces`, `update_indices`, `measurement_residuals` |
| `UpdateGSFWeights.m` | Bayes update of mixture weights from innovations. | `w_prior`, `y` (cell of innovations), `S` (cell of innovation covariances) | `w_updated` (normalized) |
| `MakeEstimates.m` | Single-component (linearized) KF update step. | `mu_prior`, `P_prior`, `z` (meas), `R`, `z_hat` (pred meas), `H` (Jacobian) | `mu_est`, `P_est`, `S` (innovation cov), `y` (innovation) |
| `jacobianECItoAER.m` | Jacobian of AER wrt ECI position. | `targetECI`, `observerECI`, `lat`, `lon`, `gst` | `J` (3×3) |
| `groundStationECI.m` | Ground station geodetic → ECI. | `lat_deg`, `lon_deg`, `alt`, `gst_deg` | `r_eci` |
| `convertOMMToEquinoctialElements.m` | OMM/TLE → equinoctial + initial covariance. | `ommData` (struct) | `nominal` (equinoctial), `uncertainty` (struct), `P` (cov) |
| `generate_sigma_points.m` | Unscented transform sigma points & weights. | `mu`, `P`, `alpha`, `beta`, `kappa` | `sigma_points`, `Wm`, `Wc` |
| `equinoctialToKeplerian.m` | Equinoctial → Keplerian elements. | `equinoctial` | `keplerian` |
| `convertKeplerianToEquinoctialElements.m` | Keplerian → Equinoctial elements. | `keplerian` | `equinoctial` |

> **Note:** Variable prefixes like `ISS_` indicate outputs produced along a time history for the ISS demo case.

---

## Script Index (Experiments / Workflows)

| Script | Purpose | **Primary Inputs** | **Primary Outputs** |
|---|---|---|---|
| `Updating_only_GSF_using_Measurement3.m` | **Latest** experiment: update **only GSF weights** from measurements (no Gaussian collapse), to adjust the mixture without state merging. | Scenario/time settings, ground-station location, measurement stream; initial `mu`, `P`, `N` | Updated mixture **weights**, updated mixture statistics over time; diagnostics (innovations, cov traces, timing); figures/logs as configured |
| `Updating_only_GSF_using_Measurement2.m` | Complete experiment: update weights **and** collapse/resample Gaussians at each update. | Same as above | Updated mixture **means/covs** with periodic collapse, resampled components/weights; diagnostics and plots |
| `Updating_both_GSF_and_UKF_using_Measurement2.m` | Joint experiment: update **UKF sigma-point weights** and **GSF weights**. | Same as above plus UKF UT params (`alpha`,`beta`,`kappa`) | Updated UKF/GSF weights, mixture stats, sigma-point logs; diagnostics and plots |
| `ISS_Positions_Monte_Carlo.m` | Monte Carlo truth model: propagate ISS positions to obtain reference mean/cov for comparisons. | Scenario initialization, # of samples, times | Ground-truth **mean/cov trajectories**, position ensembles; comparison plots |
| `Comparing_Different_Gaussian_Splits.m` | Ablation: how initial GSF splitting across equinoctial components affects accuracy. (Empirically: per-element splits work best.) | Initial `mu`, `P`, split strategy configs, `N` | Performance metrics across splits, updated mixture stats, plots/tables |
| `Collision_prob_mixture.m` | Explore perturbations minimizing collision probability via final **a** and **e** tuning (following literature approach). | Initial orbit, perturbation parameterization, time horizon | Best perturbation (semi-major axis, eccentricity) minimizing collision probability; probability curves, plots |

> Scripts typically configure scenario parameters internally. If a script saves MAT files/figures, the destination is set within the script (adjust paths if needed).

---

## Typical Workflow

1. **Initialize scenario & uncertainty**
   - Parse OMM/TLE → `nominal`, `P` using `convertOMMToEquinoctialElements.m`.
   - Choose mixture size `N` and sigma-point params (if using UKF).

2. **Build & propagate the mixture**
   - Call/inspect `GenerateGSF.m` outputs over `[startTime, stopTime]`.

3. **Ingest measurements & update**
   - Use `UpdateGSF.m` (full update), or run one of the scripts:
     - *Weight-only* updates (no collapse): `Updating_only_GSF_using_Measurement3.m`
     - *Collapse+resample* updates: `Updating_only_GSF_using_Measurement2.m`
     - *UKF+GSF* weight updates: `Updating_both_GSF_and_UKF_using_Measurement2.m`

4. **Compare against ground truth**
   - Run `ISS_Positions_Monte_Carlo.m` and compare mixture mean/cov vs Monte Carlo.

5. **Research utilities**
   - Split strategy study: `Comparing_Different_Gaussian_Splits.m`
   - Collision probability minimization: `Collision_prob_mixture.m`

---

## Conventions
- **State**: Equinoctial elements for estimation; conversions provided to/from Keplerian.
- **Measurements**: AER (azimuth, elevation, range) from a ground station; use `groundStationECI.m` and `jacobianECItoAER.m`.
- **Weights**: `Wm`, `Wc` denote sigma-point mean/cov weights; `gaussian_weights` are mixture component weights.

---

## Reference
- **Collision probability method** as explored in `Collision_prob_mixture.m` draws on:  
  *“Methods to minimize collision probability …”* (Acta Astronautica), [ScienceDirect](https://www.sciencedirect.com/science/article/pii/S0273117721005937).

---

## License
Specify your license (e.g., MIT) here.
