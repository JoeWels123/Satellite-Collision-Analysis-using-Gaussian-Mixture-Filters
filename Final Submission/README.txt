# Gaussian Sum Filter (GSF) — Orbital State Estimation (MATLAB)

A MATLAB toolbox for orbit determination using a Gaussian Sum Filter (GSF) with equinoctial elements, plus utilities for measurement modeling, conversions, and experiments.

---

## Table of Contents
- [Overview](#overview)
- [Function Index (APIs)](#function-index-apis)
- [Script Index (Experiments / Workflows)](#script-index-experiments--workflows)
- [Typical Workflow](#typical-workflow)
- [Conventions](#conventions)
- [Reference](#reference)
- [License](#license)

---

## Overview
The toolkit implements and propagates a Gaussian-Sum Unscented Kalman Filter (GS-UKF) using equinoctial orbital elements. Its main purposes are (1) to demonstrate how incorporating additional Gaussian components enhances filter performance when compared against ground-truth Monte Carlo simulations, and (2) to perform collision probability analyses between two real-world orbital objects.

---

## Function Index (APIs)

| File | Purpose | **Inputs** | **Outputs** |
|---|---|---|---|
| bhattacharyya_distance.m | Computes the Bhattacharyya distance between two multivariate Gaussian distributions. | mu1, Sigma1 (mean & cov of first Gaussian), mu2, Sigma2 (mean & cov of second Gaussian) | d (Bhattacharyya distance scalar) |

| calculate_mixture_moments.m | Calculate skewness and kurtosis for Gaussian mixture distribution over time. | gaussian_means (3×T×N), gaussian_covs (3×3×T×N), weights (1×N), numTimeSteps | mixture_skewness (3×T), mixture_kurtosis (3×T) |

| classical_to_equinoctial.m | Converts classical orbital elements to equinoctial coordinates with uncertainty. | elements1 (struct with meanmotion, Eccentricity, Inclination, RightAscensionOfAscendingNode, ArgumentOfPeriapsis, MeanAnomaly) | equinoctial (struct with a, h, k, p, q, l), P (6×6 covariance matrix) |

| compute_collision_probability_over_time_weighted.m | Compute collision probability between two objects over time using encounter frame transformation and 2D integration. | means1_pos, covs1_pos, vels1 (object 1 states), means2_pos, covs2_pos, vels2 (object 2 states), R (combined radius) | log_collision_probs (1×T), collision_probs (1×T), ratios (1×T), sigma_x_ratios (1×T), sigma_z_ratios (1×T) |

| equinoctial_to_classical.m | Converts equinoctial orbital coordinates back to classical Keplerian orbital elements. | sigma_point (6×1 vector: a, h, k, p, q, l) | elements1 (struct with eccentricity, inclination, raan, argp, meanAnomaly, meanmotion) |

| Generate_GSUKF.m | Initialize and propagate a Gaussian Sum Unscented Kalman Filter over time using equinoctial elements. | mu (equinoctial mean 6×1), P (covariance 6×6), N (num components), numSamples, startTime, stopTime, sampleTime, File (XML filename) | ISS_mixture_cov (6×6×T), ISS_mixture_mean (6×T), ISS_gaussian_means (6×T×N), ISS_gaussian_covs (6×6×T×N), ISS_mu_components (6×N), ISS_P_components (6×6×N), gaussian_weights (1×N), Wm_all (13×N), Wc_all (13×N) |

| Generate_Monte_Carlo.m | Generate and propagate Monte Carlo particle samples for orbital state estimation. | mu (equinoctial mean 6×1), P (covariance 6×6), startTime, stopTime, sampleTime, numSamples, File (XML filename) | next_Particle_ECI (6×T×numSamples), Weights (numSamples×1) |

| generate_sigma_points.m | Generate sigma points for the Unscented Transform using the scaled unscented transformation. | mu (mean 6×1), P (covariance 6×6), alpha, beta, kappa (UT parameters) | sigma_points (n×13), Wm (mean weights 1×13), Wc (covariance weights 1×13) |

| propagate_gaussian_mixture.m | Compute collision probability between two Gaussian mixture distributions over time. | fengyun_gaussian_means (6×T×N), cosmos_gaussian_means (6×T×N), fengyun_gaussian_covariances (6×6×T×N), cosmos_gaussian_covariances (6×6×T×N), fengyun_gaussian_weights (1×N), cosmos_gaussian_weights (1×N), numTimeSteps, N (num components) | collisionprobs_total (1×T), lognoprob (1×T), noprob (1×T) |

| updateOMMFile.m | Update Keplerian orbital elements in an OMM XML file and recalculate derived parameters. | inputFile (XML path), outputFile (XML path), elements (struct with meanMotion, eccentricity, inclination, raan, argPeriapsis, meanAnomaly, optional epoch) | None (modifies XML file) |

| COSMOS 1570.xml | Example of one of the two conjunction orbital objects, as listed on Space-Track.org |

| CZ4C DEB.xml | Example of the second conjunction orbital objects, as listed on Space-Track.org |

(To find orbital objects that will collide, navigate to the conjunctions part of the Space-Track.org website, this will give a list of colliding objects with their respective collision probabilities and other statistics)
---

## Script Index (Experiments / Workflows)

| Script | Purpose | **Primary Inputs** | **Primary Outputs** |
|---|---|---|---|
| Comparing_GSUKF_MonteCarlo.m | Benchmark script comparing GS-UKF performance against Monte Carlo ground truth across multiple Gaussian component counts (N=[1,20,40,60,80]). Evaluates convergence, accuracy, computational cost, and higher-order moment capture. | Scenario/time settings (startTime, stopTime, sampleTime), orbital element file (Object1), initial mu, P; numSamples for GMM fitting, num_Monte_Carlo_samples, array of N_values to test | Monte Carlo reference statistics (mean, covariance, skewness, kurtosis over time); for each N: GS-UKF mixture statistics, Bhattacharyya distances, Euclidean distances, computation times; comprehensive comparison figures (convergence vs N, distance evolution, moment matching, covariance volume, timing analysis) |

| Final_Predictions.m | End-to-end collision probability analysis between two space objects using GS-UKF propagation and Gaussian mixture collision assessment. Propagates both objects, finds temporal overlap, interpolates to common timebase, and computes pairwise Gaussian collision probabilities. | Scenario/time settings, orbital element files (Object1, Object2), initial conditions via classical_to_equinoctial; N (num Gaussian components), numSamples for GMM, combined hard-body radius R | Object1 and Object2: mixture means/covariances over time, Gaussian component means/covariances/weights; interpolated states on common timebase; collision probability time series (collisionprobs_total), relative distance statistics; figures showing position magnitudes with uncertainty, 3D Gaussian component visualization, relative distance evolution, collision probability over time with risk thresholds |

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

