# Three Compartment Model for Brain Amyloid Clearance

This repository contains code and data for simulating and analyzing a three-compartment model of amyloid-beta clearance in the brain. The compartments represent:

- **Brain (B)**
- **Cerebrospinal Fluid (CSF, C)**
- **Plasma (P)**

The model is used to fit and interpret experimental data from multiple studies, exploring the effects of sleep, pressure, and other physiological factors on amyloid clearance.

---

## Repository Structure

```
├── dmsd.py                       # Python script for combining means and SDs
├── model_dynamics.m              # MATLAB: Simulate and plot the 3-compartment model
├── parameter_variation_analysis.m # MATLAB: Analyze effect of parameter variation
├── Lumbar_puncture_fitting/      # Model fitting to lumbar puncture data
│   ├── data/                     # Experimental data (CSV)
│   ├── sleep_hypothesis/         # Model fits under sleep hypothesis
│   ├── combined_hypothesis/      # Model fits under combined hypothesis
│   ├── pressure_hypothesis/      # Model fits under pressure hypothesis
│   └── ...
├── Normal_sleep_wake_fitting/    # Model fitting to normal sleep-wake data
│   ├── Global_fit/               # Scripts for global model fitting
│   ├── Individual_fits/          # Individual fits for Blattner, Liu, Lucey datasets
│   └── ...
└── README.md                     # This file
```

---

## Model Description

The model simulates amyloid-beta production and clearance between the brain, CSF, and plasma, with parameters for transfer rates and their modulation by sleep/wake states. Key parameters include:

- **A**: Amyloid production during wake
- **sigma_A**: Scaling factor for production during sleep
- **r_bc, r_bp, r_cp, r_p**: Transfer rates between compartments during wake
- **sigma_bc, sigma_bp, sigma_cp, sigma_p**: Scaling factors for transfer rates during sleep

See the parameter tables below for optimized values from different fits.

---

## Main Scripts and Usage

### Python
- **dmsd.py**: Combines means and standard deviations from multiple samples. Run with:
  ```bash
  python3 dmsd.py
  ```

### MATLAB
- **model_dynamics.m**: Simulates the three-compartment model and plots amyloid concentrations over time.
- **parameter_variation_analysis.m**: Analyzes the effect of varying a parameter (e.g., sigma_p) on model output.

To run these scripts, open them in MATLAB and execute (e.g., `model_dynamics` or `parameter_variation_analysis`).

#### Fitting Scripts
- **Lumbar_puncture_fitting/**: Contains scripts and data for fitting the model to lumbar puncture datasets under different hypotheses. Use the provided `run_all.sh` scripts to automate batch runs (requires MATLAB and bash).
- **Normal_sleep_wake_fitting/Global_fit/**: Scripts for global fitting to normal sleep-wake data.
- **Normal_sleep_wake_fitting/Individual_fits/**: Scripts for fitting individual datasets (Blattner, Liu, Lucey).

---

## Data

Experimental data from published studies are provided in CSV format under `Normal_sleep_wake_fitting/data/` for normal sleep-wake fitting
- `blattner_wake_conc.csv`
- `lucey_wake_conc.csv`
- `liu_csf_wake_conc.csv`
- `liu_plasma_wake_conc.csv`

Experimental data from published studies are provided in CSV format under `Lumbar_puncture_fitting/data/` for lumbar puncture fitting
- `blattner2020_csf_concentration.csv`
- `lucey2018_csf_concentration.csv`
- `liu2022_csf_concentration.csv`
- `liu2022_plasma_concentration.csv`
---

## Error Metrics

### Normalized RMSE (NRMSE)
Used when all data points are treated equally (e.g., Huang et al., 2012). Formula:

\[
\text{NRMSE} = \frac{\text{RMSE}}{\text{range}(y)} \text{ or } \frac{\text{RMSE}}{\mu_y}
\]

### Weighted RMSE (WRMSE)
Used when data points have associated uncertainties (e.g., Liu et al., 2022). Formula:

\[
\text{WRMSE} = \sqrt{\frac{\sum_{i=1}^{n}w_{i}.\left( y_{i} - \hat{y}_i\right)^2}{\sum_{i=1}^{n}w_{i}}}
\]

---

## Optimized Model Parameters

### Three Compartment Model (Brain, CSF, Plasma)

| **Parameter** | **Liu Fit** | **Huang Fit** | **Combined Fit** |
|---------------|-------------|---------------|------------------|
| a             | 1.00        | 2.95          | 1.49             |
| b             | 1.60        | 4.39          | 1.49             |
| c             | 1.53        | 4.01          | 1.47             |
| A_wake        | 22.06       | 21.07         | 23.5             |
| A_sleep       | 0.8*A_wake  | 0.8*A_wake    | 0.8*A_wake       |
| a12_wake      | 1.00        | 0.12          | 0.99             |
| a12_sleep     | 2.5*a12_wake| 2.5*a12_wake  | 2.5*a12_wake     |
| a13_wake      | 0.10        | 0.01          | 0.09             |
| a13_sleep     | a*a13_wake  | a*a13_wake    | a*a13_wake       |
| a23_sleep     | 0.066       | 0.066         | 0.066            |
| a23_wake      | a23_sleep/b | a23_sleep/b   | a23_sleep/b      |
| k_wake        | 0.23        | 0.23          | 0.28             |
| k_sleep       | c*k_wake    | c*k_wake      | c*k_wake         |

---

## References
- Huang et al., 2012
- Liu et al., 2022
- Blattner et al., 2020
- Lucey et al., 2018

---

## Contact
For questions or contributions, please open an issue or contact the repository maintainer.
