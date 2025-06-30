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
