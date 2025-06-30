# Individual Model Fitting: Blattner et al., 2020

This directory contains scripts for fitting the three-compartment amyloid clearance model specifically to the Blattner et al., 2020 dataset.

## Overview

The `blattner_fit.m` script performs individual model fitting to CSF concentration data from Blattner et al., 2020. This approach allows for dataset-specific parameter optimization, which can be compared with global fitting results to understand dataset-specific characteristics.

## Dataset

- **Study**: Blattner et al., 2020
- **Data Type**: CSF concentration measurements
- **Data File**: `data/blattner_wake_conc.csv`
- **Format**: CSV with columns: `Time`, `Concentration`, `LSD`, `USD`

## Script: `blattner_fit.m`

### Purpose
Performs optimization of the three-compartment model parameters using only the Blattner dataset.

### Key Features
- **Single Dataset Fitting**: Optimizes parameters using only Blattner CSF data
- **Multiple Optimization Runs**: Performs 20 optimization runs from different starting points
- **Parameter Space Visualization**: Generates plots showing parameter sensitivity
- **NRMSE Objective Function**: Uses normalized RMSE for error calculation

### Model Parameters
The script optimizes 10 parameters with the same bounds as global fitting:

| Parameter | Description | Units | Bounds |
|-----------|-------------|-------|--------|
| `r_bc` | Brain-to-CSF transfer rate (wake) | /hr | [0.75, 2.5] |
| `r_bp` | Brain-to-plasma transfer rate (wake) | /hr | [0.01, 1] |
| `r_cp` | CSF-to-plasma transfer rate (wake) | /hr | [0.001, 3] |
| `sigma_bc` | Sleep scaling for brain-to-CSF | - | [1, 4] |
| `sigma_bp` | Sleep scaling for brain-to-plasma | - | [1, 7] |
| `sigma_cp` | Sleep scaling for CSF-to-plasma | - | [1, 7] |
| `sigma_p` | Sleep scaling for plasma clearance | - | [1, 7] |
| `A` | Amyloid production rate (wake) | pg/ml/hr | [9, 14] |
| `sigma_A` | Sleep scaling for production | - | [0.7, 0.9] |
| `r_p` | Plasma clearance rate (wake) | /hr | [0.23, 0.34] |

## Usage

### Prerequisites
1. **Data File**: Ensure `data/blattner_wake_conc.csv` is present
2. **MATLAB Requirements**: Optimization Toolbox for `fmincon`

### Running the Fit
```matlab
cd Normal_sleep_wake_fitting/Individual_fits/Blattner
blattner_fit
```

## Output

### Files Generated
- **`global_all_params/`** directory containing:
  - `run_X_results.mat`: Individual optimization results for each run
  - `parameter_space_X.csv`: Parameter space exploration data for each parameter
  - Parameter space visualization plots

### Console Output
- Progress of each optimization run
- Best parameter values found
- Final objective function value (NRMSE)

## Error Metric

Uses **Normalized Root Mean Square Error (NRMSE)**:
\[
\text{NRMSE} = \frac{\text{RMSE}}{\text{range}(y)}
\]

This normalizes the error by the range of experimental data, making it suitable for comparing fits across different scales.

## Comparison with Global Fitting

Individual fitting results can be compared with global fitting to:
- Identify dataset-specific parameter characteristics
- Assess the robustness of global parameter estimates
- Understand how individual datasets influence global fits

## Notes

- Simulates 100 days to reach steady state before extracting experimental window
- Sleep/wake cycles: 16-hour wake (8:00-24:00), 8-hour sleep (0:00-8:00)
- Initial conditions: Brain=0, CSF=600, Plasma=15.5 pg/ml
- Time step: 0.01 hours
