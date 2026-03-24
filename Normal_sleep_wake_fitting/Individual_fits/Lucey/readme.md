# Individual Model Fitting: Lucey et al., 2018

This directory contains scripts for fitting the three-compartment amyloid clearance model specifically to the Lucey et al., 2018 dataset.

## Overview

The `lucey_fit.m` script performs individual model fitting to CSF concentration data from Lucey et al., 2018. This approach allows for dataset-specific parameter optimization, which can be compared with global fitting results to understand dataset-specific characteristics.

## Dataset

- **Study**: Lucey et al., 2018
- **Data Type**: CSF concentration measurements
- **Data File**: `lucey_wake_conc.csv` (note: relative path, not in data/ subdirectory)
- **Format**: CSV with columns: `Time`, `Concentration`, `LSD`, `USD`

## Script: `lucey_fit.m`

### Purpose
Performs optimization of the three-compartment model parameters using only the Lucey dataset.

### Key Features
- **Single Dataset Fitting**: Optimizes parameters using only Lucey CSF data
- **Multiple Optimization Runs**: Performs 20 optimization runs from different starting points
- **Parameter Space Visualization**: Generates plots showing parameter sensitivity
- **NRMSE Objective Function**: Uses normalized RMSE for error calculation

### Model Parameters
The script optimizes 10 parameters with the same bounds as global fitting:

| Parameter | Description | Units | Bounds |
|-----------|-------------|-------|--------|
| `r_bc` | Brain-to-CSF transfer rate (wake) | /hr | [0.01, 0.25] |
| `r_bp` | Brain-to-plasma transfer rate (wake) | /hr | [0.01, 0.25] |
| `r_cp` | CSF-to-plasma transfer rate (wake) | /hr | [0, 0.1] |
| `sigma_bc` | Sleep scaling for brain-to-CSF | - | [1, 7] |
| `sigma_bp` | Sleep scaling for brain-to-plasma | - | [1, 7] |
| `sigma_cp` | Sleep scaling for CSF-to-plasma | - | [1, 7] |
| `sigma_p` | Sleep scaling for plasma clearance | - | [1, 7] |
| `A` | Amyloid production rate (wake) | pg/ml/hr | [0, 111] |
| `sigma_A` | Sleep scaling for production | - | [0, 0.99] |
| `r_p` | Plasma clearance rate (wake) | /hr | [0, 0.6] |

## Usage

### Prerequisites
1. **Data File**: Ensure `lucey_wake_conc.csv` is present in the current directory
2. **MATLAB Requirements**: Optimization Toolbox for `fmincon`

### Running the Fit
```matlab
cd Normal_sleep_wake_fitting/Individual_fits/Lucey
lucey_fit
```

## Output

### Files Generated
- **`global_all_params/`** directory containing:
  - `run_X_results.mat`: Individual optimization results for each run
  - `parameter_space_X.csv`: Parameter space exploration data for each parameter
  - Parameter space visualization plots (saved as high-resolution PNG)

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

## Data Selection

- Uses indices [1:9, 14:20] from the 36-hour simulation window
- This selection pattern differs slightly from other individual fits
- Extracts CSF concentration data from the 100-day simulation

## Comparison with Global Fitting

Individual Lucey fitting results can be compared with global fitting to:
- Identify dataset-specific parameter characteristics
- Assess the robustness of global parameter estimates
- Understand how the Lucey dataset influences global fits

## Notes

- Simulates 100 days to reach steady state before extracting experimental window
- Sleep/wake cycles: 16-hour wake (8:00-24:00), 8-hour sleep (0:00-8:00)
- Initial conditions: Brain=0, CSF=600, Plasma=15.5 pg/ml
- Time step: 0.01 hours
- **Important**: Data file path is relative to the script directory (not in data/ subdirectory)
