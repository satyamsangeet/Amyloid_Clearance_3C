# Individual Model Fitting: Liu et al., 2022

This directory contains scripts for fitting the three-compartment amyloid clearance model specifically to the Liu et al., 2022 dataset.

## Overview

The `liu_fit.m` script performs individual model fitting to both CSF and plasma concentration data from Liu et al., 2022. This is unique among individual fits as it uses both CSF and plasma measurements, providing a more comprehensive view of amyloid clearance dynamics.

## Dataset

- **Study**: Liu et al., 2022
- **Data Types**: CSF and plasma concentration measurements
- **Data Files**: 
  - `data/liu_csf_wake_conc.csv`
  - `data/liu_plasma_wake_conc.csv`
- **Format**: CSV with columns: `Time`, `Concentration`, `LSD`, `USD`

## Script: `liu_fit.m`

### Purpose
Performs optimization of the three-compartment model parameters using only the Liu dataset (both CSF and plasma data).

### Key Features
- **Dual Dataset Fitting**: Optimizes parameters using both Liu CSF and plasma data
- **Multiple Optimization Runs**: Performs 20 optimization runs from different starting points
- **Parameter Space Visualization**: Generates plots showing parameter sensitivity
- **Combined NRMSE Objective Function**: Uses average of CSF and plasma NRMSE

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
1. **Data Files**: Ensure both Liu CSV files are present in `data/` directory
2. **MATLAB Requirements**: Optimization Toolbox for `fmincon`

### Running the Fit
```matlab
cd Normal_sleep_wake_fitting/Individual_fits/Liu
liu_fit
```

## Output

### Files Generated
- **`global_all_params/`** directory containing:
  - `run_X_results.mat`: Individual optimization results for each run
  - `parameter_space_X.csv`: Parameter space exploration data for each parameter
  - Parameter space visualization plots

### Console Output
- Progress of each optimization run
- Individual errors for CSF and plasma fits
- Best parameter values found
- Final combined objective function value

## Error Metric

Uses **Combined Normalized Root Mean Square Error (NRMSE)**:
\[
\text{Total NRMSE} = \frac{\text{NRMSE}_{CSF} + \text{NRMSE}_{Plasma}}{2}
\]

Where each NRMSE is calculated as:
\[
\text{NRMSE} = \frac{\text{RMSE}}{\text{range}(y)}
\]

This approach gives equal weight to both CSF and plasma measurements.

## Unique Aspects

### Dual Compartment Fitting
Unlike other individual fits that use only CSF data, Liu fitting:
- Fits to both CSF and plasma measurements simultaneously
- Provides more constraints on the model parameters
- Allows validation of the three-compartment model structure

### Data Selection
- CSF data: Uses indices [1:7, 14:20] from the 36-hour simulation window
- Plasma data: Uses the same indices for consistency
- Both datasets are normalized independently before combining

## Comparison with Global Fitting

Individual Liu fitting results can be compared with global fitting to:
- Assess how dual-compartment data influences parameter estimates
- Validate the three-compartment model structure
- Understand the contribution of plasma measurements to global fits

## Notes

- Simulates 100 days to reach steady state before extracting experimental window
- Sleep/wake cycles: 16-hour wake (8:00-24:00), 8-hour sleep (0:00-8:00)
- Initial conditions: Brain=0, CSF=600, Plasma=15.5 pg/ml
- Time step: 0.01 hours
- Equal weighting of CSF and plasma errors in objective function
