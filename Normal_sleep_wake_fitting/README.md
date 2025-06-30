# Normal Sleep-Wake Model Fitting

This directory contains scripts for fitting the three-compartment amyloid clearance model to normal sleep-wake cycle data from multiple experimental studies.

## Overview

The normal sleep-wake fitting approach analyzes how amyloid-beta concentrations change during natural sleep-wake cycles, as opposed to experimental interventions like lumbar puncture. This provides insights into the physiological regulation of amyloid clearance during normal conditions.

## Directory Structure

```
Normal_sleep_wake_fitting/
├── Global_fit/                    # Global optimization across all datasets
│   ├── main.m                     # Main global fitting script
│   ├── model.m                    # Three-compartment ODE model
│   ├── objective_function.m       # NRMSE objective function
│   ├── load_data.m               # Data loading and processing
│   ├── euler_solver.m            # Numerical integration
│   ├── plot_parameter_space.m    # Parameter space visualization
│   ├── save_results.m            # Results saving and start point generation
│   └── README.md                 # Detailed documentation
├── Individual_fits/               # Individual dataset fitting
│   ├── Blattner/                 # Blattner et al., 2020 fitting
│   │   ├── blattner_fit.m        # Individual Blattner fit script
│   │   └── readme.md             # Blattner-specific documentation
│   ├── Liu/                      # Liu et al., 2022 fitting
│   │   ├── liu_fit.m             # Individual Liu fit script (CSF + Plasma)
│   │   └── readme.md             # Liu-specific documentation
│   ├── Lucey/                    # Lucey et al., 2018 fitting
│   │   ├── lucey_fit.m           # Individual Lucey fit script
│   │   └── readme.md             # Lucey-specific documentation
│   └── readme.md                 # Individual fits overview
└── README.md                     # This file
```

## Datasets

The fitting uses data from three key studies:

1. **Blattner et al., 2020**: CSF concentration measurements during normal sleep-wake cycles
2. **Lucey et al., 2018**: CSF concentration measurements during normal sleep-wake cycles  
3. **Liu et al., 2022**: Both CSF and plasma concentration measurements during normal sleep-wake cycles

## Fitting Approaches

### Global Fitting (`Global_fit/`)

**Purpose**: Optimize model parameters by fitting to all datasets simultaneously.

**Advantages**:
- More robust parameter estimates by leveraging multiple studies
- Reduces overfitting to individual dataset characteristics
- Provides a unified model for normal sleep-wake dynamics

**Method**: Uses weighted NRMSE across all datasets to find optimal parameters.

### Individual Fitting (`Individual_fits/`)

**Purpose**: Optimize model parameters for each dataset separately.

**Advantages**:
- Reveals dataset-specific parameter characteristics
- Allows comparison of parameter estimates across studies
- Helps identify potential dataset-specific biases

**Method**: Uses NRMSE for each individual dataset.

## Model Structure

Both approaches use the same three-compartment model:

```
Brain (B) ←→ CSF (C) ←→ Plasma (P)
```

**Key Features**:
- Sleep/wake cycle modulation (16h wake, 8h sleep)
- 10 parameters optimized: transfer rates and sleep scaling factors
- 100-day simulation to reach steady state
- 36-hour experimental window extraction

## Usage

### Prerequisites

1. **Data Files**: Ensure all required CSV files are present:
   - `data/blattner_wake_conc.csv`
   - `data/lucey_wake_conc.csv` 
   - `data/liu_csf_wake_conc.csv`
   - `data/liu_plasma_wake_conc.csv`

2. **MATLAB Requirements**:
   - Optimization Toolbox (for `fmincon`)
   - Statistics and Machine Learning Toolbox

### Running Global Fitting

```matlab
cd Normal_sleep_wake_fitting/Global_fit
main
```

### Running Individual Fits

```matlab
% Blattner fit
cd Normal_sleep_wake_fitting/Individual_fits/Blattner
blattner_fit

% Liu fit (CSF + Plasma)
cd Normal_sleep_wake_fitting/Individual_fits/Liu
liu_fit

% Lucey fit
cd Normal_sleep_wake_fitting/Individual_fits/Lucey
lucey_fit
```

## Output Structure

### Global Fitting Output
- `global_all_params/` directory with:
  - Individual run results (`run_X_results.mat`)
  - Final summary (`final_results.mat`)
  - Parameter space data (`parameter_space_X.csv`)
  - Visualization plots

### Individual Fitting Output
Each individual fit generates:
- `global_all_params/` directory with similar structure
- Dataset-specific parameter estimates
- Parameter space exploration plots

## Error Metrics

### Normalized RMSE (NRMSE)
Used for all normal sleep-wake fitting:
\[
\text{NRMSE} = \frac{\text{RMSE}}{\text{range}(y)}
\]

This metric normalizes errors by the data range, making comparisons across datasets meaningful.

## Parameter Bounds

All fitting approaches use the same parameter bounds:

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

## Comparison and Analysis

### Global vs Individual Results
Compare global and individual fitting results to:
- Assess parameter consistency across datasets
- Identify dataset-specific characteristics
- Validate the global fitting approach

### Key Insights
- **Liu Dataset**: Unique dual-compartment fitting provides additional constraints
- **Parameter Robustness**: Global fitting typically provides more stable estimates
- **Dataset Compatibility**: Individual fits reveal how well datasets align

## Notes

- All scripts use the same random seed (42) for reproducibility
- 20 optimization runs from different starting points for robustness
- Sleep/wake cycles: 16-hour wake (8:00-24:00), 8-hour sleep (0:00-8:00)
- Initial conditions: Brain=0, CSF=600, Plasma=15.5 pg/ml
- Time step: 0.01 hours 
