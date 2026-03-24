# Normal Sleep-Wake Model Fitting

This directory contains scripts for fitting the three-compartment amyloid clearance model to normal sleep-wake cycle data from multiple experimental studies.

## Overview

The normal sleep-wake fitting approach analyzes how amyloid-beta concentrations change during natural sleep-wake cycles, as opposed to experimental interventions like lumbar puncture. This provides insights into the physiological regulation of amyloid clearance during normal conditions.

## Directory Structure

```
Normal_sleep_wake_fitting/
├── Global_fit/                    # Global optimization across all datasets
│   ├── main.m                     # Main global fitting script
│   ├── config.m               # Data loading
│   ├── euler.m            # Numerical integration
│   ├── nonlinear_constraints.m    # Non linear constraints
│   └── README.md                 # Detailed documentation
├── Individual_fits/               # Individual dataset fitting
│   ├── Blattner/                 # Blattner et al., 2020 fitting
│   │   ├── main.m        # Individual Blattner fit script
|   |   ├── config.m        
|   |   ├── euler.m        
|   |   ├── nonlinear_constraints.m
│   │   └── readme.md             
│   ├── Lucey/                 # Lucey et al., 2018 fitting
│   │   ├── main.m        # Individual Lucey fit script
|   |   ├── config.m        
|   |   ├── euler.m        
|   |   ├── nonlinear_constraints.m
│   │   └── readme.md             
│   ├── Liu/                 # Liu et al., 2023 fitting
│   │   ├── main.m        # Individual Liu fit script
|   |   ├── config.m        
|   |   ├── euler.m        
|   |   ├── nonlinear_constraints.m
│   │   └── readme.md             
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

### Individual Fitting (`Individual_fits/`)

**Purpose**: Optimize model parameters for each dataset separately.

## Model Structure

Both approaches use the same three-compartment model:

```
Brain (B) ←→ CSF (C) ←→ Plasma (P)
```

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
- `global_fit1/` directory with:
  - Individual run results (`run_X_results.mat`)
  - Final summary (`final_results.mat`)
  - Parameter space data (`parameter_space_X.csv`)
  - Visualization plots

### Individual Fitting Output
Each individual fit generates:
- `{individual_fit_name}/` directory with similar structure
- Dataset-specific parameter estimates
- Parameter space exploration plots

## Parameter Bounds

All fitting approaches use the same parameter bounds:

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
