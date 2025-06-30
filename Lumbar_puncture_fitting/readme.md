
# Lumbar Puncture Model Fitting

This directory contains scripts for fitting the three-compartment amyloid clearance model to experimental data from lumbar puncture studies, testing different hypotheses about the mechanisms underlying amyloid clearance changes.

## Overview

Lumbar puncture studies provide experimental data on how amyloid-beta concentrations change in response to interventions that alter CSF dynamics. This directory contains systematic testing of different hypotheses about which parameters are most affected by these interventions.

## Directory Structure

```
Lumbar_puncture_fitting/
├── data/                          # Experimental data files
│   ├── blattner2020_csf_concentration.csv
│   ├── lucey2018_csf_concentration.csv
│   ├── liu2022_csf_concentration.csv
│   └── liu2022_plasma_concentration.csv
├── sleep_hypothesis/              # Testing sleep-related parameter changes
│   ├── all_sigma/                 # All sleep scaling parameters
│   ├── sA/                        # Sleep scaling for production only
│   ├── sbc_scp/                   # Sleep scaling for brain-CSF and CSF-plasma
│   ├── sbp/                       # Sleep scaling for brain-plasma only
│   ├── sp/                        # Sleep scaling for plasma clearance only
│   └── ...                        # Other parameter combinations
│   └── run_all.sh                 # Batch processing script
├── pressure_hypothesis/           # Testing pressure-related parameter changes
│   ├── rbc_rcp/                   # Brain-CSF and CSF-plasma transfer rates
│   ├── rbc_rcp_rbp/               # All brain transfer rates
│   ├── rbp_rp/                    # Brain-plasma and plasma clearance rates
│   └── run_all.sh                 # Batch processing script
├── combined_hypothesis/           # Testing combined parameter changes
│   ├── rbc_rcp_all_sigma/         # Transfer rates + all sleep scaling
│   ├── rbc_rcp_rbp_sbc_scp_sbp/   # All transfer rates + sleep scaling
│   └── run_all.sh                 # Batch processing script
├── sleep_hypothesis_plot.m        # Plotting script for sleep hypothesis results
├── pressure_hypothesis_plot.m     # Plotting script for pressure hypothesis results
├── combined_hypothesis_plot.m     # Plotting script for combined hypothesis results
└── README.md                      # This file
```

## Experimental Context

Lumbar puncture studies measure amyloid-beta concentrations in CSF and plasma before and after interventions that:
- Alter CSF pressure
- Modify sleep patterns
- Change CSF flow dynamics

These studies help identify which physiological mechanisms are most important for amyloid clearance.

## Hypothesis Testing Approach

### Sleep Hypothesis (`sleep_hypothesis/`)

**Rationale**: Sleep is known to affect amyloid clearance. This hypothesis tests which sleep-related parameters are most important.

**Parameters Tested**:
- `sigma_A`: Sleep scaling for amyloid production
- `sigma_bc`: Sleep scaling for brain-to-CSF transfer
- `sigma_bp`: Sleep scaling for brain-to-plasma transfer
- `sigma_cp`: Sleep scaling for CSF-to-plasma transfer
- `sigma_p`: Sleep scaling for plasma clearance

**Subdirectories**:
- `all_sigma/`: All sleep scaling parameters
- `sA/`: Production scaling only
- `sbc_scp/`: Brain-CSF and CSF-plasma scaling
- `sbp/`: Brain-plasma scaling only
- `sp/`: Plasma clearance scaling only

### Pressure Hypothesis (`pressure_hypothesis/`)

**Rationale**: Lumbar puncture affects CSF pressure, which may alter transfer rates between compartments.

**Parameters Tested**:
- `r_bc`: Brain-to-CSF transfer rate
- `r_bp`: Brain-to-plasma transfer rate
- `r_cp`: CSF-to-plasma transfer rate
- `r_p`: Plasma clearance rate

**Subdirectories**:
- `rbc_rcp/`: Brain-CSF and CSF-plasma transfer rates
- `rbc_rcp_rbp/`: All brain transfer rates
- `rbp_rp/`: Brain-plasma and plasma clearance rates

### Combined Hypothesis (`combined_hypothesis/`)

**Rationale**: Lumbar puncture may affect both sleep patterns and pressure dynamics simultaneously.

**Parameters Tested**: Combinations of transfer rates and sleep scaling parameters.

**Subdirectories**:
- `rbc_rcp_all_sigma/`: Transfer rates + all sleep scaling
- `rbc_rcp_rbp_sbc_scp_sbp/`: All transfer rates + sleep scaling

## Datasets

The fitting uses data from three studies:

1. **Blattner et al., 2020**: CSF concentration measurements
2. **Lucey et al., 2018**: CSF concentration measurements
3. **Liu et al., 2022**: Both CSF and plasma concentration measurements

## Usage

### Prerequisites

1. **Data Files**: Ensure all CSV files are present in `data/` directory
2. **MATLAB**: Required for model fitting
3. **Bash**: Required for batch processing scripts

### Running Individual Hypothesis Tests

Each subdirectory contains its own fitting scripts. Navigate to the desired hypothesis and parameter combination:

```bash
cd Lumbar_puncture_fitting/sleep_hypothesis/all_sigma
# Run the fitting script for this parameter combination
```

### Batch Processing

Use the `run_all.sh` scripts to process all parameter combinations within a hypothesis:

```bash
# Run all sleep hypothesis tests
cd Lumbar_puncture_fitting/sleep_hypothesis
./run_all.sh

# Run all pressure hypothesis tests
cd Lumbar_puncture_fitting/pressure_hypothesis
./run_all.sh

# Run all combined hypothesis tests
cd Lumbar_puncture_fitting/combined_hypothesis
./run_all.sh
```

### Visualization

Use the plotting scripts to visualize results:

```matlab
% Sleep hypothesis results
cd Lumbar_puncture_fitting
sleep_hypothesis_plot

% Pressure hypothesis results
pressure_hypothesis_plot

% Combined hypothesis results
combined_hypothesis_plot
```

## Model Structure

All hypothesis tests use the same three-compartment model:

```
Brain (B) ←→ CSF (C) ←→ Plasma (P)
```

**Key Features**:
- Sleep/wake cycle modulation
- Parameter optimization for specific hypotheses
- Comparison with experimental data
- Error metric: Weighted RMSE (WRMSE) for data with uncertainties

## Error Metrics

### Weighted RMSE (WRMSE)
Used for lumbar puncture fitting due to experimental uncertainties:

\[
\text{WRMSE} = \sqrt{\frac{\sum_{i=1}^{n}w_{i}.\left( y_{i} - \hat{y}_i\right)^2}{\sum_{i=1}^{n}w_{i}}}
\]

Where weights are inversely proportional to experimental uncertainties.

## Output Structure

Each hypothesis test generates:
- Parameter optimization results
- Model fits to experimental data
- Error metrics and statistics
- Visualization plots

## Analysis Approach

### Hypothesis Ranking
Compare different hypotheses by:
- Overall fit quality (WRMSE)
- Parameter significance
- Biological plausibility
- Consistency across datasets

### Key Questions
1. Which parameters are most affected by lumbar puncture?
2. Are sleep-related or pressure-related changes more important?
3. How do different datasets support different hypotheses?
4. What are the implications for amyloid clearance mechanisms?

## Notes

- Each hypothesis test uses multiple optimization runs for robustness
- Results are compared across different parameter combinations
- Biological interpretation guides hypothesis selection
- Statistical significance is assessed through error metrics 
