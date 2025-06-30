
# Combined Hypothesis Testing

This directory contains systematic testing of combined hypotheses that incorporate both sleep-related and pressure-related parameter changes for amyloid clearance alterations observed in lumbar puncture studies.

## Overview

The combined hypothesis posits that lumbar puncture interventions affect both sleep patterns and pressure dynamics simultaneously, leading to changes in both sleep scaling parameters and transfer rates. This approach tests the most comprehensive models that account for multiple physiological mechanisms.

## Directory Structure

```
combined_hypothesis/
├── rbc_rcp_all_sigma/             # Transfer rates + all sleep scaling
├── rbc_rcp_rbp_sbc_scp_sbp/       # All transfer rates + sleep scaling
├── rbc_rcp_rbp_sbc_scp_sbp_sA/    # All parameters (most comprehensive)
└── run_all.sh                     # Batch processing script
```

## Parameters Tested

### Combined Parameter Sets

Each subdirectory tests different combinations of transfer rates and sleep scaling parameters:

| Parameter Type | Parameters | Description |
|----------------|------------|-------------|
| **Transfer Rates** | `r_bc`, `r_bp`, `r_cp`, `r_p` | Pressure-dependent transfer rates |
| **Sleep Scaling** | `sigma_A`, `sigma_bc`, `sigma_bp`, `sigma_cp`, `sigma_p` | Sleep-dependent scaling factors |

### Parameter Combinations

- **`rbc_rcp_all_sigma/`**: 
  - Transfer rates: `r_bc`, `r_cp` (CSF pathway)
  - Sleep scaling: All 5 sleep parameters
  - Tests CSF pathway + all sleep effects

- **`rbc_rcp_rbp_sbc_scp_sbp/`**:
  - Transfer rates: `r_bc`, `r_cp`, `r_bp` (all brain transfers)
  - Sleep scaling: `sigma_bc`, `sigma_cp`, `sigma_bp` (transfer scaling)
  - Tests all brain transfers + transfer sleep scaling

- **`rbc_rcp_rbp_sbc_scp_sbp_sA/`**:
  - Transfer rates: `r_bc`, `r_cp`, `r_bp` (all brain transfers)
  - Sleep scaling: `sigma_bc`, `sigma_cp`, `sigma_bp`, `sigma_A` (transfers + production)
  - Most comprehensive test (all brain parameters)

## Usage

### Prerequisites

1. **Data Files**: Ensure all CSV files are present in `../data/` directory
2. **MATLAB**: Required for model fitting
3. **Bash**: Required for batch processing

### Running Individual Tests

Navigate to a specific parameter combination:

```bash
cd Lumbar_puncture_fitting/combined_hypothesis/rbc_rcp_all_sigma
# Run the fitting script for transfer rates + all sleep scaling
```

### Batch Processing

Run all combined hypothesis tests:

```bash
cd Lumbar_puncture_fitting/combined_hypothesis
./run_all.sh
```

This will process all parameter combinations systematically.

## Model Structure

All tests use the same three-compartment model with both pressure and sleep effects:

```
Brain (B) ←→ CSF (C) ←→ Plasma (P)
```

## Expected Results

### Best-Fitting Hypotheses

Based on biological knowledge, we expect:

1. **`rbc_rcp_all_sigma`**: Good fit (CSF pathway + sleep effects)
2. **`rbc_rcp_rbp_sbc_scp_sbp`**: Better fit (all brain transfers + sleep scaling)
3. **`rbc_rcp_rbp_sbc_scp_sbp_sA`**: Best fit (most comprehensive)

### Model Comparison

- **Akaike Information Criterion (AIC)**: Balance fit quality and complexity
- **Bayesian Information Criterion (BIC)**: Penalize complex models more heavily
