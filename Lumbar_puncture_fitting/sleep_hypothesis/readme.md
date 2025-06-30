
# Sleep Hypothesis Testing

This directory contains systematic testing of sleep-related hypotheses for amyloid clearance changes observed in lumbar puncture studies.

## Overview

The sleep hypothesis posits that lumbar puncture interventions primarily affect sleep-related parameters in the three-compartment amyloid clearance model. This approach tests which sleep scaling parameters are most important for explaining the experimental observations.

## Hypothesis Rationale

Sleep is known to significantly affect amyloid-beta clearance through:
- Changes in CSF flow dynamics
- Altered brain metabolic activity
- Modified glymphatic system function
- Changes in blood-brain barrier permeability

Lumbar puncture may affect sleep patterns or sleep quality, leading to changes in these sleep-dependent processes.

## Directory Structure

```
sleep_hypothesis/
├── all_sigma/                     # All sleep scaling parameters
├── sA/                           # Sleep scaling for production only
├── sbc_scp/                      # Sleep scaling for brain-CSF and CSF-plasma
├── sbp/                          # Sleep scaling for brain-plasma only
├── sp/                           # Sleep scaling for plasma clearance only
├── sA_sp/                        # Production + plasma clearance scaling
├── sbp_sp/                       # Brain-plasma + plasma clearance scaling
├── sbp_sA/                       # Brain-plasma + production scaling
├── sbp_sA_sp/                    # Brain-plasma + production + plasma clearance
├── sbc_scp_sp/                   # Brain-CSF + CSF-plasma + plasma clearance
├── sbc_scp_sbp_sp/               # All transfer scaling + plasma clearance
├── sbc_scp_sbp_sA/               # All transfer scaling + production
├── sbc_scp_sA_sp/                # Brain-CSF + CSF-plasma + production + plasma clearance
├── sbc_scp_sA/                   # Brain-CSF + CSF-plasma + production
├── sbc_sbp_scp/                  # All transfer scaling (no production/clearance)
└── run_all.sh                    # Batch processing script
```

## Parameters Tested

### Sleep Scaling Parameters

| Parameter | Description | Biological Meaning |
|-----------|-------------|-------------------|
| `sigma_A` | Sleep scaling for amyloid production | Sleep affects amyloid synthesis |
| `sigma_bc` | Sleep scaling for brain-to-CSF transfer | Sleep affects brain-CSF clearance |
| `sigma_bp` | Sleep scaling for brain-to-plasma transfer | Sleep affects brain-plasma clearance |
| `sigma_cp` | Sleep scaling for CSF-to-plasma transfer | Sleep affects CSF-plasma exchange |
| `sigma_p` | Sleep scaling for plasma clearance | Sleep affects systemic clearance |

### Parameter Combinations

Each subdirectory tests different combinations of these parameters:

- **`all_sigma/`**: Tests all 5 sleep scaling parameters simultaneously
- **`sA/`**: Tests only production scaling (most conservative)
- **`sbc_scp/`**: Tests brain-CSF and CSF-plasma transfer scaling
- **`sbp/`**: Tests brain-plasma transfer scaling only
- **`sp/`**: Tests plasma clearance scaling only
- **Combined directories**: Test various combinations of 2-4 parameters

## Usage

### Prerequisites

1. **Data Files**: Ensure all CSV files are present in `../data/` directory
2. **MATLAB**: Required for model fitting
3. **Bash**: Required for batch processing

### Running Individual Tests

Navigate to a specific parameter combination:

```bash
cd Lumbar_puncture_fitting/sleep_hypothesis/all_sigma
# Run the fitting script for all sleep scaling parameters
```

### Batch Processing

Run all sleep hypothesis tests:

```bash
cd Lumbar_puncture_fitting/sleep_hypothesis
./run_all.sh
```

This will process all parameter combinations systematically.

## Model Structure

All tests use the same three-compartment model with sleep/wake modulation:

```
Brain (B) ←→ CSF (C) ←→ Plasma (P)
```

**Sleep/Wake Cycle**:
- Wake: 16 hours (8:00-24:00)
- Sleep: 8 hours (0:00-8:00)

**Parameter Modulation**:
- Wake parameters: Base values
- Sleep parameters: Base values × sleep scaling factors

## Error Metric

Uses **Weighted RMSE (WRMSE)** to account for experimental uncertainties:

\[
\text{WRMSE} = \sqrt{\frac{\sum_{i=1}^{n}w_{i}.\left( y_{i} - \hat{y}_i\right)^2}{\sum_{i=1}^{n}w_{i}}}
\]

Where weights are inversely proportional to experimental standard errors.

## Analysis Approach

### Hypothesis Ranking

Compare different parameter combinations by:

1. **Overall Fit Quality**: Lower WRMSE indicates better fit
2. **Parameter Significance**: Which parameters show consistent changes
3. **Biological Plausibility**: Are the changes consistent with known sleep effects?
4. **Dataset Consistency**: Do different studies support the same conclusions?

### Key Questions

1. **Which sleep parameters are most affected by lumbar puncture?**
   - Production scaling (`sigma_A`)
   - Transfer rate scaling (`sigma_bc`, `sigma_bp`, `sigma_cp`)
   - Clearance scaling (`sigma_p`)

2. **Are the changes consistent across datasets?**
   - Blattner et al., 2020
   - Lucey et al., 2018
   - Liu et al., 2022

3. **What are the biological implications?**
   - Sleep quality effects
   - CSF flow changes
   - Metabolic alterations

## Expected Results

### Best-Fitting Hypotheses

Based on biological knowledge, we expect:

1. **`sigma_bc` and `sigma_cp`**: Most likely to be affected (CSF flow changes)
2. **`sigma_A`**: May be affected (sleep affects metabolism)
3. **`sigma_p`**: Less likely to be affected (systemic clearance)

### Parameter Patterns

- **Consistent changes**: Parameters that show similar effects across datasets
- **Dataset-specific effects**: Parameters that vary between studies
- **Interaction effects**: How different parameters work together

## Output Structure

Each test generates:

- **Parameter estimates**: Optimized values for tested parameters
- **Fit quality**: WRMSE and other error metrics
- **Model predictions**: Simulated concentration profiles
- **Visualization**: Plots comparing model to experimental data

## Comparison with Other Hypotheses

Compare sleep hypothesis results with:

- **Pressure hypothesis**: Are sleep effects more important than pressure effects?
- **Combined hypothesis**: Do sleep and pressure effects interact?
- **Normal sleep-wake fitting**: How do lumbar puncture effects compare to normal sleep?

## Notes

- Each parameter combination is tested with multiple optimization runs
- Results are compared systematically across all combinations
- Biological interpretation guides hypothesis selection
- Statistical significance is assessed through error metrics and parameter consistency 
