
# Combined Hypothesis Testing

This directory contains systematic testing of combined hypotheses that incorporate both sleep-related and pressure-related parameter changes for amyloid clearance alterations observed in lumbar puncture studies.

## Overview

The combined hypothesis posits that lumbar puncture interventions affect both sleep patterns and pressure dynamics simultaneously, leading to changes in both sleep scaling parameters and transfer rates. This approach tests the most comprehensive models that account for multiple physiological mechanisms.

## Hypothesis Rationale

Lumbar puncture may affect amyloid clearance through multiple mechanisms:

### Sleep Effects
- Discomfort or pain affecting sleep quality
- Changes in sleep patterns due to the procedure
- Altered sleep-dependent clearance processes

### Pressure Effects
- Direct changes in CSF pressure and flow
- Modified pressure gradients between compartments
- Altered blood-brain barrier dynamics

### Combined Effects
- Sleep and pressure effects may interact
- Multiple mechanisms may work simultaneously
- Comprehensive models may provide better fits

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

**Combined Effects**:
- Transfer rates modified by pressure changes
- Sleep scaling factors applied to sleep periods
- Both effects may interact and amplify each other

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
3. **Biological Plausibility**: Are the changes consistent with known mechanisms?
4. **Dataset Consistency**: Do different studies support the same conclusions?
5. **Model Complexity**: Balance between fit quality and parameter parsimony

### Key Questions

1. **Which combined effects are most important?**
   - Sleep + pressure interactions
   - Multiple transfer rate changes
   - Comprehensive parameter modifications

2. **Do combined models provide better fits than single hypotheses?**
   - Comparison with sleep-only models
   - Comparison with pressure-only models
   - Assessment of model improvement

3. **What are the biological implications?**
   - Multiple mechanism interactions
   - Synergistic effects
   - Comprehensive physiological understanding

## Expected Results

### Best-Fitting Hypotheses

Based on biological knowledge, we expect:

1. **`rbc_rcp_all_sigma`**: Good fit (CSF pathway + sleep effects)
2. **`rbc_rcp_rbp_sbc_scp_sbp`**: Better fit (all brain transfers + sleep scaling)
3. **`rbc_rcp_rbp_sbc_scp_sbp_sA`**: Best fit (most comprehensive)

### Parameter Patterns

- **Transfer rate changes**: Pressure effects on clearance pathways
- **Sleep scaling changes**: Sleep effects on clearance efficiency
- **Interaction effects**: How pressure and sleep effects combine
- **Synergistic effects**: Amplification of individual effects

## Biological Mechanisms

### Combined Effects

- **Sleep-Pressure Interactions**: Sleep quality may affect pressure recovery
- **Multiple Pathways**: Both CSF and blood-brain barrier pathways affected
- **Temporal Dynamics**: Pressure effects immediate, sleep effects delayed
- **Synergistic Mechanisms**: Combined effects may exceed sum of individual effects

### Physiological Implications

- **CSF Dynamics**: Pressure changes + sleep-dependent flow alterations
- **Blood-Brain Barrier**: Pressure effects + sleep-dependent permeability
- **Metabolic Changes**: Sleep effects on amyloid production + pressure effects on clearance
- **Systemic Effects**: Combined effects on plasma clearance

## Output Structure

Each test generates:

- **Parameter estimates**: Optimized values for all tested parameters
- **Fit quality**: WRMSE and other error metrics
- **Model predictions**: Simulated concentration profiles
- **Visualization**: Plots comparing model to experimental data
- **Comparison metrics**: Relative performance vs. single hypotheses

## Comparison with Other Hypotheses

Compare combined hypothesis results with:

- **Sleep hypothesis**: Do combined models improve on sleep-only fits?
- **Pressure hypothesis**: Do combined models improve on pressure-only fits?
- **Normal sleep-wake fitting**: How do lumbar puncture effects compare to normal conditions?

## Model Selection

### Criteria for Best Model

1. **Fit Quality**: Lowest WRMSE across datasets
2. **Parameter Significance**: Consistent parameter changes
3. **Biological Plausibility**: Reasonable parameter values
4. **Model Parsimony**: Balance between complexity and fit quality
5. **Prediction Accuracy**: Good fit to validation data

### Model Comparison

- **Akaike Information Criterion (AIC)**: Balance fit quality and complexity
- **Bayesian Information Criterion (BIC)**: Penalize complex models more heavily
- **Cross-validation**: Test model robustness
- **Parameter uncertainty**: Assess parameter confidence intervals

## Notes

- Each parameter combination is tested with multiple optimization runs
- Results are compared systematically across all combinations
- Biological interpretation guides hypothesis selection
- Statistical significance is assessed through error metrics and parameter consistency
- Combined models may provide the most realistic representation of lumbar puncture effects 
