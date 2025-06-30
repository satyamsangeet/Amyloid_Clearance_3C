# Pressure Hypothesis Testing

This directory contains systematic testing of pressure-related hypotheses for amyloid clearance changes observed in lumbar puncture studies.

## Overview

The pressure hypothesis posits that lumbar puncture interventions primarily affect transfer rates between compartments due to changes in CSF pressure and flow dynamics. This approach tests which transfer rate parameters are most important for explaining the experimental observations.

## Directory Structure

```
pressure_hypothesis/
├── rbc_rcp/                       # Brain-CSF and CSF-plasma transfer rates
├── rbc_rcp_rbp/                   # All brain transfer rates
├── rbp_rp/                        # Brain-plasma and plasma clearance rates
├── rbc_rcp_rp/                    # Brain-CSF, CSF-plasma, and plasma clearance
└── run_all.sh                     # Batch processing script
```

## Parameters Tested

### Transfer Rate Parameters

| Parameter | Description | Biological Meaning |
|-----------|-------------|-------------------|
| `r_bc` | Brain-to-CSF transfer rate | Amyloid clearance from brain to CSF |
| `r_bp` | Brain-to-plasma transfer rate | Amyloid clearance from brain to plasma |
| `r_cp` | CSF-to-plasma transfer rate | Amyloid exchange between CSF and plasma |
| `r_p` | Plasma clearance rate | Systemic amyloid clearance |

### Parameter Combinations

Each subdirectory tests different combinations of these parameters:

- **`rbc_rcp/`**: Tests brain-CSF and CSF-plasma transfer rates
  - Focuses on the main CSF clearance pathway
  - Most directly affected by CSF pressure changes

- **`rbc_rcp_rbp/`**: Tests all brain transfer rates
  - Includes both brain-CSF and brain-plasma pathways
  - Tests whether pressure affects all brain clearance routes

- **`rbp_rp/`**: Tests brain-plasma and plasma clearance rates
  - Focuses on systemic clearance pathways
  - Tests pressure effects on blood-brain barrier and systemic clearance

- **`rbc_rcp_rp/`**: Tests brain-CSF, CSF-plasma, and plasma clearance
  - Comprehensive test of all transfer rates
  - Most complete pressure hypothesis test

## Usage

### Prerequisites

1. **Data Files**: Ensure all CSV files are present in `../data/` directory
2. **MATLAB**: Required for model fitting
3. **Bash**: Required for batch processing

### Running Individual Tests

Navigate to a specific parameter combination:

```bash
cd Lumbar_puncture_fitting/pressure_hypothesis/rbc_rcp
# Run the fitting script for brain-CSF and CSF-plasma transfer rates
```

### Batch Processing

Run all pressure hypothesis tests:

```bash
cd Lumbar_puncture_fitting/pressure_hypothesis
./run_all.sh
```

This will process all parameter combinations systematically.

## Model Structure

All tests use the same three-compartment model with pressure-dependent transfer rates:

```
Brain (B) ←→ CSF (C) ←→ Plasma (P)
```

**Pressure Effects**:
- Transfer rates may be modified by pressure changes
- Sleep/wake cycle modulation is maintained
- Production and clearance parameters remain fixed

## Analysis Approach

### Hypothesis Ranking

Compare different parameter combinations by:

1. **Overall Fit Quality**: Lower WRMSE indicates better fit
2. **Parameter Significance**: Which transfer rates show consistent changes
3. **Biological Plausibility**: Are the changes consistent with pressure effects?
4. **Dataset Consistency**: Do different studies support the same conclusions?

## Expected Results

### Best-Fitting Hypotheses

Based on pressure dynamics, we expect:

1. **`r_bc` and `r_cp`**: Most likely to be affected (direct CSF pressure effects)
2. **`r_bp`**: May be affected (blood-brain barrier pressure effects)
3. **`r_p`**: Less likely to be affected (systemic clearance)

### Parameter Patterns

- **CSF pathway effects**: Changes in `r_bc` and `r_cp`
- **Blood-brain barrier effects**: Changes in `r_bp`
- **Systemic effects**: Changes in `r_p`
- **Interaction effects**: How different transfer rates work together

## Output Structure

Each test generates:

- **Parameter estimates**: Optimized values for tested transfer rates
- **Fit quality**: WRMSE and other error metrics
- **Model predictions**: Simulated concentration profiles
- **Visualization**: Plots comparing model to experimental data
