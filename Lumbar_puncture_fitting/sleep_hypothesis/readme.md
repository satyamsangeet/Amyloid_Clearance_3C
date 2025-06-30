# Sleep Hypothesis Testing

This directory contains systematic testing of sleep-related hypotheses for amyloid clearance changes observed in lumbar puncture studies.

## Overview

The sleep hypothesis posits that lumbar puncture interventions primarily affect sleep-related parameters in the three-compartment amyloid clearance model. This approach tests which sleep scaling parameters are most important for explaining the experimental observations.

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

## Output Structure

Each test generates:

- **Parameter estimates**: Optimized values for tested parameters
- **Fit quality**: NRMSE and other error metrics
- **Model predictions**: Simulated concentration profiles
- **Visualization**: Plots comparing model to experimental data
