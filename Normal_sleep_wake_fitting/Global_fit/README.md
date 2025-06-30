# Global Model Fitting for Three-Compartment Amyloid Clearance

This directory contains scripts for performing global optimization of the three-compartment amyloid clearance model using multiple experimental datasets simultaneously.

## Overview

The global fitting approach optimizes model parameters by fitting to multiple experimental datasets at once, including:
- Blattner et al., 2020 (CSF concentration data)
- Lucey et al., 2018 (CSF concentration data) 
- Liu et al., 2022 (CSF and plasma concentration data)

This approach provides more robust parameter estimates by leveraging information from multiple studies.

## Scripts

### Main Scripts

- **`main.m`**: Main entry point that orchestrates the global fitting process
- **`model.m`**: Implements the three-compartment ODE system
- **`objective_function.m`**: Calculates the objective function (NRMSE) for optimization
- **`load_data.m`**: Loads and processes experimental data from CSV files

### Supporting Scripts

- **`euler_solver.m`**: Numerical solver using Euler method for ODE integration
- **`plot_parameter_space.m`**: Visualizes parameter space exploration results
- **`save_results.m`**: Saves optimization results and generates start points
- **`output_function.m`**: Callback function for optimization progress tracking

## Usage

### Prerequisites

1. **Data Files**: Ensure the following CSV files are present in a `data/` subdirectory:
   - `blattner_wake_conc.csv`
   - `lucey_wake_conc.csv` 
   - `liu_csf_wake_conc.csv`
   - `liu_plasma_wake_conc.csv`

   Each CSV should contain columns: `Time`, `Concentration`, `LSD`, `USD`

2. **MATLAB Requirements**:
   - Optimization Toolbox (for `fmincon`)
   - Statistics and Machine Learning Toolbox (for data processing)

### Running the Global Fit

1. **Navigate to the directory**:
   ```matlab
   cd Normal_sleep_wake_fitting/Global_fit
   ```

2. **Run the main script**:
   ```matlab
   main
   ```

### Configuration

The optimization can be configured in `main.m`:

- **`num_starts`**: Number of optimization runs from different starting points (default: 20)
- **`initial_guess`**: Initial parameter guess
- **`lb`, `ub`**: Lower and upper bounds for parameters
- **Optimization options**: Algorithm, tolerances, maximum iterations

## Model Parameters

The model optimizes 10 parameters:

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

## Error Metric

The optimization uses **Normalized Root Mean Square Error (NRMSE)**:

\[
\text{NRMSE} = \frac{\text{RMSE}}{\text{range}(y)}
\]

This metric normalizes the error by the range of experimental data, making it suitable for comparing fits across different datasets with different scales.

## Output

### Files Generated

1. **`global_all_params/`** directory containing:
   - `run_X_results.mat`: Individual optimization results for each run
   - `final_results.mat`: Summary of all optimization runs
   - `parameter_space_X.csv`: Parameter space exploration data for each parameter

2. **Console Output**:
   - Progress of each optimization run
   - Best parameter values found
   - Final objective function value

### Results Structure

The `final_results.mat` file contains:
- `solutions`: Array of optimization results for each run
- `best_params`: Best parameter set found
- `best_fval`: Best objective function value
- `all_trajectories`: Optimization history for each run

## Visualization

The script automatically generates:
- Parameter space exploration plots showing the relationship between each parameter and the objective function
- CSV files for further analysis in other tools

## Notes

- The model simulates 100 days to reach steady state before extracting the 36-hour experimental window
- Sleep/wake cycles are modeled as 16-hour wake (8:00-24:00) and 8-hour sleep (0:00-8:00)
- Initial conditions: Brain=0, CSF=600, Plasma=15.5 pg/ml
- Time step for integration: 0.01 hours

## Troubleshooting

- **Data not found**: Ensure CSV files are in the correct `data/` subdirectory
- **Optimization fails**: Check parameter bounds and initial guesses
- **Memory issues**: Reduce `num_starts` for fewer optimization runs
