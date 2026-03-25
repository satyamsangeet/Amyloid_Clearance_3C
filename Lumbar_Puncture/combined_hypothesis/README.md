# Combined hypothesis

Sub-hypotheses live as immediate subfolders of this directory (e.g. `rbc_rcp_all_sigma`, `rbc_rcp_sbc_scp`, … — the exact set is defined in `run_all.sh`). Each subfolder typically holds MATLAB plotting scripts (`global_plot.m`, `blattner_plot.m`, `lucey_plot.m`, `liu_plot.m`) and fit outputs.

## Scripts in this folder

| File | Purpose |
|------|---------|
| `run_all.sh` | Runs the full hypothesis-testing batch for **every** sub-hypothesis folder listed inside it |
| `hypothesis_replace_params1.py` | Replaces parameters in **model1** with **baseline / default** values in all relevant `.m` files |
| `hypothesis_replace_params2.py` | Replaces parameters in **model2** with the **best-fit** values for each sub-hypothesis |
| `run_plot_scripts.sh` | Runs all plotting `.m` scripts in each sub-hypothesis folder (needs MATLAB) |
| `copy_model_files_to_sim_files.py` | Copies generated `{fit}_model2_{subhypothesis}.csv` files into `sim_files/` and its subdirectories |
| `param_replace.py` | Lower-level helpers for updating parameter assignments in MATLAB source |

## Recommended order (run from **this** directory)

All commands below assume your shell’s current working directory is **`combined_hypothesis/`** (the folder that contains these scripts).

1. `./run_all.sh`  
2. `python hypothesis_replace_params1.py`  
3. `python hypothesis_replace_params2.py`  
4. `./run_plot_scripts.sh`  
5. `python copy_model_files_to_sim_files.py`  

After step 5, the `sim_files/` tree here is ready for downstream use.

## AICc

Run **`AICC.ipynb`** from the **parent** directory (`final_hypothesis_testing_with_updated_bounds/`), once `sim_files/` for this family (and any others you include) reflects the latest CSVs.

## Notes

- **model1** = baseline / default model in the plotting code; **model2** = hypothesis / best-fit model.
- Plotting scripts emit CSVs with compartment-wise model trajectories for **global** and **individual** fits; `copy_model_files_to_sim_files.py` stages **model2** CSVs into `sim_files/` per fit type.
- You do **not** need to run `param_replace.py` by hand for the standard pipeline; the `hypothesis_replace_params*.py` scripts drive the replacements.
