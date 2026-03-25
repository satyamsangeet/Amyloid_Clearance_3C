# Pressure hypothesis

Sub-hypotheses are immediate subfolders of this directory (`rbc`, `rcp`, `rbc_rcp`, `rbc_rcp_rbp`, `rbc_rcp_rp`, `rbp_rp` — see `run_all.sh`). Each contains MATLAB plotting scripts and fit-related outputs.

## Scripts in this folder

| File | Purpose |
|------|---------|
| `run_all.sh` | Runs the full batch for **every** pressure sub-hypothesis |
| `hypothesis_replace_params1.py` | Injects **baseline / default** parameters into **model1** in all plotting `.m` files |
| `hypothesis_replace_params2.py` | Injects **best-fit** parameters into **model2** for each sub-hypothesis |
| `run_plot_scripts.sh` | Runs all plotting MATLAB scripts under each subfolder |
| `copy_model_files_to_sim_files.py` | Copies CSV outputs into `sim_files/` for downstream AICc / preprocessing |
| `param_replace.py` | Shared helpers for editing MATLAB parameter lines |

## Recommended order (run from **this** directory)

Current directory should be **`pressure_hypothesis/`**.

1. `./run_all.sh`  
2. `python hypothesis_replace_params1.py`  
3. `python hypothesis_replace_params2.py`  
4. `./run_plot_scripts.sh`  
5. `python copy_model_files_to_sim_files.py`  

After step 5, compartment-level model CSVs for **global** and **individual** fits are available under `sim_files/` as configured by the copy script.

## AICc

Run **`AICC.ipynb`** from **`../`** (the `final_hypothesis_testing_with_updated_bounds` root), after `sim_files/` here is current.

## Notes

- **model1** = default baseline; **model2** = best-fit hypothesis model in the plotting code.
- The copy script is tuned to this family’s folder names; keep sub-hypothesis directory names aligned with `FOLDERS` / `run_all.sh` lists when adding new cases.
