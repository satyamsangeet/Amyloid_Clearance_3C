# Sleep hypothesis

Sub-hypotheses live as immediate subfolders of this directory (e.g. `all_sigma`, `sbc`, `sbc_scp`, … — the full list is in `run_all.sh`). Each subfolder typically contains MATLAB plotting scripts (`global_plot.m`, `blattner_plot.m`, `lucey_plot.m`, `liu_plot.m`) and outputs from fitting.

## Scripts in this folder

| File | Purpose |
|------|---------|
| `run_all.sh` | Runs hypothesis testing for **all** sleep sub-hypothesis folders enumerated in the script |
| `hypothesis_replace_params1.py` | Sets **model1** to **baseline / default** parameters in the plotting `.m` files |
| `hypothesis_replace_params2.py` | Sets **model2** to **best-fit** parameters for each sub-hypothesis |
| `run_plot_scripts.sh` | Executes every plotting `.m` in each subfolder (MATLAB required) |
| `copy_model_files_to_sim_files.py` | Copies generated model-output CSVs into local `sim_files/` |
| `param_replace.py` | Utilities for programmatic parameter edits in MATLAB files |

## Recommended order (run from **this** directory)

Use **`sleep_hypothesis/`** as the current working directory.

1. `./run_all.sh`  
2. `python hypothesis_replace_params1.py`  
3. `python hypothesis_replace_params2.py`  
4. `./run_plot_scripts.sh`  
5. `python copy_model_files_to_sim_files.py`  

This produces CSVs with model responses per compartment for **global** and **individual** fits, then stages **model2** outputs under `sim_files/` for notebooks and AICc.

## AICc

Open **`AICC.ipynb`** in the **root** project folder (`final_hypothesis_testing_with_updated_bounds/`) after updating `sim_files/` here (and in other families if you use them).

## Notes

- Run steps 2–5 only after `run_all.sh` has finished so replacements match the latest fits.
- `param_replace.py` is support code; the standard workflow uses the two `hypothesis_replace_params*.py` scripts only.
