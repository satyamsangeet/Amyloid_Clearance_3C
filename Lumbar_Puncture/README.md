# Final hypothesis testing

This directory bundles three parallel hypothesis. Each hypothesis is in its own folder and contains multiple **sub-hypothesis** folders (one ablation / variant each). The workflow is the same in every hypothesis: fit everything, refresh plotting MATLAB with baseline and best-fit parameters, regenerate simulation CSVs, copy them into `sim_files/`, then compute AICc at the repository root.

## Layout

| Folder | Role |
|--------|------|
| `combined_hypothesis/` | Combined hypotheses |
| `sleep_hypothesis/` | Sleep-related hypotheses |
| `pressure_hypothesis/` | Pressure hypotheses |
| `AICC.ipynb` | AICc calculation (run after CSVs are in each `sim_files/`) |

Each hypothesis folder includes:

- `run_all.sh` — runs hypothesis testing across **all** sub-hypothesis folders
- `run_plot_scripts.sh` — runs every MATLAB plotting script in those folders
- `hypothesis_replace_params1.py` — writes **baseline / default** parameters into **model1** in the `.m` plotting scripts
- `hypothesis_replace_params2.py` — writes **best-fit** parameters into **model2** for each sub-hypothesis
- `copy_model_files_to_sim_files.py` — collects model-output CSVs into `sim_files/` (see that folder’s README for details)
- `param_replace.py` — helper utilities for editing parameter lines in MATLAB files (supporting code; not a required manual step in the list below)

## End-to-end workflow

Repeat the steps below **inside each** of `combined_hypothesis`, `sleep_hypothesis`, and `pressure_hypothesis`, always working from **that folder** as the current directory.

You need to make the bash script (.sh) executable by doing
```
chmod +x run_all.sh
chmod +x run_plot_scripts.sh
```

1. **`./run_all.sh`**  
   Runs every sub-hypothesis’s fitting / batch pipeline. Wait until it finishes.

2. **`python hypothesis_replace_params1.py`**  
   Updates **model1** in the plotting `.m` files to the correct **baseline / default** parameter sets (global vs dataset-specific scripts, depending on the file).

3. **`python hypothesis_replace_params2.py`**  
   Updates **model2** to the **best-fit** parameters for **each** sub-hypothesis.

4. **`./run_plot_scripts.sh`**  
   Requires MATLAB on your `PATH`. Executes each `*_plot.m` script in every sub-hypothesis folder. This produces CSVs with the model response for **each compartment** and fit type (**global** and **individual** fits: blattner, lucey, liu, etc., as named by the scripts).

5. **`python copy_model_files_to_sim_files.py`**  
   Copies the generated CSVs from the sub-hypothesis directories into the local **`sim_files/`** tree (with the subfolders that notebook preprocessing expects).

6. **AICc**  
   Open and run **`AICC.ipynb`** in this **root** directory. It expects the aggregated CSVs under each hypothesis folder’s `sim_files/` (and related preprocess paths) to be up to date.

## Tips

- Run steps 2–5 **after** `run_all.sh` completes so fitted values and outputs match the latest runs.
- If a step fails for one hypothesis, you can still process the others; AICc may need all intended `sim_files/` trees populated consistently.
- For command syntax, prerequisites (Python, MATLAB), and family-specific folder lists, see **`README.md` inside each hypothesis folder**.
