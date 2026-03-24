"""
replace_model1_params.py

Replaces the default parameter values inside the model1 function in every
.m file across the ablation study folder structure.

Folder structure expected (all relative to this script's location):
    rbc/
    rbc_rcp/
    rbc_rcp_rp/
    rbc_rcp_rbp/
    rbp_rp/
    rcp/

Each folder contains: global_plot.m, blattner_plot.m, lucey_plot.m, liu_plot.m

Rules:
  - global_plot.m   → global default values
  - blattner_plot.m → blattner default values
  - lucey_plot.m    → lucey default values
  - liu_plot.m      → liu default values

Only the model1 function block is modified; model2 is left untouched.
"""

import re
import os
import shutil
from pathlib import Path

# ── Default parameter sets ─────────────────────────────────────────────────────

DEFAULTS = {
    "global": {
        "r_bc":     "0.038",
        "r_bp":     "0.014",
        "r_cp":     "0.00537",
        "sigma_bc": "1.131",
        "sigma_bp": "1.768",
        "sigma_cp": "6.100",
        "sigma_p":  "4.253",
        "A":        "16.203",
        "sigma_A":  "0.772",
        "r_p":      "0.427",
    },
    "blattner": {
        "r_bc":     "0.019",
        "r_bp":     "0.034",
        "r_cp":     "0.0154",
        "sigma_bc": "1.660",
        "sigma_bp": "1.816",
        "sigma_cp": "5.740",
        "sigma_p":  "3.610",
        "A":        "84.523",
        "sigma_A":  "0.633",
        "r_p":      "0.298",
    },
    "lucey": {
        "r_bc":     "0.062",
        "r_bp":     "0.040",
        "r_cp":     "0.0156",
        "sigma_bc": "1.002",
        "sigma_bp": "2.479",
        "sigma_cp": "4.087",
        "sigma_p":  "4",
        "A":        "55.780",
        "sigma_A":  "0.485",
        "r_p":      "0.300",
    },
    "liu": {
        "r_bc":     "0.015",
        "r_bp":     "0.014",
        "r_cp":     "0.00320",
        "sigma_bc": "1.110",
        "sigma_bp": "1.297",
        "sigma_cp": "6.552",
        "sigma_p":  "3.055",
        "A":        "14.450",
        "sigma_A":  "0.750",
        "r_p":      "0.475",
    },
}

# Map filename stem → default key
FILE_TO_DEFAULT = {
    "global_plot":   "global",
    "blattner_plot": "blattner",
    "lucey_plot":    "lucey",
    "liu_plot":      "liu",
}

FOLDERS = ["rbc_rcp_all_sigma","rbc_rcp_rbp_sbc_scp_sbp","rbc_rcp_sbc_scp","rbc_sbc_scp","rbc_sbc_scp_sp","rbc_scp_sp","rcp_sbc_scp"]

# Parameters in the order they appear in model1 (used for replacement)
PARAM_NAMES = ["r_bc", "r_bp", "r_cp", "sigma_bc", "sigma_bp", "sigma_cp",
               "sigma_p", "A", "sigma_A", "r_p"]


# ── Core helpers ───────────────────────────────────────────────────────────────

def extract_model1_block(content: str):
    """
    Return (start_idx, end_idx) of the model1 function body in the file.
    Looks for 'function ... = model1(' and ends at the next 'function ' header
    or end of file. Returns None if model1 is not found.
    """
    # Match the start of model1 (handles any return signature)
    start_match = re.search(r'\bfunction\b[^\n]*\bmodel1\b[^\n]*\n', content)
    if start_match is None:
        return None

    start = start_match.start()

    # Find the next 'function ' keyword after model1 starts (= end of model1)
    next_func = re.search(r'\bfunction\b', content[start_match.end():])
    if next_func:
        end = start_match.end() + next_func.start()
    else:
        end = len(content)

    return start, end


def replace_param_in_block(block: str, param: str, new_value: str) -> str:
    """
    Replace   param = <old_value>;
    with      param = <new_value>;
    inside the given text block, using a word-boundary-aware regex.
    Handles optional leading whitespace.
    """
    pattern = r'(?m)(^\s*' + re.escape(param) + r'\s*=\s*)([^;]+)(;)'
    replacement = r'\g<1>' + new_value + r'\3'
    new_block, n = re.subn(pattern, replacement, block)
    if n == 0:
        print(f"    [WARN] Parameter '{param}' not found in block.")
    return new_block


def process_file(filepath: Path, defaults: dict) -> bool:
    """
    Read the .m file, locate model1, replace all parameter values, write back.
    Creates a .bak backup before modifying. Returns True if file was changed.
    """
    content = filepath.read_text(encoding="utf-8")
    bounds = extract_model1_block(content)

    if bounds is None:
        print(f"  [SKIP] model1 not found in {filepath}")
        return False

    start, end = bounds
    model1_block = content[start:end]
    modified_block = model1_block

    for param in PARAM_NAMES:
        if param in defaults:
            modified_block = replace_param_in_block(
                modified_block, param, defaults[param]
            )

    if modified_block == model1_block:
        print(f"  [NO CHANGE] {filepath}")
        return False

    # Backup original
    shutil.copy2(filepath, filepath.with_suffix(".m.bak"))

    new_content = content[:start] + modified_block + content[end:]
    filepath.write_text(new_content, encoding="utf-8")
    print(f"  [UPDATED]  {filepath}")
    return True


# ── Entry point ────────────────────────────────────────────────────────────────

def main():
    script_dir = Path(__file__).parent
    total_updated = 0
    total_skipped = 0

    for folder_name in FOLDERS:
        folder = script_dir / folder_name
        if not folder.is_dir():
            print(f"[MISSING FOLDER] {folder} — skipping.")
            continue

        print(f"\n── {folder_name} ──")

        for stem, default_key in FILE_TO_DEFAULT.items():
            filepath = folder / f"{stem}.m"
            if not filepath.exists():
                print(f"  [MISSING FILE] {filepath} — skipping.")
                total_skipped += 1
                continue

            defaults = DEFAULTS[default_key]
            changed = process_file(filepath, defaults)
            total_updated += int(changed)
            total_skipped += int(not changed)

    print(f"\n{'='*50}")
    print(f"Done. Files updated: {total_updated} | Unchanged/skipped: {total_skipped}")


if __name__ == "__main__":
    main()
