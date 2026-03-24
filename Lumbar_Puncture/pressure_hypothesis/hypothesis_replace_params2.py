"""
replace_model2_params.py

For each model2 function:
  1. Top-level param declarations  → replaced with DEFAULTS for that fit type
  2. if-block params               → replaced with IF_VALUES for that folder/fit
  3. else-block params             → replaced with DEFAULTS for that fit type
  4. model1 and everything else    → untouched
"""

import re
import shutil
from pathlib import Path

# ── Default parameter sets (used for top-level + else-block) ──────────────────

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

# ── IF-block values (only the params that change during sleep deprivation) ─────

IF_VALUES = {
    "rbc": {
        "global":   {"r_bc": "0.056"},
        "blattner": {"r_bc": "0.026"},
        "lucey":    {"r_bc": "0.085"},
        "liu":      {"r_bc": "0.018"},
    },
    "rcp": {
        "global":   {"r_cp": "0.00526"},
        "blattner": {"r_cp": "0.01028"},
        "lucey":    {"r_cp": "0.01136"},
        "liu":      {"r_cp": "0.00320"},
    },
    "rbc_rcp": {
        "global":   {"r_bc": "0.056", "r_cp": "0.00537"},
        "blattner": {"r_bc": "0.019", "r_cp": "0.01024"},
        "lucey":    {"r_bc": "0.062", "r_cp": "0.01138"},
        "liu":      {"r_bc": "0.018", "r_cp": "0.00320"},
    },
    "rbc_rcp_rbp": {
        "global":   {"r_bc": "0.063", "r_cp": "0.00537", "r_bp": "0.015"},
        "blattner": {"r_bc": "0.020", "r_cp": "0.000001", "r_bp": "0.229"},
        "lucey":    {"r_bc": "0.062", "r_cp": "0.00313", "r_bp": "0.176"},
        "liu":      {"r_bc": "0.021", "r_cp": "0.00320", "r_bp": "0.0149"},
    },
    "rbc_rcp_rp": {
        "global":   {"r_bc": "0.059", "r_cp": "0.00425", "r_p": "0.361"},
        "blattner": {"r_bc": "0.019", "r_cp": "0.0102", "r_p": "0.298"},
        "lucey":    {"r_bc": "0.062", "r_cp": "0.0113", "r_p": "0.327"},
        "liu":      {"r_bc": "0.023", "r_cp": "0.00320", "r_p": "0.452"},
    },
    "rbp_rp": {
        "global":   {"r_bp": "0.010", "r_p": "0.392"},
        "blattner": {"r_bp": "0.010", "r_p": "0.295"},
        "lucey":    {"r_bp": "0.013", "r_p": "0.300"},
        "liu":      {"r_bp": "0.0115", "r_p": "0.427"},
    },
}

FILE_TO_FIT = {
    "global_plot":   "global",
    "blattner_plot": "blattner",
    "lucey_plot":    "lucey",
    "liu_plot":      "liu",
}

FOLDERS = list(IF_VALUES.keys())


# ── Low-level helpers ──────────────────────────────────────────────────────────

def replace_param(text: str, param: str, value: str) -> str:
    """Replace `param = <old>;` anywhere in text. Word-boundary safe."""
    pattern = r'(?m)(^\s*' + re.escape(param) + r'\s*=\s*)([^;]+)(;)'
    new_text, n = re.subn(pattern, r'\g<1>' + value + r'\3', text)
    if n == 0:
        print(f"      [WARN] '{param}' not found — check the file.")
    return new_text


def extract_model2_bounds(content: str):
    m = re.search(r'\bfunction\b[^\n]*\bmodel2\b[^\n]*\n', content)
    if m is None:
        return None
    start = m.start()
    nxt = re.search(r'\bfunction\b', content[m.end():])
    end = m.end() + nxt.start() if nxt else len(content)
    return start, end


def extract_if_else_bounds(block: str):
    """
    Returns (toplevel_end, if_body_start, if_body_end, else_body_start, else_body_end)
    all relative to `block`.

    Layout inside model2:
        <function header line>
        <top-level declarations>      <- block[header_end : toplevel_end]
        if(...)                       <- the if-condition line
            <if-body params>          <- block[if_body_start : if_body_end]
        else
            <else-body params>        <- block[else_body_start : else_body_end]
        end
        <rest of function>
    """
    if_line = re.search(r'\bif\s*\([^)]*\)\s*\n', block)
    if if_line is None:
        return None

    toplevel_end  = if_line.start()
    if_body_start = if_line.end()

    tail = block[if_body_start:]

    else_match = re.search(r'\belse\b', tail)
    if else_match is None:
        return None
    if_body_end     = if_body_start + else_match.start()
    else_body_start = if_body_start + else_match.end() + 1  # skip newline after else

    end_match = re.search(r'(?m)^\s*end\b', block[else_body_start:])
    if end_match is None:
        return None
    else_body_end = else_body_start + end_match.start()

    return toplevel_end, if_body_start, if_body_end, else_body_start, else_body_end


# ── File processor ─────────────────────────────────────────────────────────────

def process_file(filepath: Path, defaults: dict, if_params: dict) -> bool:
    content = filepath.read_text(encoding="utf-8")

    m2_bounds = extract_model2_bounds(content)
    if m2_bounds is None:
        print(f"  [SKIP] model2 not found in {filepath.name}")
        return False

    m2_start, m2_end = m2_bounds
    block = content[m2_start:m2_end]

    bounds = extract_if_else_bounds(block)
    if bounds is None:
        print(f"  [SKIP] if/else structure not found in {filepath.name}")
        return False

    toplevel_end, if_body_start, if_body_end, else_body_start, else_body_end = bounds

    toplevel  = block[:toplevel_end]
    if_body   = block[if_body_start:if_body_end]
    else_body = block[else_body_start:else_body_end]
    remainder = block[else_body_end:]

    # 1. Defaults → top-level declarations
    new_toplevel = toplevel
    for param, val in defaults.items():
        new_toplevel = replace_param(new_toplevel, param, val)

    # 2. IF_VALUES → if-body
    new_if = if_body
    for param, val in if_params.items():
        new_if = replace_param(new_if, param, val)

    # 3. Defaults → else-body
    new_else = else_body
    for param, val in defaults.items():
        new_else = replace_param(new_else, param, val)

    new_block = (
        new_toplevel
        + block[toplevel_end:if_body_start]   # the if(...) condition line
        + new_if
        + block[if_body_end:else_body_start]  # the `else` keyword line
        + new_else
        + remainder
    )

    if new_block == block:
        print(f"  [NO CHANGE] {filepath.name}")
        return False

    new_content = content[:m2_start] + new_block + content[m2_end:]

    shutil.copy2(filepath, filepath.with_suffix(".m.bak"))
    filepath.write_text(new_content, encoding="utf-8")
    print(f"  [UPDATED]  {filepath.name}")
    return True


# ── Validation ─────────────────────────────────────────────────────────────────

def validate_no_placeholders():
    bad = []
    for folder, fits in IF_VALUES.items():
        for fit, params in fits.items():
            for param, val in params.items():
                if val == "???":
                    bad.append(f"{folder}/{fit}/{param}")
    if bad:
        print("\n[ABORT] Placeholders still present:")
        for b in bad:
            print(f"  {b}")
        raise SystemExit(1)


# ── Entry point ────────────────────────────────────────────────────────────────

def main():
    validate_no_placeholders()

    script_dir = Path(__file__).parent
    total_updated = total_skipped = 0

    for folder_name in FOLDERS:
        folder = script_dir / folder_name
        if not folder.is_dir():
            print(f"\n[MISSING FOLDER] {folder} — skipping.")
            continue

        print(f"\n── {folder_name} ──")

        for stem, fit_type in FILE_TO_FIT.items():
            filepath = folder / f"{stem}.m"
            if not filepath.exists():
                print(f"  [MISSING FILE] {filepath.name} — skipping.")
                total_skipped += 1
                continue

            changed = process_file(
                filepath,
                DEFAULTS[fit_type],
                IF_VALUES[folder_name][fit_type]
            )
            total_updated += int(changed)
            total_skipped += int(not changed)

    print(f"\n{'='*55}")
    print(f"Done.  Files updated: {total_updated} | Unchanged/skipped: {total_skipped}")


if __name__ == "__main__":
    main()
