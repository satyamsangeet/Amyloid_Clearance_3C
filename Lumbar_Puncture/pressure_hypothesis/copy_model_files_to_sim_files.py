"""
copy_to_sim_files.py

For each ablation folder (rbc, rcp, rbc_rcp, etc.):
  - Copies {fit}_model2_{folder}.csv  →  sim_files/{fit}/

Run from the parent directory containing both the ablation folders and sim_files/.
"""

import shutil
from pathlib import Path

FOLDERS = ["rbc", "rcp", "rbc_rcp", "rbc_rcp_rp", "rbc_rcp_rbp", "rbp_rp"]
FIT_TYPES = ["global", "blattner", "lucey", "liu"]

def main():
    base = Path(__file__).parent
    sim_files = base / "sim_files"

    total_copied = 0
    total_missing = 0

    for folder_name in FOLDERS:
        src_folder = base / folder_name
        if not src_folder.is_dir():
            print(f"[MISSING FOLDER] {src_folder} — skipping.")
            continue

        print(f"\n── {folder_name} ──")

        for fit in FIT_TYPES:
            filename = f"{fit}_model2_{folder_name}.csv"
            src = src_folder / filename
            dst_dir = sim_files / fit
            dst = dst_dir / filename

            if not src.exists():
                print(f"  [MISSING] {src}")
                total_missing += 1
                continue

            if not dst_dir.exists():
                print(f"  [MISSING DEST DIR] {dst_dir} — skipping {filename}.")
                total_missing += 1
                continue

            shutil.copy2(src, dst)
            print(f"  [COPIED] {filename}  →  sim_files/{fit}/")
            total_copied += 1

    print(f"\n{'='*55}")
    print(f"Done.  Copied: {total_copied} | Missing/skipped: {total_missing}")

if __name__ == "__main__":
    main()
