#!/bin/bash
# run_all_matlab.sh
# Runs all 4 MATLAB scripts in each of the 6 ablation folders sequentially.
# Run from the parent directory containing the ablation folders.

FOLDERS=("rbc_rcp_all_sigma" "rbc_rcp_rbp_sbc_scp_sbp" "rbc_rcp_sbc_scp" "rbc_sbc_scp" "rbc_sbc_scp_sp" "rbc_scp_sp" "rcp_sbc_scp")
FILES=("global_plot" "blattner_plot" "lucey_plot" "liu_plot")

PARENT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

total=0
failed=0
succeeded=0

echo "======================================================"
echo " MATLAB batch runner"
echo " Parent dir: $PARENT_DIR"
echo "======================================================"

for folder in "${FOLDERS[@]}"; do
    folder_path="$PARENT_DIR/$folder"

    if [ ! -d "$folder_path" ]; then
        echo ""
        echo "[MISSING FOLDER] $folder_path — skipping."
        continue
    fi

    echo ""
    echo "── $folder ──"

    for script in "${FILES[@]}"; do
        m_file="$folder_path/${script}.m"

        if [ ! -f "$m_file" ]; then
            echo "  [MISSING FILE] ${script}.m — skipping."
            ((failed++))
            ((total++))
            continue
        fi

        echo "  [RUNNING] $folder/${script}.m"
        total=$((total + 1))

        # Run MATLAB in no-display, no-splash, batch mode.
        # cd into the folder first so any relative file paths in the script resolve correctly.
        matlab -nodisplay -nosplash -nodesktop -batch \
            "cd('$folder_path'); run('${script}.m'); exit;" \
            2>&1 | sed 's/^/    /'   # indent MATLAB output for readability

        exit_code=${PIPESTATUS[0]}

        if [ $exit_code -eq 0 ]; then
            echo "  [DONE] $folder/${script}.m"
            ((succeeded++))
        else
            echo "  [FAILED] $folder/${script}.m  (exit code: $exit_code)"
            ((failed++))
        fi
    done
done

echo ""
echo "======================================================"
echo " Finished.  Succeeded: $succeeded / $total  |  Failed: $failed"
echo "======================================================"
