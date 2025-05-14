#!/bin/bash

# List of hypothesis folders
folders=("hypothesis1" "hypothesis2" "hypothesis3" "hypothesis4" "hypothesis5")

# Loop through folders and run the existing bash script in each
for folder in "${folders[@]}"; do
    echo "============================================"
    echo "Running MATLAB jobs in $folder..."
    (
        cd "$folder" || exit
        chmod +x *.sh  # Ensure the script is executable
        ./$(basename *.sh)  # Run the only .sh script in the folder
    )
    echo "Finished $folder"
    echo "============================================"
done

echo "All folders processed."

