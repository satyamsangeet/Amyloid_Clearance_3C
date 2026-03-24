#!/bin/bash

# Define the MATLAB scripts to run in sequence
SCRIPT1="global_fit1.m"

# Function to run a MATLAB script
run_matlab_script() {
    script=$1
    echo "Running MATLAB script: $script"
    matlab -nodisplay -nosplash -nodesktop -r "try, run('$script'), catch e, disp(getReport(e)), end, exit"
    echo "Completed: $script"
    echo "----------------------------------------"
}

# Run each script in sequence
echo "Starting MATLAB batch processing..."
run_matlab_script "$SCRIPT1"
echo "All MATLAB scripts completed."
