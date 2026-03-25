#!/bin/bash

# Define the MATLAB scripts to run in sequence
SCRIPT1="global_fit.m"
SCRIPT2="blattner_fit.m"
SCRIPT3="lucey_fit.m"
SCRIPT4="liu_fit.m"

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
run_matlab_script "$SCRIPT2"
run_matlab_script "$SCRIPT3"
run_matlab_script "$SCRIPT4"
echo "All MATLAB scripts completed."
