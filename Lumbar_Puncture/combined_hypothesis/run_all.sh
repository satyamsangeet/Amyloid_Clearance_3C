#!/bin/bash

# List of hypothesis folders
folders=("rbc_rcp_all_sigma" "rbc_rcp_rbp_sbc_scp_sbp" "rbc_rcp_sbc_scp" "rbc_sbc_scp" "rbc_sbc_scp_sp" "rbc_scp_sp" "rcp_sbc_scp")

# Total number of folders
total_folders=${#folders[@]}
current_folder=0

# Spinner animation frames
spinner_frames=('⠋' '⠙' '⠹' '⠸' '⠼' '⠴' '⠦' '⠧' '⠇' '⠏')
spinner_index=0

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
PURPLE='\033[0;35m'
BOLD='\033[1m'
NC='\033[0m' # No Color

# Function to draw modern progress bar
draw_progress() {
    local current=$1
    local total=$2
    local folder_name=$3
    local width=40
    local percentage=$((current * 100 / total))
    local filled=$((current * width / total))
    local empty=$((width - filled))
    
    printf "\r\033[K"
    
    local spinner="${spinner_frames[$spinner_index]}"
    spinner_index=$(( (spinner_index + 1) % ${#spinner_frames[@]} ))
    
    printf "${CYAN}${spinner}${NC} ${BOLD}Processing:${NC} ${YELLOW}%s${NC}\n" "$folder_name"
    printf "${BLUE}╭─────────────────────────────────────────────────────────────╮${NC}\n"
    printf "${BLUE}│${NC} Progress: ["
    
    if [ $percentage -lt 33 ]; then
        printf "${RED}%*s${NC}" $filled | tr ' ' '█'
    elif [ $percentage -lt 66 ]; then
        printf "${YELLOW}%*s${NC}" $filled | tr ' ' '█'
    else
        printf "${GREEN}%*s${NC}" $filled | tr ' ' '█'
    fi
    
    printf "%*s" $empty | tr ' ' '░'
    printf "] ${BOLD}%3d%%${NC} ${BLUE}│${NC}\n" $percentage
    printf "${BLUE}│${NC} Completed: ${GREEN}%2d${NC}/${BOLD}%2d${NC} folders                                ${BLUE}│${NC}\n" $current $total
    printf "${BLUE}╰─────────────────────────────────────────────────────────────╯${NC}\n"
    
    printf "\033[4A"
}

# Function to show completion status
show_completion() {
    local folder_name=$1
    local status=$2
    
    printf "\r\033[K"
    if [ "$status" = "success" ]; then
        printf "${GREEN}✓${NC} ${BOLD}Completed:${NC} ${folder_name}\n"
    else
        printf "${RED}✗${NC} ${BOLD}Failed:${NC} ${folder_name}\n"
    fi
}

# Function to run MATLAB job in folder
run_matlab_job() {
    local folder=$1
    local start_time=$(date +%s)
    
    cd "$folder" || return 1
    chmod +x *.sh 2>/dev/null
    
    # Find and run the .sh script
    local script_file=$(ls *.sh 2>/dev/null | head -n1)
    if [[ -n "$script_file" ]]; then
        ./"$script_file" > /dev/null 2>&1 &
        local job_pid=$!

        while kill -0 $job_pid 2>/dev/null; do
            draw_progress $current_folder $total_folders "$folder"
            sleep 0.1
        done
        
        wait $job_pid
        local exit_status=$?
        
        cd - > /dev/null
        return $exit_status
    else
        cd - > /dev/null
        return 1
    fi
}

clear
printf "${BOLD}${PURPLE}╔═══════════════════════════════════════════════════════════════╗${NC}\n"
printf "${BOLD}${PURPLE}║${NC}                    ${BOLD}MATLAB Jobs Processor${NC}                    ${BOLD}${PURPLE}║${NC}\n"
printf "${BOLD}${PURPLE}╚═══════════════════════════════════════════════════════════════╝${NC}\n\n"

printf "${BOLD}Starting processing of ${CYAN}%d${NC} folders...${NC}\n\n" $total_folders

printf "\n\n\n\n"

# Looping through folders
for folder in "${folders[@]}"; do
    ((current_folder++))
    
    # Running job
    if run_matlab_job "$folder"; then
        printf "\r\033[K\033[4B"
        show_completion "$folder" "success"
    else
        printf "\r\033[K\033[4B"
        show_completion "$folder" "failed"
    fi
    
    sleep 0.5
done

# Final completion message
printf "\n${BOLD}${GREEN}╔═══════════════════════════════════════════════════════════════╗${NC}\n"
printf "${BOLD}${GREEN}║${NC}                        ${BOLD}PROCESSING COMPLETE!${NC}                   ${BOLD}${GREEN}║${NC}\n"
printf "${BOLD}${GREEN}║${NC}              Successfully processed ${CYAN}%2d${NC}/${BOLD}%2d${NC} folders              ${BOLD}${GREEN}║${NC}\n" $total_folders $total_folders
printf "${BOLD}${GREEN}╚═══════════════════════════════════════════════════════════════╝${NC}\n"

printf "\n${BOLD}All folders have been processed!${NC} 🎉\n"
