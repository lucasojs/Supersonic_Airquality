#!/bin/bash

# Define the directories to search
dirs=("CH4_Fixed" "CH4_Flux")

# Loop through each directory
for dir in "${dirs[@]}"; do
  # Find subdirectories containing 'Std' in their names
  find "$dir" -type d -name "*Std*" | while read -r subdir; do
    # Check if the 'restarts' directory exists within the found subdirectory
    if [ -d "$subdir/restarts" ]; then
      # Loop through 'cspec' and 'tracers' directories within 'restarts'
      for type in "cspec" "tracers"; do
        if [ -d "$subdir/restarts/$type" ]; then
          # Find and remove files except January data
          if [ "$type" == "cspec" ]; then
            find "$subdir/restarts/$type" -type f -name "spec_rst.merra_4x5_UCX.*" ! -name "*010100" -exec rm -f {} \;
          elif [ "$type" == "tracers" ]; then
            find "$subdir/restarts/$type" -type f -name "trac_rst.merra_L72.T93.*" ! -name "*01010000" -exec rm -f {} \;
          fi
        fi
      done
    fi
  done
done

