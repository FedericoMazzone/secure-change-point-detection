#!/bin/bash

# Root directories
input_root="data/generated"
output_root="logs"

# Find all .csv files
find "$input_root" -type f -name "*.csv" | while read -r csv_file; do
    # Compute the relative path from input_root
    relative_path="${csv_file#$input_root/}"
    
    # Replace .csv with .out and prepend logs/ path
    output_file="$output_root/${relative_path%.csv}.out"
    
    # Create the output directory if it doesn't exist
    output_dir=$(dirname "$output_file")
    mkdir -p "$output_dir"

    # Run the command and save stdout and stderr
    ./build/cpd -f "$csv_file" -l -v > "$output_file" 2>&1
done