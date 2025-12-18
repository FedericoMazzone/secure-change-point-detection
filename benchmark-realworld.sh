#!/bin/bash

executable="./build/cpd"                    # Path to the CPD executable
output_file="logs/benchmark-realworld.out"  # Output CSV file
datasets=(                                  # Real-world datasets to test
    "data/meditation.csv"
    "data/patient5-S2toS3.csv"
    "data/traffic_network.csv"
)
type="frequency"                            # Type of change points to test



# Create logs directory if it doesn't exist
if [ ! -d "logs" ]; then
    mkdir logs
fi

# Create output file with header if it doesn't exist
if [ ! -s "$output_file" ]; then
    echo "dataset,type,length,encrypt_runtime(s),turning_runtime(s),partials_runtime(s),argmax_runtime(s),cpd_runtime(s),memory(MB),base_acc,comp_acc" > "$output_file"
fi

# Run experiments for different datasets
for filepath in "${datasets[@]}"; do

    echo "Experiment: dataset=$filepath, type=$type"

    dataset_name=$(basename "$filepath" .csv)
    log_file="logs/${dataset_name}-${type}.log"
    echo "Log file: $log_file"

    command="$executable -f $filepath -l -v -t $type"
    echo "Running command: $command"

    # Run the command and measure memory usage
    $command > "$log_file" 2>&1 &
    pid=$!
    max_res=0
    while kill -0 "$pid" 2>/dev/null; do
        current_res=$(ps -o rss= -p "$pid" 2>/dev/null)
        current_res=${current_res:-0}
        if [ "$current_res" -gt "$max_res" ]; then
            max_res="$current_res"
        fi
        sleep 0.5
    done
    final_max_res=$(awk "BEGIN {printf \"%.2f\", $max_res/1024}")

    # Extract runtimes, length, and accuracy from log file
    length=$(grep -m1 -oP 'Number of time points:\s*\K[0-9]+' "$log_file")
    encrypt_runtime=$(grep -m1 -oP 'Encrypting time-series data\.\.\.DONE \(\K[0-9]+\.[0-9]+(?=s\))' "$log_file")
    turning_runtime=$(grep -m1 -oP 'Turning rates computed in \K[0-9]+\.[0-9]+(?=s)' "$log_file")
    partials_runtime=$(grep -m1 -oP 'Partial sums computed in \K[0-9]+\.[0-9]+(?=s)' "$log_file")
    argmax_runtime=$(grep -m1 -oP 'Argmax computed in \K[0-9]+\.[0-9]+(?=s)' "$log_file")
    cpd_runtime=$(grep -m1 -oP 'CPD runtime:\s*\K[0-9]+\.[0-9]+(?=s)' "$log_file")
    true_index=$(grep -m1 -oP 'Label:\s*\K[0-9]+' "$log_file")
    base_index=$(grep -m1 -oP 'Expected change-point:\s*\K[0-9]+' "$log_file")
    comp_index=$(grep -m1 -oP 'Computed change-point:\s*\K[0-9]+' "$log_file")
    base_acc=$(awk -v b="$base_index" -v t="$true_index" 'BEGIN{printf "%.4f", (b>t?b-t:t-b)/t}')
    comp_acc=$(awk -v c="$comp_index" -v t="$true_index" 'BEGIN{printf "%.4f", (c>t?c-t:t-c)/t}')
    
    # Append results to output file
    echo "${dataset_name},${type},${length},${encrypt_runtime},${turning_runtime},${partials_runtime},${argmax_runtime},${cpd_runtime},${final_max_res},${base_acc},${comp_acc}" | tee -a "$output_file"

    echo "==============================="

done

echo "Benchmarking completed. Results saved in $output_file."
