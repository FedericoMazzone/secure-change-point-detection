#!/bin/bash

executable="./build/cpd"                    # Path to the CPD executable
output_file="logs/benchmark-synthetic.out"  # Output CSV file
n_values=(1000 10000 100000 1000000)        # Lengths of time series to test
types=("mean" "variance" "frequency")       # Types of change points to test



# Create logs directory if it doesn't exist
if [ ! -d "logs" ]; then
    mkdir logs
fi

# Create output file with header if it doesn't exist
if [ ! -s "$output_file" ]; then
    echo "dataset,type,length,encrypt_runtime(s),turning_runtime(s),partials_runtime(s),argmax_runtime(s),cpd_runtime(s),memory(MB)" > "$output_file"
fi

# Run experiments for different n values and types
for type in "${types[@]}"; do

    for n in "${n_values[@]}"; do

        echo "Experiment: type=$type, n=$n"

        # Generate synthetic data
        case $type in
            "mean")
                filepath=$(python3 ./src/data_utils/ts_gen_distr_change.py --dist1 normal --dist2 normal --n "$n" --mu1 0.0 --mu2 1.0 --sigma1 1.0 --sigma2 1.0 | grep -m1 -oP 'Time series saved to \K.*')
                ;;
            "variance")
                filepath=$(python3 ./src/data_utils/ts_gen_distr_change.py --dist1 normal --dist2 normal --n "$n" --mu1 0.0 --mu2 0.0 --sigma1 1.0 --sigma2 2.0 | grep -m1 -oP 'Time series saved to \K.*')
                ;;
            "frequency")
                filepath=$(python3 ./src/data_utils/ts_gen_freq_change.py --distribution normal --n "$n" | grep -m1 -oP 'Time series saved to \K.*')
                ;;
            *)
                echo "Unknown type: $type"
                exit 1
                ;;
        esac
        echo "Generated file: $filepath"

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

        # Extract runtimes and length from log file
        length=$(grep -m1 -oP 'Number of time points:\s*\K[0-9]+' "$log_file")
        encrypt_runtime=$(grep -m1 -oP 'Encrypting time-series data\.\.\.DONE \(\K[0-9]+\.[0-9]+(?=s\))' "$log_file")
        turning_runtime=$(grep -m1 -oP 'Turning rates computed in \K[0-9]+\.[0-9]+(?=s)' "$log_file")
        partials_runtime=$(grep -m1 -oP 'Partial sums computed in \K[0-9]+\.[0-9]+(?=s)' "$log_file")
        argmax_runtime=$(grep -m1 -oP 'Argmax computed in \K[0-9]+\.[0-9]+(?=s)' "$log_file")
        cpd_runtime=$(grep -m1 -oP 'CPD runtime:\s*\K[0-9]+\.[0-9]+(?=s)' "$log_file")
        echo "${dataset_name},${type},${length},${encrypt_runtime},${turning_runtime},${partials_runtime},${argmax_runtime},${cpd_runtime},${final_max_res}" | tee -a "$output_file"

        # Append results to output file
        echo "==============================="

    done

done

echo "Benchmarking completed. Results saved in $output_file."
