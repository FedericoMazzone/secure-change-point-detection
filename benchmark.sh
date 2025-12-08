#!/bin/bash

executable="./build/cpd"

output_file="benchmark_cpd.out"

if [ ! -s "$output_file" ]; then
    echo "dataset,type,length,encrypt_runtime(s),turning_runtime(s),partials_runtime(s),argmax_runtime(s),cpd_runtime(s),memory(MB)" > "$output_file"
fi

if [ ! -d "logs" ]; then
    mkdir logs
fi

files=(
    "data/meditation.csv"
    "data/patient5-S2toS3.csv"
    "data/traffic_network.csv"

    "data/generated/mean-change/time_series_n1000_normal_0_1_to_normal_1_1.csv"
    "data/generated/mean-change/time_series_n10000_normal_0_1_to_normal_1_1.csv"
    "data/generated/mean-change/time_series_n20000_normal_0_1_to_normal_1_1.csv"
    "data/generated/mean-change/time_series_n40000_normal_0_1_to_normal_1_1.csv"
    "data/generated/mean-change/time_series_n100000_normal_0_1_to_normal_1_1.csv"
    "data/generated/mean-change/time_series_n500000_normal_0_1_to_normal_1_1.csv"
    "data/generated/mean-change/time_series_n1000000_normal_0_1_to_normal_1_1.csv"

    "data/generated/variance-change/time_series_n1000_normal_0_1_to_normal_0_2.csv"
    "data/generated/variance-change/time_series_n10000_normal_0_1_to_normal_0_2.csv"
    "data/generated/variance-change/time_series_n20000_normal_0_1_to_normal_0_2.csv"
    "data/generated/variance-change/time_series_n40000_normal_0_1_to_normal_0_2.csv"
    "data/generated/variance-change/time_series_n100000_normal_0_1_to_normal_0_2.csv"
    "data/generated/variance-change/time_series_n500000_normal_0_1_to_normal_0_2.csv"
    "data/generated/variance-change/time_series_n1000000_normal_0_1_to_normal_0_2.csv"

    "data/generated/frequency-change/time_series_gaussian_n1000.csv"
    "data/generated/frequency-change/time_series_gaussian_n10000.csv"
    "data/generated/frequency-change/time_series_gaussian_n20000.csv"
    "data/generated/frequency-change/time_series_gaussian_n40000.csv"
    "data/generated/frequency-change/time_series_gaussian_n100000.csv"
    "data/generated/frequency-change/time_series_gaussian_n500000.csv"
    "data/generated/frequency-change/time_series_gaussian_n1000000.csv"
)

types=(
    "frequency"
    "frequency"
    "frequency"

    "mean"
    "mean"
    "mean"
    "mean"
    "mean"
    "mean"
    "mean"

    "variance"
    "variance"
    "variance"
    "variance"
    "variance"
    "variance"
    "variance"

    "frequency"
    "frequency"
    "frequency"
    "frequency"
    "frequency"
    "frequency"
    "frequency"
)

if [ "${#files[@]}" -ne "${#types[@]}" ]; then
    echo "Error: files[] and types[] have different lengths."
    exit 1
fi

for i in "${!files[@]}"; do
    file="${files[$i]}"
    type="${types[$i]}"

    if [ ! -f "$file" ]; then
        echo "Warning: file $file not found, skipping."
        continue
    fi

    dataset_name=$(basename "$file")
    echo "Running CPD on $dataset_name (type=$type)..."

    log_file="logs/${dataset_name}-${type}.log"

    command="$executable -f $file -l -v -t $type"
    echo "$command"

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

    length=$(grep -oP 'Number of time points:\s*\K[0-9]+' "$log_file")
    encrypt_runtime=$(grep -oP 'Encrypting time-series data\.\.\.DONE \(\K[0-9]+\.[0-9]+(?=s\))' "$log_file")
    turning_runtime=$(grep -oP 'Turning rates computed in \K[0-9]+\.[0-9]+(?=s)' "$log_file")
    partials_runtime=$(grep -oP 'Partial sums computed in \K[0-9]+\.[0-9]+(?=s)' "$log_file")
    argmax_runtime=$(grep -oP 'Argmax computed in \K[0-9]+\.[0-9]+(?=s)' "$log_file")
    cpd_runtime=$(grep -oP 'CPD runtime:\s*\K[0-9]+\.[0-9]+(?=s)' "$log_file")

    echo "${dataset_name},${type},${length},${encrypt_runtime},${turning_runtime},${partials_runtime},${argmax_runtime},${cpd_runtime},${final_max_res}" | tee -a "$output_file"

done

echo "Benchmarking completed. Results saved in $output_file."



# ./build/cpd -f data/meditation.csv -l -v -t frequency
# ./build/cpd -f data/patient5-S2toS3.csv -l -v -t frequency
# ./build/cpd -f data/traffic_network.csv -l -v -t frequency

# ./build/cpd -f data/generated/mean-change/time_series_n1000_normal_0_1_to_normal_1_1.csv -l -v -t mean
# ./build/cpd -f data/generated/mean-change/time_series_n10000_normal_0_1_to_normal_1_1.csv -l -v -t mean
# ./build/cpd -f data/generated/mean-change/time_series_n20000_normal_0_1_to_normal_1_1.csv -l -v -t mean
# ./build/cpd -f data/generated/mean-change/time_series_n40000_normal_0_1_to_normal_1_1.csv -l -v -t mean
# ./build/cpd -f data/generated/mean-change/time_series_n100000_normal_0_1_to_normal_1_1.csv -l -v -t mean
# ./build/cpd -f data/generated/mean-change/time_series_n500000_normal_0_1_to_normal_1_1.csv -l -v -t mean
# ./build/cpd -f data/generated/mean-change/time_series_n1000000_normal_0_1_to_normal_1_1.csv -l -v -t mean

# ./build/cpd -f data/generated/variance-change/time_series_n1000_normal_0_1_to_normal_0_2.csv -l -v -t variance
# ./build/cpd -f data/generated/variance-change/time_series_n10000_normal_0_1_to_normal_0_2.csv -l -v -t variance
# ./build/cpd -f data/generated/variance-change/time_series_n20000_normal_0_1_to_normal_0_2.csv -l -v -t variance
# ./build/cpd -f data/generated/variance-change/time_series_n40000_normal_0_1_to_normal_0_2.csv -l -v -t variance
# ./build/cpd -f data/generated/variance-change/time_series_n100000_normal_0_1_to_normal_0_2.csv -l -v -t variance
# ./build/cpd -f data/generated/variance-change/time_series_n500000_normal_0_1_to_normal_0_2.csv -l -v -t variance
# ./build/cpd -f data/generated/variance-change/time_series_n1000000_normal_0_1_to_normal_0_2.csv -l -v -t variance

# ./build/cpd -f data/generated/frequency-change/time_series_gaussian_n1000.csv -l -v -t frequency
# ./build/cpd -f data/generated/frequency-change/time_series_gaussian_n10000.csv -l -v -t frequency
# ./build/cpd -f data/generated/frequency-change/time_series_gaussian_n20000.csv -l -v -t frequency
# ./build/cpd -f data/generated/frequency-change/time_series_gaussian_n40000.csv -l -v -t frequency
# ./build/cpd -f data/generated/frequency-change/time_series_gaussian_n100000.csv -l -v -t frequency
# ./build/cpd -f data/generated/frequency-change/time_series_gaussian_n500000.csv -l -v -t frequency
# ./build/cpd -f data/generated/frequency-change/time_series_gaussian_n1000000.csv -l -v -t frequency
