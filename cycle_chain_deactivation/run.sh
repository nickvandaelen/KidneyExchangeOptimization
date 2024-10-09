#!/bin/bash

script_dir="$(cd "$(dirname "$0")" && pwd)"
instance_dir="$script_dir/../Instance Files"
solutions_dir="$script_dir/Solutions"

mkdir -p "$solutions_dir"

file=""
all_flag=false
skip_patterns=("A.*")
generalized_flag=false
max_cycle_length=""
max_chain_length=""

while getopts ":f:as:c:h:g" opt; do
  case ${opt} in
    f )
      file=$OPTARG
      ;;
    a )
      all_flag=true
      ;;
    s )
      skip_patterns+=("$OPTARG")
      ;;
    g )
      generalized_flag=true
      ;;
    c )
      max_cycle_length=$OPTARG
      ;;
    h )
      max_chain_length=$OPTARG
      ;;
    \? )
      echo "Invalid option: -$OPTARG" 1>&2
      echo "Usage: $0 -f <filename>, -a, -g, -c <max_cycle_length>, -h <max_chain_length>, or -s <skip_pattern>"
      exit 1
      ;;
    : )
      echo "Option -$OPTARG requires an argument." 1>&2
      exit 1
      ;;
  esac
done

should_skip_file() {
    local base_file="$1"
    for pattern in "${skip_patterns[@]}"; do
        if [[ "$base_file" =~ $pattern ]]; then
            return 0
        fi
    done
    return 1
}

run_instance() {
    local input_file="$1"
    local output_file="$2"

    if [ "$generalized_flag" = true ]; then
        if [ -z "$max_cycle_length" ] || [ -z "$max_chain_length" ]; then
            echo "Error: When using the -g flag, you must specify both -c <max_cycle_length> and -h <max_chain_length>."
            exit 1
        fi
        echo "Running generalized version with max cycle length $max_cycle_length and max chain length $max_chain_length for $input_file"
        python3 "$script_dir/run.py" "$input_file" "$output_file" -g -c "$max_cycle_length" -h "$max_chain_length"
    else
        echo "Running normal version for $input_file"
        python3 "$script_dir/run.py" "$input_file" "$output_file"
    fi
}

if [ -n "$file" ]; then
    base_file=$(basename -- "$file")
    filename_no_ext="${base_file%.*}"
    output_file="$solutions_dir/${filename_no_ext}_SOLUTION.txt"
    run_instance "$instance_dir/$file" "$output_file"
elif [ "$all_flag" = true ]; then
    for file in "$instance_dir"/*.txt; do
        base_file=$(basename "$file")

        if should_skip_file "$base_file"; then
            echo "Skipping $base_file"
            continue
        fi

        filename_no_ext="${base_file%.*}"
        output_file="$solutions_dir/${filename_no_ext}_SOLUTION.txt"
        run_instance "$file" "$output_file"
    done
else
    echo "Error: No arguments provided."
    echo "Usage: $0 -f <filename>, -a, -g, -c <max_cycle_length>, -h <max_chain_length>, or -s <skip_pattern>"
    exit 1
fi
