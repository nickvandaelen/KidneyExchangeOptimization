#!/bin/bash

script_dir="$(cd "$(dirname "$0")" && pwd)"
instance_dir="$script_dir/../Instance Files"
solutions_dir="$script_dir/Solutions"

mkdir -p "$solutions_dir"

file=""
all_flag=false
skip_patterns=("A.*")

while getopts ":f:as:" opt; do
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
    \? )
      echo "Invalid option: -$OPTARG" 1>&2
      echo "Usage: $0 -f <filename>, -a, or -s <skip_pattern>"
      exit 1
      ;;
    : )
      echo "Option -$OPTARG requires an argument." 1>&2
      exit 1
      ;;
  esac
done

if [ -z "$file" ] && [ "$all_flag" = false ]; then
    echo "Error: No arguments provided."
    echo "Usage: $0 -f <filename>, -a, or -s <skip_pattern>"
    exit 1
fi

if [ -n "$file" ]; then
    filename=$(basename -- "$file")
    filename_no_ext="${filename%.*}"
    output_file="$solutions_dir/${filename_no_ext}_SOLUTION.txt"
fi

should_skip_file() {
    local base_file="$1"
    for pattern in "${skip_patterns[@]}"; do
        if [[ "$base_file" =~ $pattern ]]; then
            return 0
        fi
    done
    return 1
}

if [ "$all_flag" = true ]; then
    echo "Running all instances..."
    for file in "$instance_dir"/*.txt; do
        base_file=$(basename "$file")

        if should_skip_file "$base_file"; then
            echo "Skipping $base_file"
            continue
        fi

        filename_no_ext="${base_file%.*}"
        output_file="$solutions_dir/${filename_no_ext}_SOLUTION.txt"
        echo "Processing $file..."
        python3 "$script_dir/run.py" "$file" "$output_file"
    done
elif [ -f "$instance_dir/$file" ]; then
    echo "Processing $file ..."
    python3 "$script_dir/run.py" "$instance_dir/$file" "$output_file"
else
    echo "Error: '$file' is not a valid filename."
    echo "Usage: $0 -f <filename>, -a, or -s <skip_pattern>"
    exit 1
fi
