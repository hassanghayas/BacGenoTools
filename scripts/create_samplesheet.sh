#!/bin/bash
set -euo pipefail

usage() {
    echo "Usage: create_samplesheet.sh <path/to/directory>"
}

# to see if the argument for directory is provided
if [ "$#" -ne 1 ]; then
    usage
    echo "Error: no file directory provided!! Please provide the fastq file directory"
    exit 1
fi

# to see if directory is provided
if [ ! -d "$1" ]; then
    usage
    echo "Error: $1 is not a directory! Please provide the fastq file directory"
    exit 1
fi

# to check if provided directory is empty
if [ -z "$(find "$1" -mindepth 1 -print -quit)" ]; then
    echo "Error: Directory '$1' is empty."
    exit 1
fi

# to check if the directory contains .fastq.gz
if ! find "$1" -maxdepth 1 -type f -name "*.fastq.gz" | grep -q .; then
    echo "Error: Directory '$1' does not contain fastq.gz"
    exit 1
fi

# The first argument is the base directory where fastq is located
fastq_directory=$(basename $1)

echo 'sample,R1,R2'
for id in $(basename -a $fastq_directory/*.fastq.gz |cut -d_ -f1| uniq); do
    R1=$(find $PWD/$fastq_directory/${id}_*R1*fastq.gz)
    R2=$(find $PWD/$fastq_directory/${id}_*R2*fastq.gz)

    echo "$id,$R1,$R2"
done
