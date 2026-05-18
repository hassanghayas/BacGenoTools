#!/bin/bash

# Description
# This script converts a combined AMR gene file into a presence/absence matrix.
# It reads a sample ID from column 1 and gene names from a specified column,
# then outputs a matrix where rows = samples and columns = genes (1 = present, 0 = absent).

# Usage function
usage() {
    echo "Usage:"
    echo "  $0 -i <input_file> -o <output_file> -c <gene_column_number>"
    echo ""
    echo "Author: Hassan Ghayas"
    echo ""
    echo "Description:"
    echo "This script converts a combined AMR gene file into a presence/absence matrix."
    echo "It reads a sample ID from column 1 and gene names from a specified column."
    echo ""
    echo "Options:"
    echo "  -i    Input file (combined samples AMR genes)"
    echo "  -o    Output file (matrix)"
    echo "  -c    Column number containing gene names"
    echo ""
    echo "Example:"
    echo "  $0 -i amr_combined.tsv -o matrix.tsv -c 5"
    exit 1
}

# If no arguments provided, show help
if [ $# -eq 0 ]; then
    usage
fi

# Parse arguments
while getopts "i:o:c:h" opt; do
    case $opt in
        i) input="$OPTARG" ;;
        o) output="$OPTARG" ;;
        c) col_num="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# If any required argument missing, ask interactively
[ -z "$input" ] && read -p "Enter input file name: " input
[ -z "$output" ] && read -p "Enter output file name: " output
[ -z "$col_num" ] && read -p "Enter gene column number: " col_num

# Main awk logic
awk -v col="$col_num" '
    NR>1 {
        samples[$1]=1
        genes[$col]=1
        data[$1 FS $col]=1
    }
    END {
        # print header
        printf "sample"
        for (g in genes) printf "\t%s", g
        print ""

        # print rows
        for (s in samples) {
            printf "%s", s
            for (g in genes) {
                printf "\t%d", ( (s FS g) in data ? 1 : 0 )
            }
            print ""
        }
    }
' $input > $output
