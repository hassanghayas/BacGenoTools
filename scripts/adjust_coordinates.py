#!/usr/bin/env python3

import argparse
import sys
from Bio import SeqIO


def convert_coordinates(fasta_file, annotation_file, output_file):
    # Step 1: Calculate the cumulative offset for each contig
    contig_offsets = {}
    current_offset = 0
    
    print("Processing reference genome contigs...")
    for record in SeqIO.parse(fasta_file, "fasta"):
        # We strip weights/spaces to match typical tab-file naming conventions
        contig_id = record.id 
        contig_offsets[contig_id] = current_offset
        print(f"  {contig_id}: Global start position = {current_offset}")
        current_offset += len(record.seq)

    # Step 2: Read original annotations and shift coordinates
    print(f"\nTranslating coordinates from {annotation_file}...")
    with open(annotation_file, 'r') as infile, open(output_file, 'w') as outfile:
        # Read header if it exists
        header = infile.readline()
        # BRIG custom feature files typically expect: Start\tEnd\tLabel
        # Or you can keep your own header format:
        outfile.write("Start\tEnd\tLabel\n")
        
        for line in infile:
            if not line.strip():
                continue
            parts = line.strip().split('\t')
            
            contig_id = parts[0]
            start_pos = int(parts[1])
            end_pos = int(parts[2])
            gene_name = parts[3]
            
            if contig_id in contig_offsets:
                # Add the offset to the local coordinates
                global_start = start_pos + contig_offsets[contig_id]
                global_end = end_pos + contig_offsets[contig_id]
                
                # Write to the new BRIG-compatible file
                outfile.write(f"{global_start}\t{global_end}\t{gene_name}\n")
            else:
                print(f"Warning: Contig '{contig_id}' found in annotations but not in FASTA!")

    print(f"Success! BRIG-compatible file saved to: {output_file}")

def main():
    parser = argparse.ArgumentParser(
        description="adjust_coordinates.py v0.1"
            "\nConvert local contig coordinates to global coordinates for BRIG annotation visualization."
            "\nAuthor: Hassan Ghayas (https://github.com/hassanghayas) \n"
            ,
        add_help=False,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "-f", "--fasta", required=True,
        help="Reference genome FASTA file with contigs."
    )
    parser.add_argument(
        "-a", "--annotation", required=True,
        help="Original annotation file with local coordinates."
    )
    parser.add_argument(
        "-o", "--output", default="brig_annotation.txt",
        help="Output BRIG-compatible annotation file (default: %(default)s)."
    )

    # Re-adding help for a cleaner look
    parser.add_argument(
        '-h', '--help', action='help',
        help='show this help message'
    )

    args = parser.parse_args()

    convert_coordinates(args.fasta, args.annotation, args.output)

if __name__ == "__main__":
    main()