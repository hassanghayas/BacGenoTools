#!/usr/bin/env python3
import sys
from Bio import AlignIO
from Bio.Seq import MutableSeq

def purge_recombination(alignment_file, fastgear_recomp_file, output_file):
    print(f"[*] Loading alignment: {alignment_file}")
    # Load alignment as a dictionary of mutable sequences for easy editing
    alignment = AlignIO.read(alignment_file, "fasta")
    seq_dict = {record.id: MutableSeq(str(record.seq)) for record in alignment}
    
    print(f"[*] Processing fastGEAR recombinations from: {fastgear_recomp_file}")
    
    # fastGEAR columns are typically: isolate_name, start_pos, end_pos (1-indexed)
    # Note: Check your fastGEAR output file header to ensure column positions match.
    purged_count = 0
    with open(fastgear_recomp_file, 'r') as f:
        # Skip header if it exists; adjust if your file has no header
        header = f.readline() 
        
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split()
            
            # Adjust indexes based on your fastGEAR output columns
            isolate = parts[0]
            start = int(parts[1]) - 1 # Convert 1-indexed to Python 0-indexed
            end = int(parts[2])
            
            if isolate in seq_dict:
                # Mask out the recombinant region with standard gap characters '-'
                # You can change "-" to "N" if preferred by your downstream tools
                seq_dict[isolate][start:end] = "-" * (end - start)
                purged_count += 1
            else:
                print(f"[Warning] Isolate {isolate} found in fastGEAR but not in alignment.")

    print(f"[*] Masked {purged_count} recombinant blocks.")
    print(f"[*] Writing purged alignment to: {output_file}")
    
    # Save the masked alignment back to a new FASTA file
    with open(output_file, 'w') as out_f:
        for record in alignment:
            purged_seq = seq_dict[record.id]
            out_f.write(f">{record.id}\n{purged_seq}\n")
            
    print("[+] Done!")

if __name__ == "__main__":
    # Example usage format: 
    # python purge.py core_genome_alignment.aln recombinations.txt purged_alignment.aln
    if len(sys.argv) < 4:
        print("Usage: python purge_recombination.py <alignment.aln> <recombinations.txt> <output.aln>")
        print("Description: this script purges recombinant regions from a core genome alignment based on fastGEAR output.")
        sys.exit(1)
        
    purge_recombination(sys.argv[1], sys.argv[2], sys.argv[3])