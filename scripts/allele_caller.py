#!/usr/bin/env python3

import argparse
import subprocess
from collections import defaultdict
import time



def run_blast(alleles, genome, threads):
    cmd = [
        "blastn",
        "-query", alleles,
        "-subject", genome,
        "-task", "megablast",
        "-num_threads", str(threads),
        "-outfmt", "6 qseqid sseqid pident length qlen qcovs bitscore"
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        raise RuntimeError("BLAST failed")

    return result.stdout


def parse_blast(blast_output, min_id, min_cov):

    best_hits = {}

    for line in blast_output.strip().split("\n"):
        if not line:
            continue

        qseqid, sseqid, pident, length, qlen, qcovs, bitscore = line.split()

        pident = float(pident)
        qcovs = float(qcovs)
        bitscore = float(bitscore)
        hit = {
            "allele": qseqid,
            "identity": pident,
            "coverage": qcovs,
            "bitscore": bitscore,
            "contig": sseqid
        }
        locus = qseqid.split("_")[0]

        if pident < min_id or qcovs < min_cov:
            continue

        if locus not in best_hits:
            best_hits[locus] = hit
        else:
            old = best_hits[locus]

            if (
                qcovs > old["coverage"] or
                (qcovs == old["coverage"] and pident > old["identity"]) or
                (qcovs == old["coverage"] and pident == old["identity"] and bitscore > old["bitscore"])
            ):
                best_hits[locus] = hit

    return best_hits


def get_loci(alleles_fasta):
    loci = set()

    with open(alleles_fasta) as f:
        for line in f:
            if line.startswith(">"):
                allele = line.strip()[1:]
                locus = allele.split("_")[0]
                loci.add(locus)

    return sorted(loci)


def write_output(best_hits, loci, output):

    present = 0
    absent = 0

    with open(output, "w") as out:

        out.write("Locus\tAllele\tIdentity\tCoverage\tContig\n")

        for locus in loci:

            if locus in best_hits:

                hit = best_hits[locus]

                out.write(
                    f"{locus}\t{hit['allele']}\t"
                    f"{hit['identity']:.2f}\t"
                    f"{hit['coverage']:.2f}\t"
                    f"{hit['contig']}\n"
                )

                present += 1

            else:

                out.write(f"{locus}\tAbsent\t0\t0\t-\n")

                absent += 1

    return present, absent


def main():
    start_time = time.time()

    parser = argparse.ArgumentParser(
        description="allele_caller.py v0.1"
                    "\nAuthor: Hassan Ghayas (https://github.com/hassanghayas) \n"
                    "\nAllele caller for custom MLST. "
                    "This script uses BLAST to identify best allele matches.",
        add_help=False,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Creating Custom Groups
    required = parser.add_argument_group('Required')
    optional = parser.add_argument_group('Optional')

    # Adding required arguments to the 'required' group
    required.add_argument("-g", "--genome", required=True,
                        help="Genome assembly fasta")

    required.add_argument("-a", "--alleles", required=True,
                        help="Allele fasta database")

    # Adding optional arguments to the 'optional' group
    optional.add_argument("-o", "--output", default="output.tsv",
                        help="name of output table (default: %(default)s)")

    optional.add_argument("--threads", type=int, default=1,
                        help="number of threads for BLAST (default: %(default)s)")

    optional.add_argument("--min_id", type=float, default=90,
                        help="Minimum percent identity (default: %(default)s)")

    optional.add_argument("--min_cov", type=float, default=90,
                        help="Minimum coverage (default: %(default)s)")

    # Re-adding help to the optional group for a cleaner look
    optional.add_argument('-h', '--help', action='help', help='show this help message')

    args = parser.parse_args()

    print("Reading loci...")
    loci = get_loci(args.alleles)

    print("Running BLAST...")
    blast_output = run_blast(args.alleles, args.genome, args.threads)

    print("Parsing BLAST hits...")
    best_hits = parse_blast(blast_output, args.min_id, args.min_cov)

    print("Writing results...")
    present, absent = write_output(best_hits, loci, args.output)

    print(f"\nTotal loci: {len(loci)}")
    print(f"Present: {present}")
    print(f"Absent: {absent}")

    end_time = time.time()

    runtime = end_time - start_time
    print(f"Total runtime: {runtime:.2f} seconds")

if __name__ == "__main__":
    main()
