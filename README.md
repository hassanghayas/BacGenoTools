# BacGenoTools

**BacGenoTools** is a lightweight and extensible toolkit for bacterial genomics analysis.

## Overview

BacGenoTools provides tools for analyzing bacterial whole-genome sequencing data with a focus on:

* Allele calling for custom MLST/Alelles [(allele_caller.py)](scripts/allele_caller.py)
* Assembly filtering based on size, coverage, contamination  [(assembly_qc_filter.pl)](scripts/assembly_qc_filter.pl)

Future updates will include additional modules for broader bacterial genomics analysis.

<!-- ## Features

* Fast and simple allele calling from genomic data
* MLST-style sequence typing
* Modular design for easy expansion
* Written in Python (with planned support for Perl utilities) -->

## Installation

Clone the repository:

```bash
git clone https://github.com/hassanghayas/BacGenoTools.git
cd BacGenoTools
```
create conda environment and install depedencies:

```bash
conda create -n bacgenotools -c conda-forge -c bioconda python=3.10 perl bwa blast
conda activate bacgenotools
```


<!-- Install dependencies (example):

```bash
pip install -r requirements.txt
``` -->

## Usage

Usage for allele calling:

```bash
python allele_caller.py --genome <genome.fasta> --alleles <alleles.fasta> --output <output.tsv>
```
Usage for assembly filter:

```bash
perl assembly_qc_filter.pl -c <assembly.fasta> -1 <read1> -2 <read2> -o <output prefix> -p <phiX.fasta/contamination.fasta>
```

<!-- ## Project Structure

```
BacGenoTools/
│── scripts/
│── README.md
``` -->

<!-- ## Roadmap

Planned features:

* [ ] SNP/variant calling
* [ ] Genome assembly integration
* [ ] Antimicrobial resistance (AMR) detection
* [ ] Quality control (QC) module
* [ ] Visualization tools
* [ ] Workflow automation

## Contributing

Contributions are welcome. You can:

* Open an issue for bugs or feature requests
* Submit pull requests
* Suggest improvements to documentation

## License

This project is licensed under the MIT License.

## Citation

If you use BacGenoTools in your research, please cite:

```
(Coming soon)
```

## Contact

For questions or collaboration:

* GitHub Issues
* Email: [your-email@example.com](mailto:your-email@example.com) -->

---

**BacGenoTools** aims to provide a simple tools for bacterial whole genome sequence analysis.
