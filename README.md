How to Use
This script automates the process of creating a codon-based nucleotide alignment from a set of unaligned coding sequences (CDS).

Prerequisites
Before you begin, ensure you have the following installed:

Python 3

Biopython: pip install biopython

MAFFT: Must be installed and accessible in your system's PATH.

Command-Line Arguments
-i, --input: (Required) Path to the input FASTA file containing nucleotide coding sequences.

-o, --output: (Required) Path for the output codon-aligned FASTA file.

-g, --genetic-code: (Optional) The NCBI genetic code table to use for translation. Default: 1 (The Standard Code).

--mafft-args: (Optional) A quoted string of arguments to pass to MAFFT. Default: "--auto --quiet".
