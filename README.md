Of course. Here is the guide formatted in Markdown to be copied directly.

PAL2NAL: A Codon-Based Alignment Script
This script automates the process of creating a codon-based nucleotide alignment from a set of unaligned coding sequences (CDS). It translates the sequences, aligns them at the protein level using MAFFT, and then back-translates the alignment to reflect the original codons.

Prerequisites
Before you begin, ensure you have the following installed:

Python 3

Biopython: Install using the command pip install biopython

MAFFT: Must be installed and accessible in your system's PATH.

Command-Line Arguments
-i, --input: (Required) Path to the input FASTA file containing nucleotide coding sequences.

-o, --output: (Required) Path for the output codon-aligned FASTA file.

-g, --genetic-code: (Optional) The NCBI genetic code table to use for translation. Default: 1 (The Standard Code).

--mafft-args: (Optional) A quoted string of arguments to pass to MAFFT. Default: "--auto --quiet".

Example Usage
First, create an example input file named sample_cds.fasta.

Sample Input File (sample_cds.fasta)
Fragmento de código

>gene_A
ATGGCTTCACTACTGAACGTGAAAGCTTCCGGCGTGGAACAGCCGTCAGGATGGCAAG
>gene_B
ATGGCAAGCGCCGGCGTGGAGCAGCCATCAGGGTGGCAATAA
>gene_C
ATGGCGAGCGCTGGAGTGGAGCAGCCGTCAGGAAGGTGGCAA
Note: The sequence gene_B includes a terminal stop codon (TAA), which the script will automatically detect and remove before translation.

Example 1: Basic Alignment
This command aligns the sequences using the standard genetic code and default MAFFT settings.

Bash

python pal2nal.py -i sample_cds.fasta -o aligned_basic.fasta -g 1
Explanation: The script reads sample_cds.fasta, translates the sequences using genetic code 1, aligns the resulting proteins, and writes the final codon alignment to aligned_basic.fasta.

Example 2: Using a Different Genetic Code
This is useful for sequences from organisms that don't use the standard code, such as vertebrate mitochondria.

Bash

python pal2nal.py -i sample_cds.fasta -o aligned_mito.fasta -g 2
Explanation: This performs the alignment using genetic code 2 (The Vertebrate Mitochondrial Code).

Example 3: Advanced Alignment with Custom MAFFT Parameters
For a potentially more accurate (but slower) alignment, you can override the default MAFFT arguments using the --mafft-args flag. This example uses the L-INS-i algorithm, which is suitable for sequences with one alignable domain.

Bash

python pal2nal.py -i sample_cds.fasta -o aligned_linsi.fasta --mafft-args "--linsi"
Explanation: The script passes the "--linsi" flag directly to MAFFT, overriding the default "--auto --quiet". This is useful for difficult alignments that may benefit from more sophisticated algorithms.
