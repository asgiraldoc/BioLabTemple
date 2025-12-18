#!/usr/bin/env python3

import argparse
import sys
import tempfile
import os
from subprocess import run, PIPE

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data import CodonTable

# --- Argument Parsing ---
def get_args():
    """
    Parses command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="A script to translate CDS, align proteins with MAFFT, and back-translate to a codon alignment.",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input FASTA file with nucleotide coding sequences (CDS)."
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output file for the codon-aligned nucleotide sequences."
    )
    parser.add_argument(
        "-g", "--genetic-code",
        type=int,
        default=1,
        help="NCBI genetic code table for translation. Common codes:\n"
             "1: Standard Code\n"
             "2: Vertebrate Mitochondrial\n"
             "5: Invertebrate Mitochondrial\n"
             "11: Bacterial, Archaeal and Plant Plastid\n"
             "(Default: 1)"
    )
    parser.add_argument(
        "--mafft-args",
        type=str,
        default="--auto --quiet",  # <-- CAMBIO REALIZADO AQUÍ
        help="Additional arguments to pass to MAFFT. (Default: '--auto --quiet')"
    )

    return parser.parse_args()

# --- Core Functions ---

def validate_and_translate_sequences(records, genetic_code_id):
    """
    Validates sequences for stop codons, translates them, and prepares them for alignment.
    - Removes terminal stop codons automatically.
    - Halts if an internal stop codon is found.
    """
    protein_records = []
    original_nucleotide_records = []

    try:
        codon_table = CodonTable.generic_by_id[genetic_code_id]
        stop_codons = set(codon_table.stop_codons)
    except KeyError:
        sys.exit(f"Error: Genetic code table {genetic_code_id} is not valid.")

    print("Validating sequences and checking for stop codons...", file=sys.stderr)

    for rec in records:
        seq_str = str(rec.seq).upper()

        # Check for terminal stop codon
        if len(seq_str) >= 3 and seq_str[-3:] in stop_codons:
            print(f"INFO: Terminal stop codon '{seq_str[-3:]}' found and removed from sequence '{rec.id}'.", file=sys.stderr)
            seq_str = seq_str[:-3]
            rec.seq = Seq(seq_str) # Update sequence record

        # Check for internal stop codons
        for i in range(0, len(seq_str), 3):
            codon = seq_str[i:i+3]
            if len(codon) == 3 and codon in stop_codons:
                sys.exit(f"ERROR: Internal stop codon '{codon}' found in sequence '{rec.id}' at position {i+1}. Halting execution.")

        # Final length check before translation
        if len(rec.seq) % 3 != 0:
            sys.exit(f"Error: Sequence '{rec.id}' has a post-validation length of {len(rec.seq)}, which is not a multiple of 3. Cannot translate.")

        # Translate the cleaned sequence
        protein_seq = rec.seq.translate(table=genetic_code_id, cds=False) # cds=False because we already checked stops
        protein_records.append(SeqRecord(protein_seq, id=rec.id, description=""))
        original_nucleotide_records.append(rec) # Keep the cleaned nucleotide version

    return protein_records, original_nucleotide_records

def align_with_mafft(protein_records, mafft_args):
    """
    Aligns protein sequences using MAFFT via a temporary file.
    """
    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix=".fasta") as temp_in:
        SeqIO.write(protein_records, temp_in, 'fasta')
        temp_in_name = temp_in.name

    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix=".fasta") as temp_out:
        temp_out_name = temp_out.name

    print("Running MAFFT for protein alignment...", file=sys.stderr)

    command = ['mafft'] + mafft_args.split() + [temp_in_name]

    try:
        with open(temp_out_name, 'w') as f_out:
            process = run(command, check=True, stdout=f_out, stderr=PIPE, universal_newlines=True)

        if process.stderr:
            # Solo muestra mensajes de MAFFT si son errores, ya que --quiet suprime el progreso
            print("MAFFT messages (stderr):\n" + process.stderr, file=sys.stderr)

    except FileNotFoundError:
        sys.exit("Error: 'mafft' command not found. Is MAFFT installed and in your PATH?")
    except Exception as e:
        sys.exit(f"An error occurred during MAFFT execution: {e}")

    aligned_protein_records = list(SeqIO.parse(temp_out_name, 'fasta'))

    os.remove(temp_in_name)
    os.remove(temp_out_name)

    return aligned_protein_records

def back_translate(aligned_proteins, original_nucleotides):
    """
    Back-translates a protein alignment to a codon alignment.
    """
    print("Back-translating protein alignment to codons...", file=sys.stderr)

    nuc_dict = {rec.id: rec.seq for rec in original_nucleotides}
    codon_aligned_records = []

    for prot_rec in aligned_proteins:
        if prot_rec.id not in nuc_dict:
            print(f"Warning: ID '{prot_rec.id}' from alignment not found in original nucleotide file. Skipping.", file=sys.stderr)
            continue

        original_nuc_seq = nuc_dict[prot_rec.id]
        codon_aligned_seq = ""
        nuc_pos = 0

        for amino_acid in prot_rec.seq:
            if amino_acid == '-':
                codon_aligned_seq += '---'
            else:
                codon = original_nuc_seq[nuc_pos:nuc_pos+3]
                codon_aligned_seq += codon
                nuc_pos += 3

        codon_aligned_records.append(SeqRecord(Seq(codon_aligned_seq), id=prot_rec.id, description=""))

    return codon_aligned_records

# --- Main Execution Block ---
def main():
    """
    Main function to orchestrate the pipeline.
    """
    args = get_args()

    try:
        print(f"Reading nucleotide sequences from: {args.input}", file=sys.stderr)
        initial_nuc_records = list(SeqIO.parse(args.input, 'fasta'))
        if not initial_nuc_records:
            sys.exit(f"Error: No sequences found in {args.input}. Is the file empty or in the wrong format?")
    except FileNotFoundError:
        sys.exit(f"Error: Input file not found at {args.input}")
    except Exception as e:
        sys.exit(f"An error occurred reading the input file: {e}")

    # 2. Validate, clean, and translate sequences
    protein_records, cleaned_nuc_records = validate_and_translate_sequences(initial_nuc_records, args.genetic_code)

    # 3. Align proteins with MAFFT
    aligned_protein_records = align_with_mafft(protein_records, args.mafft_args)

    # 4. Back-translate to codon alignment
    codon_alignment = back_translate(aligned_protein_records, cleaned_nuc_records)

    # 5. Write final nucleotide alignment to output
    try:
        SeqIO.write(codon_alignment, args.output, 'fasta')
        print(f"\n✅ Success! Codon alignment saved to: {args.output}", file=sys.stderr)
    except Exception as e:
        sys.exit(f"An error occurred writing the output file: {e}")


if __name__ == '__main__':
    main()
