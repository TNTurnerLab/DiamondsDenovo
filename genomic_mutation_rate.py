#!/bin/python3

import argparse
import subprocess
import tempfile
from pyliftover import LiftOver
from Bio import SeqIO
import numpy as np
import pyBigWig

def get_sequence_from_fasta(fasta_dict, chrom, start, end):
    start -= 1
    if chrom in fasta_dict:
        print (f"The chromosome is: {chrom}")
        sequence = fasta_dict[chrom].seq[start:end]
        return str(sequence).upper()
        print(f"Fetched sequence: {sequence[:50]}...{sequence[-50:]}")
    return None

def run_mafft(human_seq, other_seq):
    with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix=".fasta") as fasta_file:
        fasta_file.write(">human\n")
        fasta_file.write(human_seq + "\n")
        fasta_file.write(">other\n")
        fasta_file.write(other_seq + "\n")
        fasta_file_name = fasta_file.name

    output_file_name = tempfile.mktemp(suffix=".fasta")

    mafft_command = [
        "mafft",
        "--auto",
        fasta_file_name
    ]
    with open(output_file_name, "w") as output_file:
        subprocess.run(mafft_command, stdout=output_file)

    aligned_seqs = list(SeqIO.parse(output_file_name, "fasta"))
    human_aligned_seq = str(aligned_seqs[0].seq)
    other_aligned_seq = str(aligned_seqs[1].seq)

    return human_aligned_seq, other_aligned_seq

def count_mutations(seq1, seq2):
    mutations = 0
    gaps_seq1 = 0
    gaps_seq2 = 0
    for a, b in zip(seq1, seq2):
        if a != b:
            if a == '-':
                gaps_seq1 += 1
            elif b == '-':
                gaps_seq2 += 1
            else:
                mutations += 1
    total_mutations = mutations + gaps_seq1 + gaps_seq2
    return total_mutations

def calculate_gc_content(seq):
    return (seq.count('G') + seq.count('C')) / len(seq)

def get_cadd_scores(bw_file, chrom, start, end):
    with pyBigWig.open(bw_file) as bw:
        scores = bw.values(chrom, start, end)
    return scores

def display_alignment(human_seq, other_seq):
    alignment_display = []
    for a, b in zip(human_seq, other_seq):
        if a == b:
            alignment_display.append('|')
        else:
            alignment_display.append('*' if a != '-' and b != '-' else ' ')
    return ''.join(alignment_display)

def process_region(region, lo, human_fasta_dict, other_fasta_dict, cadd_file, divergence_time, gc_reference):
    chrom, start, end = region[0], int(region[1]), int(region[2])
    print(f"Processing region: {chrom}:{start}-{end}")

    other_coords = lo.convert_coordinate(chrom, start, '+')
    if not other_coords:
        print(f"Could not lift over {chrom}:{start}-{end}")
        return None

    other_chrom = other_coords[0][0]
    other_start = int(other_coords[0][1])
    other_end = other_start + (end - start)
    print(f"{other_coords}")

    print(f"Liftover {chrom}:{start}-{end} to other {other_coords[0][0]}:{other_start}-{other_end}")

    human_seq = get_sequence_from_fasta(human_fasta_dict, chrom, start, end)
    other_seq = get_sequence_from_fasta(other_fasta_dict, other_chrom, other_start, other_end)

    if not human_seq or not other_seq:
        print(f"Could not get sequences for {chrom}:{start}-{end}")
        return None

    human_aligned_seq, other_aligned_seq = run_mafft(human_seq, other_seq)

    mutations = count_mutations(human_aligned_seq, other_aligned_seq)
    aligned_length = len(human_aligned_seq)
    mutation_rate = mutations / (aligned_length * divergence_time)

    gc_human = calculate_gc_content(human_seq)
    gc_other = calculate_gc_content(other_seq)
    average_gc_content = (gc_human + gc_other) / 2.0
    gc_adjusted_mutation_rate = mutation_rate * (1 + abs(average_gc_content - gc_reference) / gc_reference)

    cadd_scores = get_cadd_scores(cadd_file, chrom, start, end)
    average_conservation_score = np.mean(cadd_scores)
    constrained_mutation_rate = mutation_rate * (1/average_conservation_score)

    return (chrom, start, end, *region[3:], average_gc_content, mutations, aligned_length, average_conservation_score, mutation_rate, gc_adjusted_mutation_rate, constrained_mutation_rate)

def main(bed_file, chain_file, human_fasta_path, other_fasta_path, divergence_time, cadd_file, output_file, combined_output_file=None):

    lo = LiftOver(chain_file)

    human_fasta_dict = SeqIO.to_dict(SeqIO.parse(human_fasta_path, "fasta"))
    other_fasta_dict = SeqIO.to_dict(SeqIO.parse(other_fasta_path, "fasta"))

    with open(bed_file) as f:
        regions = [line.strip().split() for line in f.readlines()]

    results = []

    total_mutations = 0
    total_aligned_length = 0
    total_gc_content = 0
    total_conservation_score = 0

    divergence_time = args.divergence_time
    print(f"Divergence Time: {divergence_time}")
 
    gc_reference = 0.42

    for region in regions:
        result = process_region(region, lo, human_fasta_dict, other_fasta_dict, cadd_file, divergence_time, gc_reference)
        if result:
            chrom, start, end, *bed_columns, average_gc_content, mutations, aligned_length, average_conservation_score, mutation_rate, gc_adjusted_mutation_rate, constrained_mutation_rate = result
            results.append((chrom, start, end, *bed_columns, average_gc_content, mutations, aligned_length, average_conservation_score, mutation_rate, gc_adjusted_mutation_rate, constrained_mutation_rate))

            total_mutations += mutations
            total_aligned_length += aligned_length
            total_gc_content += average_gc_content * aligned_length
            total_conservation_score += average_conservation_score * aligned_length

    combined_gc_content = total_gc_content / total_aligned_length
    combined_conservation_score = total_conservation_score / total_aligned_length
    combined_mutation_rate = total_mutations / (total_aligned_length * divergence_time)
    combined_gc_adjusted_mutation_rate = combined_mutation_rate * (1 + abs(combined_gc_content - gc_reference) / gc_reference)
    combined_constrained_mutation_rate = combined_mutation_rate * (1/combined_conservation_score)

    combined_headers = ["combined_gc_content", "total_mutations", "total_aligned_length", "combined_conservation_score", "combined_mutation_rate", "combined_gc_adjusted_mutation_rate", "combined_constrained_mutation_rate"]

    bed_column_headers = [f"inputBedColumn{i+1}" for i in range(3, len(regions[0]))]
    headers = ["chrom", "start", "end"] + bed_column_headers + ["average_gc_content", "mutations", "aligned_length", "average_conservation_score", "mutation_rate", "gc_adjusted_mutation_rate", "constrained_mutation_rate"]
    with open(output_file, "w") as f:
        f.write("\t".join(headers) + "\n")
        for result in results:
            f.write("\t".join(map(str, result)) + "\n")

    if combined_output_file:
        combined_results = [combined_gc_content, total_mutations, total_aligned_length, combined_conservation_score, combined_mutation_rate, combined_gc_adjusted_mutation_rate, combined_constrained_mutation_rate]
        
        with open(combined_output_file, "w") as f:
            f.write("\t".join(combined_headers) + "\n")
            f.write("\t".join(map(str, combined_results)) + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate mutation rates between human and other genomic regions")
    parser.add_argument("-b", "--bed_file", required=True, help="Input BED file with human genomic coordinates")
    parser.add_argument("-c", "--chain_file", required=True, help="Chain file for LiftOver")
    parser.add_argument("-hf", "--human_fasta", required=True, help="Human reference genome FASTA file")
    parser.add_argument("-cf", "--other_fasta", required=True, help="Other reference genome FASTA file")
    parser.add_argument("-d", "--divergence_time", required=True, type=float, help="Divergence time")
    parser.add_argument("-p", "--cadd_file", required=True, help="CADD bigWig file")
    parser.add_argument("-o", "--output_file", required=True, help="Output file to save mutation rates")
    parser.add_argument("--combined_output_file", help="Optional output file to save combined metrics", required=False)
    args = parser.parse_args()
    main(args.bed_file, args.chain_file, args.human_fasta, args.other_fasta, args.divergence_time, args.cadd_file, args.output_file, args.combined_output_file)