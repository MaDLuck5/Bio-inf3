#!/usr/bin/env python3
"""
Assignment 1
"""
# imports
import os
import argparse
from Bio import Seq, SeqIO, AlignIO
from Bio.Seq import Seq
import sys
from Bio.Align.Applications import ClustalwCommandline


def alignment_score(msa_alignment, nuc_fasta, location, change):
    scoring = list()
    msa = AlignIO.read(msa_alignment, "fasta")

    print(msa)

    mutated_aa_seq = nuc_translator_to_aa(single_point_mutator(nuc_fasta, location, change))

    for i in range(msa.get_alignment_length()):
        temporary_dict = {}

        unique_aa = set(msa[:, i])
        for amino_acid in unique_aa:
            temporary_dict[amino_acid] = alignment[:, i].count(amino_acid) / len(alignment[:, i])
        temporary_dict = {x: a for x, a in sorted(temporary_dict.items(), key=lambda item: item[1], reverse=True)}
        scoring.append(temporary_dict)

    return alignment


def alignment(in_file, out_file_dir, num):
    """

    :param num:
    :param out_file_dir:
    :param in_file:
    :return:
    """

    out_file = out_file_dir + f"\\MSA{num}.fasta"
    clustalw_exe = r"C:\Program Files (x86)\ClustalW2\clustalw2.exe"
    clustalw_cline = ClustalwCommandline(clustalw_exe, infile=in_file, outfile=out_file,
                                         align=True, output="fasta")
    clustalw_cline()

    return out_file


def nuc_translator_to_aa(nuc_seq):
    """

    :param nuc_seq:
    :return:
    """

    protein_seq = Seq.translate(nuc_seq)

    return protein_seq


def single_point_mutator(nuc_fasta, point, mutation):
    nuc_sequence = SeqIO.read(nuc_fasta, "fasta")
    temp_nuc_seq = list(str(nuc_sequence.seq))
    temp_nuc_seq[(point - 1)] = mutation
    mutated_aa_seq = Seq("".join(temp_nuc_seq))
    return mutated_aa_seq


def command_line_parsing():
    parser = argparse.ArgumentParser(description="Calculate the severity scores for SNPs"
                                                 "in a MSA.")

    parser.add_argument("-s", "--sequence", type=str,
                        help="the path to the coding Amino acid sequence fasta file, Please "
                             "specify the path in this format'<path>' ")

    parser.add_argument("-f", "--file", type=str,
                        help="The file path to the Multiple sequence alignment, only multi Fasta files. Please "
                             "specify the path in '<path>'")

    parser.add_argument("-l", "--location", type=int,
                        help="A single location for the SNP to calculate the severity")

    parser.add_argument("-r", "--replacement", type=str,
                        help="the Amino acid that replaces the on the location given")

    args = parser.parse_args()
    return args


def main():
    args = command_line_parsing()
    working_dir = os.getcwd()
    output_file_dir = working_dir + '\\output'

    msa_01 = alignment(args.file, output_file_dir, 1)
    alignment_score(msa_01, args.sequence, args.location, args.replacement)

    return 0


if __name__ == "__main__":
    main()
    sys.exit()
