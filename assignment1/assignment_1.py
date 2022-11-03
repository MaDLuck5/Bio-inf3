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


def alignment_score(msa_alignment):
    """

    :param msa_alignment:
    :return:
    """
    scoring = list()
    msa = AlignIO.read(msa_alignment, "fasta")

    print(msa)

    for i in range(msa.get_alignment_length()):
        temporary_dict = dict()

        unique_aa = set(msa[:, i])
        for amino_acid in unique_aa:
            temporary_dict[amino_acid] = msa[:, i].count(amino_acid) / len(msa[:, i])
        temporary_dict = {x: a for x, a in sorted(temporary_dict.items(), key=lambda item: item[1], reverse=True)}
        scoring.append(temporary_dict)

    return scoring


def generate_cost(mutated_aa_seq, scoring):
    """

    :param mutated_aa_seq:
    :param scoring:
    :return:
    """
    cost = 0
    position_cost = list()

    for aa, d in zip(mutated_aa_seq, scoring):
        keys = list(d.keys())

        pre_cost = cost

        if aa != keys[0]:
            cost += d[keys[0]]
            if len(d) > 1:
                for key in keys:
                    if aa == key:
                        cost -= d[key]

        position_cost.append(cost - pre_cost)

    for i, e in enumerate(position_cost):
        print("Position {} has a cost of {}".format(i, e))
    print("The total cost for the mutated sequence is: {}".format(cost))


def alignment(in_file, out_file_dir, num):
    """
    Performs alignment of a multifasta amino acid file, creates a output file with a given indexing number.
    :param num: indexing number for output file
    :param out_file_dir:
    :param in_file: path to input multi fast file
    :return:out_file: path to the output file
    """

    out_file = out_file_dir + f"\\MSA{num}.fasta"
    clustalw_exe = r"C:\Program Files (x86)\ClustalW2\clustalw2.exe"
    clustalw_cline = ClustalwCommandline(clustalw_exe, infile=in_file, outfile=out_file,
                                         align=True, output="fasta")
    clustalw_cline()

    return out_file


def nuc_translator_to_aa(nuc_seq):
    """
    Translates a Nucleotide sequence to a amino acid sequence
    :param nuc_seq: Biopython Seq object
    :return: protein_seq: Biopython Seq object
    """

    protein_seq = Seq.translate(nuc_seq)

    return protein_seq


def single_point_mutator(nuc_fasta, point, mutation):
    """
    Takes a Biopython Seq object and inserts a SNP mutation at the specified location
    :param nuc_fasta: Biopython Seq object
    :param point: point in the sequence for the mutation
    :param mutation: the mutation wanted
    :return:mutated_aa_seq: Biopython Seq object
    """
    temp_nuc_seq = list(str(seq_reader(nuc_fasta)))
    temp_nuc_seq[(point - 1)] = mutation
    mutated_aa_seq = Seq("".join(temp_nuc_seq))
    return mutated_aa_seq


def seq_reader(path_to_file):
    """
    sequence reader, reads a single fasta styled sequence and returns a Biopython Seq object
    :param path_to_file:
    :return: sequence: a  Biopython Seq object
    """
    fasta_object = SeqIO.read(path_to_file, "fasta")
    sequence = fasta_object.seq
    return sequence


def command_line_parsing():
    """
    command line parser
    :return:
    """
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
    """
    Main function
    :return:
    """
    args = command_line_parsing()
    working_dir = os.getcwd()
    output_file_dir = working_dir + '\\output'

    msa_01 = alignment(args.file, output_file_dir, 1)

    aa_seq = nuc_translator_to_aa(seq_reader(args.sequence))
    mutated_aa_seq = nuc_translator_to_aa(single_point_mutator(args.sequence, args.location, args.replacement))

    scoring = alignment_score(msa_01)

    generate_cost(aa_seq, scoring)

    return 0


if __name__ == "__main__":
    main()
    sys.exit()
