#!/usr/bin/env python3
"""
Assignment 1
"""
# imports
import argparse
from Bio import Seq, SeqIO, AlignIO
import sys
from Bio.Align.Applications import ClustalOmegaCommandline


def alignment_score(msa_alignment, nuc_seq):
    scoring = list()
    msa = AlignIO.read(msa_alignment, "fasta")

    print(msa)

    for i in range(msa.get_alignment_length()):
        temporary_dict = {}

        unique_aa = set(msa[:, i])
        for amino_acid in unique_aa:
            temporary_dict[amino_acid] = alignment[:, i].count(amino_acid) / len(alignment[:, i])
        temporary_dict = {x: a for x, a in sorted(temporary_dict.items(), key=lambda item: item[1], reverse=True)}
        scoring.append(temporary_dict)

    return alignment


def alignment(in_file, out_file):
    """

    :param in_file:
    :param out_file:
    :return:
    """
    clustalomega = ClustalOmegaCommandline(infile=in_file, outfile=out_file, verbose=True, auto=False,
                                           align=True, output="FASTA")
    clustalomega()


def nuc_translator_to_aa(nuc_seq):
    """

    :param nuc_seq:
    :return:
    """
    aa_seq = SeqIO.read(nuc_seq, "fasta")
    protein_seq = Seq.translate(aa_seq)

    return protein_seq


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
    print(args.location)

    return 0


if __name__ == "__main__":
    main()
    sys.exit()
