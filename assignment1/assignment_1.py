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
from typing import List


def alignment_score(msa_alignment, flag=False):
    """

    :param flag:
    :param msa_alignment:
    :return:
    """
    scoring = list()
    msa = AlignIO.read(msa_alignment, "fasta")

    # print(msa)

    for i in range(msa.get_alignment_length()):
        temporary_dict = dict()

        unique_aa = set(msa[:, i])
        for amino_acid in unique_aa:

            percentage = msa[:, i].count(amino_acid) / len(msa[:, i])
            if flag:
                print(f"position:{i} amino acid:{amino_acid} is conserveted in {int(percentage * 100)}%")
            temporary_dict[amino_acid] = percentage
        temporary_dict = {x: a for x, a in sorted(temporary_dict.items(), key=lambda item: item[1], reverse=True)}
        scoring.append(temporary_dict)

    return scoring


def generate_cost(mutated_aa_seq, scoring):
    """

    :param mutated_aa_seq:
    :param scoring:
    :return:
    """
    counter = 0
    cost = 0
    position_cost = list()
    change_dict = {}
    for aa, d in zip(mutated_aa_seq, scoring):
        counter += 1
        keys = list(d.keys())

        pre_cost = cost

        if aa != keys[0]:
            cost += d[keys[0]]
            if len(d) > 1:
                for key in keys:
                    if aa not in change_dict:
                        change_dict[aa] = [key, round(d[key]*100, ndigits=2)]
                    else:
                        change_dict[aa].append(key)
                        change_dict[aa].append(round(d[key]*100, ndigits=2))
                    if aa == key:
                        cost -= d[key]
            else:
                if keys[0] not in change_dict:
                    change_dict[aa] = [keys[0], d[keys[0]] * 100]
                else:
                    change_dict[aa] += [keys[0], d[keys[0]] * 100]

        position_cost.append(cost - pre_cost)

    # for i, e in enumerate(position_cost):
    # print("Position {} has a cost of {}".format(i, e))
    print(f"lengthe of given mutated Sequence ={len(mutated_aa_seq)}")
    print(f"The total cost for the mutated sequence is: {round(cost, ndigits=2)}")
    if cost >= 0.8:
        print(f"cost more than 1 a Amino acid has been changed in a highly conserved place")
        for key in change_dict:
            print(f"the new Amino Acid: {key}")
            print(f"The changed Amino Acid and conservation in % :{change_dict[key]}")


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


def single_point_mutator(sequence, point, mutation):
    """
    Takes a Biopython Seq object and inserts a SNP mutation at the specified location
    :param sequence: Biopython Seq object
    :param point: point in the sequence for the mutation
    :param mutation: the mutation wanted
    :return:mutated_aa_seq: Biopython Seq object
    """
    temp_nuc_seq = list(str(sequence))
    temp_nuc_seq[(point - 1)] = mutation
    mutated_aa_seq = Seq("".join(temp_nuc_seq))
    return mutated_aa_seq


def fasta_reader(path_to_file):
    """
    sequence reader, reads a single fasta styled sequence and returns a Biopython Seq object
    :param path_to_file:
    :return: sequence: a  Biopython Seq object
    """
    fasta_object = SeqIO.read(path_to_file, "fasta")

    return fasta_object


def fasta_altere(fasta_object, point, mutation):
    """

    :param fasta_object:
    :param point:
    :param mutation:
    :return:
    """
    sequence = fasta_object.seq
    fasta_object.seq = single_point_mutator(sequence, point, mutation)
    return fasta_object


def fasta_writer(fasta_list, output_dir):
    """

    :param fasta_list:
    :param output_dir:
    :return:
    """
    new_file_path = f'{output_dir}\\new_multi.fasta'
    SeqIO.write(fasta_list, new_file_path, "fasta")

    return new_file_path


def multi_fasta_parser(fasta_path):
    """

    :param fasta_path:
    :return:
    """
    temp_list = []
    fasta_sequences = SeqIO.parse(open(fasta_path), 'fasta')
    for item in fasta_sequences:
        temp_list.append(item)

    return temp_list


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

    parser.add_argument("-p", "--percentage_disp", type=bool, default=False,
                        help="if TRUE is given, programs shows the conservation of each amino acid in the MSA in "
                             "percentages")

    parser.add_argument("-save", "--create_new_msa", type=bool, default=False,
                        help="if True is given will creat new MSA with the mutated sequence in the out put folder")

    args = parser.parse_args()
    return args


def main():
    """
    Main function
    :return:
    """
    fasta_list = []

    args = command_line_parsing()
    working_dir = os.getcwd()
    output_file_dir = working_dir + '\\output'

    mutated_fasta = fasta_altere(fasta_reader(args.sequence), args.location, args.replacement)
    mutated_fasta.seq = nuc_translator_to_aa(mutated_fasta.seq)

    if args.create_new_msa:

        fasta_list.append(mutated_fasta)
        fasta_list.extend(multi_fasta_parser(args.file))
        new_multifasta_path = fasta_writer(fasta_list, output_file_dir)
        alignment(args.file, output_file_dir, " new+snp")

    msa_02 = alignment(args.file, output_file_dir, 2)

    scoring2 = alignment_score(msa_02, args.percentage_disp)

    generate_cost(mutated_fasta.seq, scoring2)

    return 0


if __name__ == "__main__":
    main()
    sys.exit()
