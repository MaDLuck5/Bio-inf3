#!/usr/bin/env python3
"""
Assignment 1
"""
# imports
import argparse
import Bio
import sys



def main(args):
    parser = argparse.ArgumentParser(description="Calculate the severity scores for SNPs"
                                                 "in a MSA.")

    parser.add_argument("file", metavar="F", type=str,
                        help="The file path to the Multiple sequence alignment, only multi Fasta files  ")

    parser.add_argument("-l", metavar="location", type=str,
                        help="A single location for the SNP to calculate the severity")


    parser.add_argument("-s", metavar="save", type=str, choices=["amino", "alignment", "all"],
                        help="What the program wil save and write to the destination denoted "
                             "by -out, you can choose from the following options:"
                             "amino; only saves the amino translation,"
                             "alignment; only saves the alignment information,"
                             "all; saves amino & alignment information."
                             "If no option was given, will not save data and continue to only "
                             "calculate the single SNP location.")

    parser.add_argument("-out", metavar="output", type=str,
                        help="The path and filename for saving the desired information, "
                             "make sure to denote it as an .csv file.")

    args = parser.parse_args()
    return 0


if __name__ == "__main__":
    main(sys.argv)
    sys.exit()


