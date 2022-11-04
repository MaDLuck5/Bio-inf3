Author: Mats Slik
Date: 04-11-2022
Version: 1.0


Name
Smoking gun, an assigment to determine if an SNP in a gen sequence is deleterious or not.


Description


Installation
Make sure python3.7 or higher is installed on the pc.
and ClusteloW2 is installed: source: "http://www.clustal.org/download/current/"


Usage
With the terminal, go towards the directory where the .py script resides and execute the following command to make it executable:
chmod +x SNP_main.py

the programs takes the following arguments optional arguments:
  
  -h, --help            show this help message and exit
  -s SEQUENCE, --sequence SEQUENCE
                        the path to the coding Amino acid sequence fasta file, Please specify the path in this format'<path>'
  -f FILE, --file FILE  The file path to the Multiple sequence alignment, only multi Fasta files. Please specify the path in '<path>'
  -l LOCATION, --location LOCATION
                        A single location for the SNP to calculate the severity
  -r REPLACEMENT, --replacement REPLACEMENT
                        the Amino acid that replaces the on the location given
optional arguments:
  -cw CLUSTALWLOCATION, --clustalwlocation CLUSTALWLOCATION
                        The path to the installed Clustal W2 program .exe
  -p PERCENTAGE_DISP, --percentage_disp PERCENTAGE_DISP
                        if TRUE is given, programs shows the conservation of each amino acid in the MSA in percentages
  -save CREATE_NEW_MSA, --create_new_msa CREATE_NEW_MSA
                        if True is given will creat new MSA with the mutated sequence in the out put folder


After which, you can run the following example code to get the scrip to run:
assignment_1.py -s <path to mRNA sequence> -f <path to protein multifasta> -l <position of SNP> -r <Nuc change to>

EXAMPLE:
assignment_1.py -s "C:\Users\matsp\Documents\Bio-Infromatiscs-jaar3\Bio-inf3\assignment1\sequence2.txt" -f "C:\Users\matsp\Documents\Bio-Infromatiscs-jaar3\Bio-inf3\assignment1\homologene2.fasta" -l 32 -r G -p -save


IF this doesn't work, it is always possible to run the script the following way in the terminal:
 python .\assignment_1.py -s "C:\Users\matsp\Documents\Bio-Infromatiscs-jaar3\Bio-inf3\assignment1\sequence2.txt" -f "C:\Users\matsp\Documents\Bio-Infromatiscs-jaar3\Bio-inf3\assignment1\homologene2.fasta" -l 32 -r G -p -save

If you need help with the script, you can always type the following command to get help for the program:
assignment_1.py -h


Support
m.p.slik@st.hanze.nl


Authors and acknowledgment
Thanks to the developers of python, biopython, the people at NCBI & The Hanze Hogeschool Groningen Bio-informatics staff.