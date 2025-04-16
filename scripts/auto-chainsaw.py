from chainsaw import chainsaw, mutual_info, unroot
import sys
import csv
from Bio import Phylo
import argparse

description = """
This script is used to generate the data required to produce Figures
2A and 3A.  The inputs were trees reconstructed using FastTree2.
Results are written to stdout in CSV format.
"""

parser = argparse.ArgumentParser(description)
parser.add_argument("infile", type=str, help="Input tree")
parser.add_argument("gene", choices=["full", "VP1"], help="Which gene/segment?")
parser.add_argument("outfile", type=argparse.FileType('w'), nargs="?", default=sys.stdout,
                        help="<output, optional> File to write output, defaults to stdout.")
args = parser.parse_args()

if args.gene == "full":
    cutoffs = [0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 
               0.01, 0.012, 0.014, 0.016, 0.018, 0.02, 0.022, 0.024, 0.026, 
               0.028, 0.03, 0.032, 0.034, 0.036, 0.038, 0.04]
else:
    cutoffs = [0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 
               0.175, 0.2, 0.225, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]

writer = csv.writer(args.outfile)
writer.writerow(["cutoff", "nsubtrees", "mutual.inf", "normalized"])
for cutoff in cutoffs:
    phy = Phylo.read(args.infile, "newick")
    unroot(phy)
    subtrees = chainsaw(phy, cutoff)
    minfo, norm_mi = mutual_info(subtrees, hema=(args.gene == 'full'))
    writer.writerow([cutoff, len(subtrees), minfo, norm_mi])
    args.outfile.flush()
