from chainsaw import chainsaw, mutual_info, unroot
import sys
import csv
from Bio import Phylo
import argparse
import ipdb

description = """
This script is used to generate the data required to produce Figures
2A and 3A.  The inputs were trees reconstructed using FastTree2.
Results are written to stdout in CSV format.
"""

parser = argparse.ArgumentParser(description)
parser.add_argument("infile", type=str, help="Input tree")
parser.add_argument("type", choices=["AA", "NT"], help="Which type of alignment?")
parser.add_argument("outfile", type=argparse.FileType('w'), nargs="?", default=sys.stdout,
                        help="<output, optional> File to write output, defaults to stdout.")
args = parser.parse_args()

if args.type == "NT":
    cutoffs = [0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01,
               0.011, 0.012, 0.013, 0.014, 0.015, 0.016, 0.017, 0.018, 0.019,
               0.02, 0.021, 0.022, 0.023, 0.024, 0.025, 0.026, 0.027, 0.028, 0.029,
               0.03] # taken from chainsaw output with cutoff "None"
else:
    cutoffs = [0.0005, 0.0006, 0.0007, 0.0008, 0.0009, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01,
               0.011] # not yet adapted

writer = csv.writer(args.outfile)
writer.writerow(["cutoff", "nsubtrees", "mutual.inf", "normalized"])
for cutoff in cutoffs:
    phy = Phylo.read(args.infile, "newick")
    unroot(phy)
    subtrees = chainsaw(phy, cutoff)
    ## Debugging
    # ipdb.set_trace()
    minfo, norm_mi = mutual_info(subtrees, hema=(args.type == 'NT'))
    writer.writerow([cutoff, len(subtrees), minfo, norm_mi])
    args.outfile.flush()
