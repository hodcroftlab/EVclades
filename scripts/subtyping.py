"""
node-wise clustering of phylogenetic trees

Calculate summary statistics for internal nodes of the tree, and 
select subtrees rooted at those nodes on the basis of one or more
criteria defined on these statistics.

Probably a much more effective method will be to "cut" up the tree 
on internal branches with lengths below some threshold.  A nice 
feature of this method is that the input tree doesn't have to be rooted.
"""

from Bio import Phylo
import sys
import random
import csv
import argparse
import re
import ipdb

parser = argparse.ArgumentParser(
    description="Compute subtree-level statistics throughout tree and extract clusters."
)
parser.add_argument("infile", type=argparse.FileType('r'),
                    help="<input NWK> Newick string for rooted tree.")
parser.add_argument("outfile", type=argparse.FileType('w'),
                    help="<output CSV> File to write subtyping results.")
parser.add_argument("--mindiv", type=float, default=None, required=False,
                    help="<optional> Minimum divergence between subtrees.")
parser.add_argument("--maxpat", type=float, default=None, required=False,
                    help="<optional> Maximum mean patristic distance within subtree.")
args = parser.parse_args()

phy = Phylo.read(args.infile, 'newick')

def choose2(x):
    return x*(x-1)/2

# calculate mean node-to-tip distance from tips down to root
print("Processing tips...")
for tip in phy.get_terminals():
    tip.ntips = 1
    tip.bl_tot = 0
    tip.path_tot = 0  
    tip.mean_pat = 0

# store parent links and "V"-lengths
for parent in phy.get_nonterminals():
    vlen = sum([(x.branch_length or 0) for x in parent.clades])
    for child in parent.clades:
        child.parent = parent
        child.vlen = vlen
        # TODO: this doesn't really work, figure out why....

print("Calculating subtree stats...")
for node in phy.get_nonterminals(order="postorder"):
    node.ntips = 0
    node.bl_tot = 0  # sum all branch lengths above this node
    node.path_tot = 0  # sum node-to-tip lengths
    
    numer = 0   # terms to calculate mean patristic distance
    denom = 0
    nnprod = 1
    dlsum = 0
    
    for child in node.clades:
        node.ntips += child.ntips
        node.bl_tot += child.branch_length + child.bl_tot
        node.path_tot += child.branch_length * child.ntips + child.path_tot
        numer += choose2(child.ntips)*child.mean_pat
        denom += choose2(child.ntips)
        nnprod *= child.ntips
        dlsum += child.mean_pat + child.branch_length
    
    node.mean_pat =  (numer + nnprod*dlsum) / (denom + nnprod)

# ipdb.set_trace()
# sanity check that path_tot is correct
#test_set = [node for node in phy.get_nonterminals() if node.ntips > 10 and node.ntips < 100]
#for node in phy.test_set:
#    node.correct = sum([node.distance(tip) for tip in node.get_terminals()])

# calculate means
for node in phy.get_nonterminals():
    node.bl_mean = node.bl_tot / node.ntips
    node.path_mean = node.path_tot / node.ntips


# def find_subtrees(minlen, maxlen, child=None):
#     """ recursive function, assumes outer scope variable subtrees """
#     if child is None:
#         child = phy.root
#     if (child.branch_length is not None and child.vlen > minlen and child.mean_pat < maxlen):  # child.path_mean
#         subtrees.append(child)
#     else:
#         for child in child.clades:
#             if child.is_terminal():
#                 continue
#             find_subtrees(minlen, maxlen, child)

def find_subtrees(minlen, maxlen, child=None):
    if child is None:
        child = phy.root

    # make sure the attributes exist
    if hasattr(child, 'vlen') and hasattr(child, 'mean_pat'):
        # print(f"Checking node: vlen={child.vlen:.4f}, mean_pat={child.mean_pat:.4f}")

        if (child.branch_length is not None and 
            child.vlen > minlen and 
            child.mean_pat < maxlen):
            subtrees.append(child)
            return  # do not go deeper once a subtree is collected

    if not child.is_terminal():
        for subchild in child.clades:
            find_subtrees(minlen, maxlen, subchild)


# def find_subtrees(minlen, maxlen, child=None):
#     """ recursive function, assumes outer scope variable subtrees """
#     if child is None:
#         child = phy.root
#     # print(f"child.branch_length: {child.branch_length}")
#     # print(f"child.mean_pat: {child.mean_pat}")
#     # print(f"minlen: {minlen}")
#     # print(f"maxlen: {maxlen}")
#     # print(f"child.ntips: {child.ntips}")
#     # ipdb.set_trace()

#     for sub_child in child.clades:
#         # print(f"sub_child.branch_length: {sub_child.branch_length}")
#         # print(f"sub_child.vlen: {sub_child.vlen}")
#         # print(f"sub_child.mean_pat: {sub_child.mean_pat}")
#         if (sub_child.branch_length is not None and sub_child.vlen > minlen and sub_child.mean_pat < maxlen):
#             subtrees.append(sub_child)
#         else:
#             if not sub_child.is_terminal():
#                 find_subtrees(minlen, maxlen, sub_child)


def count_labels(subtree): ## H3N8 for example; subtypes in flu
    pat = re.compile(r"_([A-Z0-9]{1,5})_") #"_(H[0-9]+)(N[0-9]+)?_")
    ## now it will only recognize strings if they have one letter and one digit
    counts = {}
    tips = subtree.get_terminals()
    for tip in tips:
        match = pat.findall(tip.name)
        if len(match) ==0: 
            continue
        serotype = match[0]
        if serotype not in counts:
            counts.update({serotype: 0})
        counts[serotype] += 1
    # print(counts) ## debuggin
    return counts


print("Searching for subtrees...")

if args.mindiv is None or args.maxpat is None:
    # generate subtrees over a grid of threshold settings
    writer = csv.writer(args.outfile)
    writer.writerow(["minlen", "maxlen", "nsubtrees", "nlabels", "mean.n.types"])
    
    tot_labels = sum(count_labels(phy.root).values())
    # ipdb.set_trace()
    # for minl in [i/100 for i in range(10)]:
    #     # varying pendant branch length
    #     for maxl in [j/10 for j in range(1, 21)]:
    #         # varying mean patristic distance
    for minl in [i / 1000 for i in range(0, 45, 5)]:  # 0.000 to 0.040 by 0.005
        for maxl in [i / 1000 for i in range(10, 61, 5)]:  # 0.010 to 0.06 by 0.005
            subtrees = []
            find_subtrees(minl, maxl)
            nsubtrees = len(subtrees)
            nlabels = 0
            ntypes = 0
            for subtree in subtrees:
                counts = count_labels(subtree)  
                ntypes += len(counts)  # mean number of label types per subtree
                nlabels += sum(counts.values())  # total number of labels in subtrees

            if nsubtrees == 0:
                print(f"No subtrees for minl={minl}, maxl={maxl}") # used for debugging and fixing minl and maxl
                mean_types = 0
            else:
                mean_types = ntypes / nsubtrees

            writer.writerow([minl, maxl, nsubtrees, nlabels/tot_labels, mean_types])

            # writer.writerow([minl, maxl, nsubtrees, nlabels/tot_labels, 
            #                 ntypes / nsubtrees])

    args.outfile.close()
else:
    # write out detailed output for one set of thresholds
    subtrees = []
    find_subtrees(args.mindiv, args.maxpat)

    print(f"Found {len(subtrees)} subtrees, writing outputs...")
    writer = csv.writer(args.outfile)
    writer.writerow(['subtree', 'serotype', 'count'])

    for i, subtree in enumerate(subtrees):
        counts = count_labels(subtree)    
        for serotype, count in counts.items():
            writer.writerow([i, serotype, count])

    args.outfile.close()
