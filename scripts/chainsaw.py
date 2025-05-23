description = """
Partition tree by cutting on internal branches with length 
exceeding threshold.
"""

from Bio import Phylo
from Bio.Phylo.BaseTree import Tree
import argparse
import bisect
import sys
import re
from math import log


def get_parents(phy):
    # Clades do not store references to parents
    parents = {}
    for parent in phy.get_nonterminals():
        for child in parent.clades:
            parents.update({child: parent})
    return parents


def unroot(phy):
    """ 
    Shorthand, designates a trifurcating Clade as root.
    Note this sets phy.rooted=True
    """
    if len(phy.root.clades) == 2:
        phy.root_with_outgroup(phy.root.clades[0])


def cuttree(phy, clade):
    """
    Cut a tree into two subtrees at the branch associated with a clade
    :param phy: BaseTree.Tree object
    :param clade: BaseTree Clade object, corresponds to internal branch to cut
    :returns: new BaseTree.Tree object rooted on clade; input `phy` is modified 
              in-place with clade removed.
    """
    if clade.is_terminal():
        sys.stderr.write("Warning: cannot cut tree on terminal branch. Use prune method. "
                         "Returning unaltered tree.\n")
        return phy
    
    parents = get_parents(phy)
    parent = parents.get(clade, None)
    if parent is None:
        sys.stderr.write("Clade has no parent in tree\n")
        return phy
    
    clade.branch_length = 0
    subtree1 = Tree(clade)
    unroot(subtree1)

    parent.clades.remove(clade)
    nchild = len(parent.clades)
    
    grandpar = parents.get(parent, None)
    if grandpar is None:
        # parent is current root    
        if nchild == 0:
            sys.stderr.write("cuttree: parent had only one child\n")
            sys.exit()
        elif nchild == 1:
            # parent is stem of one subtree, drop it
            subtree2 = Tree(parent.clades[0])
            del parent
        else:
            subtree2 = Tree(parent)
    else:
        if nchild < 2:
            for child in parent.clades:
                child.branch_length += parent.branch_length
                grandpar.clades.append(child)
            grandpar.clades.remove(parent)
            del parent
        phy.root_with_outgroup(grandpar)
        subtree2 = Tree(phy.root)
    
    subtree2.root.branch_length = 0
    unroot(subtree2)
    return subtree1, subtree2


def longest(phy):
    """ Get longest internal branch of input tree """
    nodes = phy.get_nonterminals()
    if len(nodes) < 2:
        return None  # only terminal descendants left
    intermed = [(node.branch_length, ni) for ni, node in enumerate(nodes)
                if node.branch_length]
    intermed.sort(reverse=True)
    return nodes[intermed[0][1]]


def chainsaw(phy, cutoff):
    """
    Locate the longest internal branch and cut if its length exceeds threshold
    :param phy:  Bio.Phylo BaseTree object
    :param cutoff:  float, branch length threshold
    :return:  list, BaseTree objects produced from partition into subtrees
    """
    # initialize list with original tree
    subtrees = [phy]
    cutting = True  # enter loop
    while cutting:
        cutting = False
        newtrees = []
        # for each subtree in list
        for subtree in subtrees:
            node = longest(subtree)  # locate longest internal branch
            if node and node.branch_length > cutoff:
                # if branch length exceeds threshold, cut the branch
                cutting = True
                # replace the subtree with the two new subtrees
                st1, st2 = cuttree(subtree, node)
                newtrees.extend([st1, st2])
            else:
                newtrees.append(subtree)
        subtrees = newtrees
        # continue until no internal branches are longer than threshold
    return subtrees


def mutual_info(subtrees, hema=True):
    """
    Calculate the normalized mutual information of HnNn labels and
    subtrees as different partitions of the data.
    :param subtrees:  list, BaseTree objects
    :return:  (float, float), mutual information (normalized)
    """
    # entabulate serotype labels per subtree
    tab = []
    if hema:
            pat = re.compile(r"_((?:[A-Z]\d+|[A-Z]{1,3}))_")  # match A1, A2, B1, B2, B3, US00, ... EV-D68 subgenotypes
    else:
        pat = re.compile(".+_H[0-9]+(N[0-9]+)_.+")  # match Nn

    total = 0  # number of tips with any HnNn label
    stcount = {}  # by subtree
    serocount = {}  # total HnNn counts

    for idx, subtree in enumerate(subtrees):
        stcount.update({idx: 0})
        tips = [tip.name for tip in subtree.get_terminals()]
        labels = {}  # HnNn counts in this subtree
        for tip in tips:
            matches = pat.findall(tip)
            if matches:
                label = matches[0]
                if label not in labels:
                    labels.update({label: 0})
                labels[label] += 1
                total += 1
                stcount[idx] += 1
                if label not in serocount:
                    serocount.update({label: 0})
                serocount[label] += 1
        tab.append(labels)

    # calculate mutual information
    minfo = 0
    for i, subtree in enumerate(tab):
        pi = stcount[i] / total  # marginal freq of subtree
        for serotype, count in subtree.items():
            pj =  serocount[serotype] / total  # marginal freq of serotype
            pij = count / total  # cell-specific joint prob
            if pij == 0:
                continue
            minfo += pij * log(pij / (pi*pj))

    # normalization
    # hu = -1.0 * sum([count / total * log(count / total) for _, count in stcount.items()])
    # hu calculation failed in auto-chainsaw -> added a count > 0 check
    hu = -1.0 * sum([(count / total) * log(count / total) for _, count in stcount.items() if count > 0])
    hv = -1.0 * sum([count / total * log(count / total) for _, count in serocount.items()])
    norm_minfo = 2 * minfo / (hu + hv)
    return minfo, norm_minfo


if __name__ == "__main__":
    # command line interface
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("infile", type=str,
                        help="<input> Path to Newick file to process.")
    parser.add_argument("outfile", type=argparse.FileType('w'), nargs="?", default=sys.stdout,
                        help="<output, optional> File to write output, defaults to stdout.")
    parser.add_argument("--cutoff", type=float,
                        help="<input> Maximum internal branch length. "
                        "If none given, display summary of lengths.")
    parser.add_argument("--nbin", type=int, default=20,
                        help="<option> number of bins for branch length summary.")
    parser.add_argument("-f", "--format", choices=["summary", "labels", "trees"],
                        default="summary",
                        help="<option> Format to write output. Defaults to 'summary'. "
                        "'labels' writes all tip labels for each subtree index. "
                        "'trees' writes a set of Newick tree strings.")
    args = parser.parse_args()

    phy = Phylo.read(args.infile, 'newick')
    unroot(phy)

    if args.cutoff is None:
        # display some summary statistics of input tree and quit
        bl = [node.branch_length for node in phy.get_nonterminals() if node is not phy.root]
        bl.sort()  # ascending order
        print(bl[-10:])
        # display histogram
        blmax = bl[-1]
        blstep = blmax/args.nbin
        cutval = 0
        left = 0
        for i in range(args.nbin):
            cutval += blstep
            right = bisect.bisect_left(bl, cutval)
            print(cutval, right-left)
            left = right
        sys.exit()

    # main routine
    subtrees = chainsaw(phy, args.cutoff)

    # write output
    if args.format == "labels":
        args.outfile.write("subtree,tip.label\n")
        for idx, subtree in enumerate(subtrees):
            for tip in subtree.get_terminals():
                args.outfile.write(f"{idx},\"{tip.name}\"\n")
    elif args.format == "summary":
        args.outfile.write("subtree,ntips,tip.label\n")
        for idx, subtree in enumerate(subtrees):
            tips = subtree.get_terminals()
            args.outfile.write(f"{idx},{len(tips)},{tips[0].name}\n")
    else:
        for subtree in subtrees:
            Phylo.write(subtree, args.outfile, "newick")

    if args.outfile is not sys.stdout:
        args.outfile.close()
