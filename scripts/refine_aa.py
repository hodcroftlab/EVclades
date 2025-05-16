import argparse
from Bio import AlignIO, Phylo
from Bio.Phylo.BaseTree import Tree
import numpy as np

def load_alignment_length(alignment_file: str) -> int:
    alignment = AlignIO.read(alignment_file, "fasta")
    return alignment.get_alignment_length()

def midpoint_reroot(tree: Tree):
    """
    Perform midpoint rerooting of a Bio.Phylo tree.
    """
    distances = {}
    terminals = tree.get_terminals()
    for i, t1 in enumerate(terminals):
        for t2 in terminals[i+1:]:
            dist = tree.distance(t1, t2)
            distances[(t1, t2)] = dist
    (t1, t2), _ = max(distances.items(), key=lambda x: x[1])
    tree.root_with_outgroup(t1, t2)
    tree.root_at_midpoint()
    print(f"\tMidpoint rerooted tree with outgroup {t1.name} and {t2.name}")
    return tree

def rescale_branch_lengths(tree: Tree, scale: float):
    for clade in tree.find_clades():
        if clade.branch_length is not None:
            clade.branch_length /= scale

def write_node_data(tree: Tree, output_file: str):
    with open(output_file, "w") as f:
        f.write("node_name\tparent_name\tbranch_length\n")
        for clade in tree.find_clades(order="level"):
            node_name = clade.name if clade.name else "internal"
            parent = tree.get_path(clade)
            parent_name = parent[-2].name if len(parent) >= 2 and parent[-2].name else "root"
            bl = clade.branch_length if clade.branch_length is not None else "NA"
            f.write(f"{node_name}\t{parent_name}\t{bl}\n")
    print(f"Node data written to {output_file}")

def process_tree(tree_file: str, alignment_file: str, divergence_unit: str, output_tree: str, output_node_data: str):
    print(f"\tProcessing tree file: {tree_file}")
    tree = Phylo.read(tree_file, "newick")
    tree = midpoint_reroot(tree)

    if divergence_unit == "mutations-per-site":
        if not alignment_file:
            raise ValueError("Alignment is required for --divergence-unit mutations-per-site")
        seq_len = load_alignment_length(alignment_file)
        print(f"\tAlignment length: {seq_len} sites")
        rescale_branch_lengths(tree, seq_len)

    Phylo.write(tree, output_tree, "newick")
    print(f"\tRerooted tree written to {output_tree}")

    if output_node_data:
        write_node_data(tree, output_node_data)

def main():
    parser = argparse.ArgumentParser(description="Reroot an amino acid tree using midpoint.")
    parser.add_argument("--tree", required=True, help="Input Newick tree file")
    parser.add_argument("--alignment", help="FASTA amino acid alignment")
    parser.add_argument("--divergence-unit", choices=["mutations", "mutations-per-site"], default="mutations",
                        help="Branch length unit (default: mutations)")
    parser.add_argument("--output-tree", required=True, help="Output Newick tree file")
    parser.add_argument("--output-node-data", help="Output node data file", default=None)

    args = parser.parse_args()

    process_tree(
        tree_file=args.tree,
        alignment_file=args.alignment,
        divergence_unit=args.divergence_unit,
        output_tree=args.output_tree,
        output_node_data=args.output_node_data,
    )

if __name__ == "__main__":
    main()