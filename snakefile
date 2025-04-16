# These scripts were cloned from the PoonLab FluClades repository.
# https://github.com/PoonLab/fluclades?tab=readme-ov-file
# The scripts were put in a workflow manager - Snakemake - and
# adapted for Enteroviruses

import os
from dotenv import load_dotenv, find_dotenv
load_dotenv(find_dotenv())

## Define paths and parameters
EMAIL = os.environ.get("EMAIL") # for Entrez transactions; taken from .env file
MIN_LENGTH = "6000" # min sequence length for whole genome run
ID_FIELD = "accession" # accession or strain, required for augur functions
EXCLUDE = "data/exclude.txt"
GFF_PATH = "data/genome_annotation.gff3" # annotation gff3 for EV-D68
REFERENCE_PATH = "data/reference.fasta" # EV-D68 reference sequence
ROOTING = "mid_point" # mid-point rooting of tree (using augur refine instead of midpoint.R)
INCL_PARAMS = False # include mindiv and maxpat in subtyping script -> if False, None is taken

## Run all rules in snakemake
rule all:
    input:
        "final_output.png",
        "chainsaw-ha.pdf",
        "gisaid-nseqs.pdf",
        "subtree-grid.pdf"

## get-metadata.py: replaced by ingest

rule index_sequences:
    message:
        """
        Creating an index of sequence composition for filtering
        """
    input:
        sequences = "data/sequences.fasta",
    output:
        sequence_index = "results/sequence_index.tsv"
    shell:
        """
        augur index \
            --sequences {input.sequences} \
            --output {output.sequence_index}
        """

rule filter:
    """
    Exclude sequences from before {MIN_DATE} and subsample to {MAX_SEQS} sequences.
    Only take sequences longer than {MIN_LENGTH}
    """
    input:
        sequences = "data/sequences.fasta",
        sequence_index = rules.index_sequences.output.sequence_index,
        metadata = "data/metadata.tsv",
        exclude = EXCLUDE,
    output:
        filtered_sequences = "results/filtered_sequences_raw.fasta",
        filtered_metadata = "results/filtered_metadata_raw.tsv",
    params: 
        min_length="" if MIN_LENGTH == "" else "--min-length " + MIN_LENGTH,
        strain_id_field = ID_FIELD,
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.sequence_index} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --exclude {input.exclude} \
            {params.min_length} \
            --output {output.filtered_sequences} \
            --output-metadata {output.filtered_metadata}
        """

rule relabel_fasta:
    """ 
    Relabel fasta sequences with metadata. Adding subgenotypes to fasta headers.
    """
    input:
        fasta = "results/filtered_sequences_raw.fasta",
        metadata = "results/filtered_metadata_raw.tsv",
        subgenotype = "data/RIVM_genotype.csv"
    output:
        "results/relabeled_seqs.fasta"
    shell:
        "python3 scripts/relabel-fasta.py {input.fasta} {input.metadata} {input.subgenotype} {output}"

# added subtype in relabel_fasta - using whole genome only (for now)
# rule filter_protein:
#     """ 
#     Filter protein sequences
#     """
#     input:
#         "sequences.fasta"
#     output:
#         "filtered_seqs.fasta"
#     shell:
#         "python scripts/filter-prot.py {input} {output}"


rule compress_sequences:
    """ 
    Compress sequences, remove duplicates
    """
    input:
        "results/relabeled_seqs.fasta"
    output:
        fasta = "results/compressed_seqs.fasta",
        csv = "results/duplicates.csv"
    shell:
        "python scripts/compress-seqs.py {input} {output.fasta} {output.csv}"


## concat-genes.py: skipped; specific to Influenza?

rule align:
    message:
        """
        Aligning sequences to {input.reference} using Nextclade3.
        """
    input:
        sequences = "results/compressed_seqs.fasta",
        reference = REFERENCE_PATH,
        annotation = GFF_PATH,
    output:
        alignment = "results/aligned.fasta",
        tsv = "results/nextclade.tsv",
    params:
        translation_template = lambda w: "results/translations/cds_{cds}.translation.fasta",
        #high-diversity 
        penalty_gap_extend = 1, #make longer gaps more costly - default is 0
        penalty_gap_open = 13,  #make gaps more expensive relative to mismatches - default is 13
        penalty_gap_open_in_frame = 18, #make gaps more expensive relative to mismatches - default is 7
        penalty_gap_open_out_of_frame = 23, #make out of frame gaps more expensive - default is 8 # prev was 19
        kmer_length = 6, #reduce to find more matches - default is 10
        kmer_distance = 25, #reduce to try more seeds - default is 50
        min_match_length = 30, #reduce to keep more seeds - default is 40
        allowed_mismatches = 15, #increase to keep more seeds - default is 8
        min_length = 30, # min_length - default is 100
        #cost of a mutation is 4
    shell:
        """
        nextclade3 run \
        -j {threads} \
        {input.sequences} \
        --input-ref {input.reference} \
        --input-annotation {input.annotation} \
        --penalty-gap-open {params.penalty_gap_open} \
        --penalty-gap-extend {params.penalty_gap_extend} \
        --penalty-gap-open-in-frame {params.penalty_gap_open_in_frame} \
        --penalty-gap-open-out-of-frame {params.penalty_gap_open_out_of_frame} \
        --kmer-length {params.kmer_length} \
        --kmer-distance {params.kmer_distance} \
        --min-match-length {params.min_match_length} \
        --allowed-mismatches {params.allowed_mismatches} \
        --min-length {params.min_length} \
        --include-reference false \
        --output-tsv {output.tsv} \
        --output-translations {params.translation_template} \
        --output-fasta {output.alignment} 
        """

rule tree:
    message:
        """
        Creating a maximum likelihood tree
        """
    input:
        alignment = "results/aligned.fasta",
    output:
        tree = "results/tree_raw.nwk",
    threads: 9
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --nthreads {threads}\
            --output {output.tree} \
        """

# midpoint.R: not needed; done by refine "--root mid_point".

rule refine:
    input:
        tree=rules.tree.output.tree,
        alignment="results/aligned.fasta",
    output:
        tree="results/tree.midpoint.nwk",
        node_data="results/branch_lengths.json",
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --root {ROOTING} \
            --keep-polytomies \
            --divergence-unit mutations-per-site \
            --output-node-data {output.node_data} \
            --output-tree {output.tree}
        """


rule subtyping: # output too poor
    """
    node-wise clustering of phylogenetic trees

    Calculate summary statistics for internal nodes of the tree, and 
    select subtrees rooted at those nodes on the basis of one or more
    criteria defined on these statistics.
    """
    input:
        "results/tree.midpoint.nwk"
    output:
        "results/mindiv0_01.maxpat1_2.subtypes.csv" if INCL_PARAMS else "results/subtree-grid.csv"
    params:
        yes = "--mindiv 0.01 --maxpat 1.2" if INCL_PARAMS else "" # tip-ot-tip distance cutoff = 1.2
    shell:
        """
        python scripts/subtyping.py {input} {output} {params.yes}
        """

CUTOFF = False # set to True to use cutoff
rule chainsaw:
    """
    Partition tree by cutting on internal branches with length 
    exceeding threshold.
    """
    input:
        infile="results/tree.midpoint.nwk"
    output:
        outfile="results/chainsaw_output.csv"
    params:
        cutoff="--cutoff 0.015 " if CUTOFF else "",  # or specify a float value
        nbin=20,
        format="summary" # "summary", "labels" or "trees"
    shell:
        """
        python scripts/chainsaw.py {input.infile} {output.outfile} \
        {params.cutoff} --nbin {params.nbin} --format {params.format}
        """

rule auto_chainsaw:
    """
    This script is used to generate the data required to produce Figures
    2A and 3A.  The inputs were trees reconstructed using FastTree2.
    Results are written to stdout in CSV format.
    """
    input:
        "results/tree.midpoint.nwk"
    output:
        "results/chainsaw-nsubtrees.csv"
    params:
        gene = "full" # "full" or "VP1",
    shell:
        """
        python scripts/auto-chainsaw.py {input}\
        {params.gene} \
        {output}
        """

################################
# !! R scripts are still in work !!
################################

rule plot_trees:
    """ 
    Create tree plots
    """
    input:
        "treeplots.RData"  #       tree = "results/tree_raw.nwk"
    output:
        "final_output.png"
    shell:
        "Rscript scripts/plot-trees.R"


rule chainsaw_plot:
    input:
        "results/chainsaw-nsubtrees.csv",
        "results/chainsaw-nsubtrees-na.csv",
        "results/chainsaw-HA-0.18.labels.csv",
        "results/chainsaw-NA-0.41.labels.csv",
        "results/chainsaw-nsubtrees-others.csv"
    output:
        "chainsaw-ha.pdf"
    shell:
        "Rscript scripts/chainsaw-plot.R"


rule coldates_plot:
    """ 
    Plot collection dates
    """
    input:
        "gisaid.csv"
    output:
        "gisaid-nseqs.pdf"
    shell:
        "Rscript scripts/coldates.R"


rule subtree_grid_plot:
    """ 
    Plot subtree grid
    """
    input:
        "results/subtree-grid.csv",
        "results/HA.mindiv0_08.maxpat1_2.subtypes.csv"
    output:
        "subtree-grid.pdf"
    shell:
        "Rscript scripts/subtree-grid.R"


