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
CUTOFF = 0.008 # set to None to get a list of different cutoff values; vary between 0.005 and 0.011
FORMAT = "labels" # "summary", "labels" or "trees"
MAX_DATE = "2017-01-01"
TYPE = "AA" # "AA" or "NT"

# CUTOFFs = [0.005, 0.006,0.007, 0.008,0.009,0.01,0.011, 0.012, 0.013, 0.014, 0.015]
# pdfunite results/plots/chainsaw-table_0.0*.pdf chainsaw-table_2018.pdf
    
## Run all rules in snakemake
rule all:
    input:
        f"results/plots/chainsaw-table_{CUTOFF}.pdf",
        "results/plots/treeplots.png",
        "results/plots/genbank-nseqs.pdf",
        # "results/plots/subtree-grid.pdf",
        # "results/plots/nsubtrees.pdf"


###### Alignment & Tree building ######
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
        max_date = "--max-date " + MAX_DATE,
    shell:
        """
        (echo -e filter:
        augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.sequence_index} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --exclude {input.exclude} \
            {params.min_length} \
            {params.max_date} \
            --output {output.filtered_sequences} \
            --output-metadata {output.filtered_metadata}) > log.out 2>&1
        """

rule relabel_fasta:
    """ 
    Relabel fasta sequences with metadata. 
    Adding RIVM subgenotypes to fasta headers to mimic flu strain names.
    """
    input:
        fasta = "results/filtered_sequences_raw.fasta",
        metadata = "results/filtered_metadata_raw.tsv",
        subgenotype = "data/RIVM_genotype_80.csv"
    output:
        "results/relabeled_seqs.fasta"
    shell:
        """
        echo -e "\nrelabel_fasta:" >> log.out
        python3 scripts/relabel-fasta.py {input.fasta} {input.metadata} \
        {input.subgenotype} {output} >> log.out 2>&1
        """

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

# run alignment & translation before compressing?

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
        """
        echo -e "\ncompress_sequences:" >> log.out
        python scripts/compress-seqs.py {input} {output.fasta} {output.csv} >> log.out 2>&1
        """


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
        --output-fasta {output.alignment} >> log.out 2>&1
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
        echo -e "\ntree:" >> log.out
        augur tree \
            --alignment {input.alignment} \
            --nthreads {threads}\
            --output {output.tree} >> log.out 2>&1
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
        
        input_tips=$(grep -oE '[^\(\),]+' {input.tree} | sort | uniq | wc -l)
        output_tips=$(grep -oE '[^\(\),]+' {output.tree} | sort | uniq | wc -l)
        dropped=$((input_tips - output_tips))
        echo -e "\nrefine: Dropped sequences due to clock filter: $dropped" >> log.out
        """

###### node-wise clustering ######
rule subtyping_params:
    input:
        "results/tree.midpoint.nwk"
    output:
        "results/mindiv0_01.maxpat1_2.subtypes.csv"
    shell:
        """
        echo -e "\nnode-wise clustering" >> log.out
        python scripts/subtyping.py {input} {output} --mindiv 0.01 --maxpat 1.2
        """

rule subtyping_grid:
    input:
        "results/tree.midpoint.nwk"
    output:
        "results/subtree-grid.csv"
    shell:
        """
        python scripts/subtyping.py {input} {output}
        """

rule subtree_grid_plot:
    """ 
    Plot subtree grid
    """
    input:
        grid = "results/subtree-grid.csv",
        params_set = "results/mindiv0_01.maxpat1_2.subtypes.csv"
    output:
        out_pdf = "results/plots/subtree-grid.pdf",
        nsubtrees = "results/plots/nsubtrees.pdf"
    shell:
        "Rscript scripts/subtree-grid.R {input.grid} {input.params_set} {output.out_pdf} {output.nsubtrees} >> log.out"

##### edge-wise clustering #####
rule chainsaw:
    """
    Partition tree by cutting on internal branches with length 
    exceeding threshold.
    """
    input:
        infile = "results/tree.midpoint.nwk"
    output:
        outfile = "results/chainsaw.subtrees.nwk" if FORMAT == "trees" else ("results/chainsaw_output.csv" if CUTOFF is None else f"results/chainsaw-{CUTOFF}.{FORMAT}.csv")  
    params:
        cutoff = lambda wildcards: "" if CUTOFF is None else f"--cutoff {CUTOFF}",
        nbin = 20
    shell:
        """
        echo -e "\nedge-wise clustering" >> log.out
        python scripts/chainsaw.py {input.infile} {output.outfile} \
        {params.cutoff} --nbin {params.nbin} --format {FORMAT} >> log.out 
        """

# cutoff = None -> creates a list of different cutoff values
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
        {params.gene} {output}
        """

##### Plotting #####
rule plot_trees:
    """ 
    Create tree plots
    """
    input:
        tree = "results/tree.midpoint.nwk",
        # tree = "results/chainsaw.subtrees.nwk"
    output:
        plot_png = "results/plots/treeplots.png",
        plot_pdf = "results/plots/inferred.pdf"
    params:
        replace_labels = "False",
        cutoff_subtree = {CUTOFF},
    shell:
        """
        Rscript scripts/plot-trees.R {input.tree} {output.plot_png} {output.plot_pdf} \
        {params.replace_labels} {params.cutoff_subtree} >> log.out
        """



rule chainsaw_plot:
    input:
        nsubtrees = "results/chainsaw-nsubtrees.csv",
        label = f"results/chainsaw-{CUTOFF}.{FORMAT}.csv",
        # label = "results/chainsaw-0.008.labels.csv",
    params:
        keep_NA = "False", # "True" or "False"
        cutoff_subtree = {CUTOFF}, # best number of subtrees: 0.008
    output:
        out_pdf = "results/plots/chainsaw.tiff",
        out_table = f"results/plots/chainsaw-table_{CUTOFF}.pdf" 
    shell:
        """
        Rscript scripts/chainsaw-plot.R {input.nsubtrees} {input.label} \
        {output.out_pdf} {output.out_table} \
        {params.keep_NA} {params.cutoff_subtree} \
        >> log.out
        """


rule coldates_plot: # Figure S1
    """ 
    Plot collection dates
    """
    input:
        "data/metadata.tsv"
    output:
        "results/plots/genbank-nseqs.pdf"
    shell:
        "Rscript scripts/coldates.R {input} {output} >> log.out"

rule clean:
    shell:
        """
        find results/plots/ -mindepth 1 ! -type d -delete
        find results/translations/ -mindepth 1 ! -type d -delete
        find results/ -mindepth 1 ! -type d -delete
        rm log.out
        rm Rplots.pdf
        """