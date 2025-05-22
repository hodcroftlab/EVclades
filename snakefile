# These scripts were cloned from the PoonLab FluClades repository.
# https://github.com/PoonLab/fluclades?tab=readme-ov-file
# The scripts were put in a workflow manager - Snakemake - and
# adapted for Enteroviruses

import os
from dotenv import load_dotenv, find_dotenv
load_dotenv(find_dotenv())

configfile: "config.yaml"

## Define paths
EMAIL = os.environ.get("EMAIL") # for Entrez transactions; taken from .env file
ID_FIELD = "accession" # accession or strain, required for augur functions
EXCLUDE = "data/exclude.txt"
GFF_PATH = "data/genome_annotation.gff3" # annotation gff3 for EV-D68
REFERENCE_PATH = "data/reference.fasta" # EV-D68 reference sequence

## Define parameters
CDS_LIST = ["VP4", "VP2", "VP3", "VP1", "2A", "2B", "2C", "3A", "3B", "3C", "3D"]
MIN_LENGTH = "6000" # min sequence length for whole genome run
ROOTING = "mid_point" # mid-point rooting of tree (using augur refine instead of midpoint.R)
FORMAT = "labels" # "summary", "labels" or "trees"

# cut-off and max-date
MAX_DATE = "2017-01-01"
CUTOFF = 0.001 # set to None to get a list of different cutoff values; vary between 0.005 and 0.011
# CUTOFFs = [0.005, 0.006,0.007, 0.008,0.009,0.01,0.011, 0.012, 0.013, 0.014, 0.015]
# pdfunite results/plots/chainsaw-table_0.0*.pdf chainsaw-table_2018.pdf

#params for node-wise clustering
MINDIV = "0.01"
MAXPAT = "1.2"


wildcard_constraints:
    TYPE = "AA|NT" # amino acid or nucleotide run

## Run all rules in snakemake
rule all:
    input:
        expand("results/plots/chainsaw-table_{CUTOFF}_{TYPE}.pdf", TYPE=["AA", "NT"], CUTOFF=config["cutoff"]),
        expand("results/plots/treeplots_{TYPE}.png", TYPE=["AA", "NT"]),
        "results/plots/genbank-nseqs.pdf",
        alignment = "results/aligned_NT.fasta"
        # "results/plots/subtree-grid.pdf",
        # "results/plots/nsubtrees.pdf"

rule all_NT:
    input:
        expand("results/plots/chainsaw-table_{CUTOFF}_NT.pdf", 
               CUTOFF=[f"{config['cutoff']['NT']:.3f}"]),
        "results/plots/treeplots_NT.png",
        "results/plots/genbank-nseqs.pdf",
        "results/aligned_NT.fasta",
        "results/plots/subtree-grid_NT.pdf",
        "results/plots/nsubtrees_NT.pdf"

rule all_AA:
    input:
        expand("results/plots/chainsaw-table_{CUTOFF}_AA.pdf", 
               CUTOFF=[f"{config['cutoff']['AA']:.3f}"]),
        "results/plots/treeplots_AA.png",
        "results/aligned_NT.fasta"

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
        echo -e "\nrelabel_fasta:" 
        python3 scripts/relabel-fasta.py {input.fasta} {input.metadata} \
        {input.subgenotype} {output} 
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

rule compress_sequences:
    """ 
    Compress sequences, remove duplicates
    """
    input:
        seqs = lambda wildcards: "results/relabeled_seqs.fasta" if wildcards.TYPE == "NT" else "results/aligned_AA.fasta",
    params:
        type = lambda wildcards: "NT" if wildcards.TYPE == "NT" else "AA", 
    output:
        fasta = "results/compressed_seqs_{TYPE}.fasta",
        csv = "results/duplicates_{TYPE}.csv"
    shell:
        """
        echo -e "\ncompress_sequences:" 
        python scripts/compress-seqs.py {input.seqs} {output.fasta} {output.csv} {params.type} 
        """

rule align:
    message:
        """
        Aligning sequences to {input.reference} using Nextclade3.
        """
    input:
        sequences = "results/compressed_seqs_NT.fasta",
        reference = REFERENCE_PATH,
        annotation = GFF_PATH,
    output:
        alignment = "results/aligned_NT.fasta",
        tsv = "results/nextclade_NT.tsv",
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

rule concat_genes:
    message: "Translating polyprotein for all samples"
    params:
        translation = "results/translations/cds_XX.translation.fasta", # Translation output from Nextclade3; XX is the CDS name
        cds = " ".join(CDS_LIST),                               # List of CDS to translate

    input:
        nt_alignment = "results/aligned_NT.fasta",
        annotation = GFF_PATH,                               # GFF file with CDS annotations
        reference = REFERENCE_PATH,                         # Reference genome (optional for now)

    output:
        aa_alignment = "results/aligned_AA.fasta"  # Output AA polyprotein sequences
    shell:
        """
        python scripts/concat_genes.py \
            --translation {params.translation} \
            --cds {params.cds}\
            --annotation {input.annotation} \
            --reference {input.reference} \
            --alignment {input.nt_alignment} \
            -o {output.aa_alignment} 
        """

#TODO: More robust ML tree: -n 10, --epsilon 0.01, -B 1000 (add as params)
rule tree: 
    message:
        """
        Creating a maximum likelihood tree
        """
    input:
        alignment = lambda wildcards: "results/aligned_NT.fasta" if wildcards.TYPE == "NT" else "results/compressed_seqs_AA.fasta",
    params:
        # model = "auto" -> find optimal model; NT = "GTR+F+I+R4" with 2 cores (-ninit 2 -n 2 --epsilon 0.05); AA = FLAVI+F+R4 with 2 cores (-ninit 2 -n 2 --epsilon 0.05)
        model = lambda wildcards: "FLAVI+F+R4" if wildcards.TYPE == "AA" else "GTR+F+I+R4",  # "JTT" and "GTR" were defaults
    output:
        tree = "results/tree_raw_{TYPE}.nwk",
    threads: 2
    shell:
        """
        echo -e "\nCreating tree for {wildcards.TYPE} alignment..."
        augur tree \
            --alignment {input.alignment} \
            --method iqtree \
            --substitution-model {params.model} \
            --nthreads {threads} \
            --tree-builder-args "-ninit 2 -n 2 --epsilon 0.05" \
            --output {output.tree} \
            --override-default-args
        """

# midpoint.R: only needed for AA; NT: done by refine "--root mid_point".

rule refine:
    input:
        tree=rules.tree.output.tree,
        alignment="results/aligned_{TYPE}.fasta",
    output:
        tree="results/tree.midpoint_{TYPE}.nwk",
    params:
        node_data="results/branch_lengths_{TYPE}.json",

    shell:
        """
        if [ "{wildcards.TYPE}" == "AA" ]; then
            echo -e "\nRunning refine_aa for AA sequences..."
            Rscript scripts/midpoint.R {input.tree} {output.tree}

        else
            echo -e "\nRunning augur refine for NT sequences..."
            augur refine \
                --tree {input.tree} \
                --alignment {input.alignment} \
                --root {ROOTING} \
                --keep-polytomies \
                --divergence-unit mutations-per-site \
                --output-node-data {params.node_data} \
                --output-tree {output.tree}
        fi
        
        input_tips=$(grep -oE '[^\(\),]+' {input.tree} | sort | uniq | wc -l)
        output_tips=$(grep -oE '[^\(\),]+' {output.tree} | sort | uniq | wc -l)
        dropped=$((input_tips - output_tips))
        echo -e "\nrefine: Dropped sequences due to clock filter: $dropped" 
        """

        # python scripts/refine_aa.py \
        #         --tree {input.tree} \
        #         --alignment {input.alignment} \
        #         --divergence-unit mutations-per-site \
        #         --output-node-data {output.node_data} \
        #         --output-tree {output.tree} 

###### node-wise clustering ######
rule subtyping_params:
    message:
        """
        Subtyping (node-wise clustering) for {wildcards.TYPE} tree
        """
    input:
        "results/tree.midpoint_{TYPE}.nwk"
    output:
        f"results/mindiv{MINDIV}_maxpat{MAXPAT}.subtypes_{{TYPE}}.csv"
    params:
        mindiv = MINDIV,
        maxpat = MAXPAT,
    shell:
        """
        echo -e "\nnode-wise clustering" 
        python scripts/subtyping.py {input} {output} --mindiv {params.mindiv} --maxpat {params.maxpat}
        """


rule subtyping_grid:
    input:
        "results/tree.midpoint_{TYPE}.nwk"
    output:
        "results/subtree-grid_{TYPE}.csv"
    shell:
        """
        echo -e "\nsubtree-grid" 
        python scripts/subtyping.py {input} {output}
        """

rule subtree_grid_plot:
    """ 
    Plot subtree grid
    """
    input:
        grid = "results/subtree-grid_{TYPE}.csv",
        params_set = rules.subtyping_params.output
    output:
        out_pdf = "results/plots/subtree-grid_{TYPE}.pdf",
        nsubtrees = "results/plots/nsubtrees_{TYPE}.pdf"
    shell:
        "Rscript scripts/subtree-grid.R {input.grid} {input.params_set} {output.out_pdf} {output.nsubtrees} "

##### edge-wise clustering #####
# cutoff = None -> creates a list of different cutoff values
rule auto_chainsaw:
    """
    This script is used to generate the data required to produce Figures
    2A and 3A.  The inputs were trees reconstructed using FastTree2.
    Results are written to stdout in CSV format.
    """
    input:
        "results/tree.midpoint_{TYPE}.nwk"
    output:
        nsubtrees = "results/chainsaw-nsubtrees_{TYPE}.csv",
        values = "results/chainsaw_output_{TYPE}.csv"
    params:
        type = "{TYPE}", # "AA" or "NT",
        nbin = 20

    shell:
        """
        python scripts/chainsaw.py {input} {output.values} --nbin {params.nbin} >> {output.values} 2>&1
        python scripts/auto-chainsaw.py {input} {params.type} {output.nsubtrees}
        """

rule chainsaw:
    """
    Partition tree by cutting on internal branches with length 
    exceeding threshold.
    """
    input:
        infile = "results/tree.midpoint_{TYPE}.nwk"
    output:
        outfile = "results/chainsaw.subtrees_{TYPE}.nwk" if FORMAT == "trees" else f"results/chainsaw-{CUTOFF}.{FORMAT}_{{TYPE}}.csv"
    params:
        cutoff = lambda wildcards: "" if CUTOFF is None else f"--cutoff {CUTOFF}",
        nbin = 20
    shell:
        """
        echo -e "\nedge-wise clustering" 
        python scripts/chainsaw.py {input.infile} {output.outfile} \
        {params.cutoff} --nbin {params.nbin} --format {FORMAT}  
        """


##### Plotting #####
rule plot_trees:
    """ 
    Create tree plots
    """
    input:
        tree = "results/tree.midpoint_{TYPE}.nwk",
        # tree = "results/chainsaw.subtrees.nwk"
    output:
        plot_png = "results/plots/treeplots_{TYPE}.png",
        plot_pdf = "results/plots/inferred_{TYPE}.pdf"
    params:
        replace_labels = "True",
        cutoff_subtree = CUTOFF,
    shell:
        """
        Rscript scripts/plot-trees.R {input.tree} {output.plot_png} {output.plot_pdf} \
        {params.replace_labels} {params.cutoff_subtree} 
        """



rule chainsaw_plot:
    input:
        nsubtrees = "results/chainsaw-nsubtrees_{TYPE}.csv",
        label = f"results/chainsaw-{CUTOFF}.{FORMAT}_{{TYPE}}.csv"
    params:
        keep_NA = "False",
        cutoff_subtree = CUTOFF
    output:
        out_pdf = "results/plots/chainsaw_{TYPE}.tiff",
        out_table = f"results/plots/chainsaw-table_{CUTOFF}_{{TYPE}}.pdf"
    shell:
        """
        Rscript scripts/chainsaw-plot.R {input.nsubtrees} {input.label} \
        {output.out_pdf} {output.out_table} \
        {params.keep_NA} {params.cutoff_subtree}
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
        "Rscript scripts/coldates.R {input} {output} "

rule clean:
    shell:
        """
        find results/plots/ -mindepth 1 ! -type d -delete
        find results/translations/ -mindepth 1 ! -type d -delete
        find results/ -mindepth 1 ! -type d -delete
        rm log.out
        rm Rplots.pdf
        """