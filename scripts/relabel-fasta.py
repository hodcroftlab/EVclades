import argparse
from Bio import SeqIO
import pandas as pd
import ipdb

parser = argparse.ArgumentParser(
    description="Replace labels in CDS FASTA with metadata from CSV"
)
parser.add_argument("infile", type=argparse.FileType('r'),
                    help="<input FASTA> Sequences to relabel")
parser.add_argument("tsvfile", type=argparse.FileType('r'),
                    help="<input TSV> File with metadata")
parser.add_argument("subgenotype", type=argparse.FileType('r'),
                    help="<input CSV> File with subgenotypes")
parser.add_argument("outfile", type=argparse.FileType('w'),
                    help="<output FASTA> File to write relabeled sequences.")
args = parser.parse_args()

keys = ['accession', 'strain', 'host', 'country', 'date', 'subgenotype']
meta = pd.read_csv(args.tsvfile, sep="\t")

# Add subgenotypes to metadata
subgenotypes = pd.read_csv(args.subgenotype, names=['accession', 'subgenotype'], skiprows=1)
metadata = pd.merge(meta, subgenotypes, on='accession')[keys]

# remove all XX in date
metadata['date'] = metadata['date'].str.replace('-XX','')

# add subgenotypes as labels to sequences: "accession_subgenotype" 
records = SeqIO.parse(args.infile, "fasta")
for record in records:
    accn = record.description.split(".")[0].replace('lcl|', '')
    # ipdb.set_trace()
    if accn not in metadata['accession'].values:
        args.outfile.write(f">{accn}__\n{record.seq}\n")
    else:
        md = metadata[metadata['accession'] == accn]
        label = f"{md['accession'].values[0]}_{md['subgenotype'].values[0]}"
        args.outfile.write(f">{label}\n{record.seq}\n")