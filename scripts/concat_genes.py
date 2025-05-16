import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from BCBio import GFF
import pandas as pd
import ipdb
from collections import defaultdict

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Translate nucleotide alignment to amino acid polyprotein using CDS from GFF.")
    parser.add_argument("--translation", "-t", help="Path to CDS translation FASTA file")
    parser.add_argument("--annotation","-a", help="Path to GFF3 annotation file")
    parser.add_argument("--reference", "-ref", help="Path to reference FASTA file")
    parser.add_argument("-o","--output", help="Path to output amino acid FASTA")
    parser.add_argument("--cds", nargs="+", help="List of CDS names")
    parser.add_argument("--alignment", help="Path to nucleotide alignment file")
    args = parser.parse_args()
    # ipdb.set_trace()

    # Load the reference sequence (optional)
    seq_dict = SeqIO.to_dict(SeqIO.parse(args.reference, "fasta"))

    # Copy all translation fastas into a dataframe
    records = []
    for cds in args.cds:
        # print(f"Translating {cds}...")
        t = args.translation.replace("XX", cds,1)
        # print(t)
        for record in SeqIO.parse(t, "fasta"):
            records.append({"id": record.id, "seq": str(record.seq), "cds": cds})

    df = pd.DataFrame(records)
    # print(df.head())

    # Group by the id and concatenate the sequences
    grouped = df.groupby("id")["seq"].apply(lambda x: "".join(x)).reset_index()
    # print(grouped.head())

    # Write the output to a fasta file
    with open(args.output, "w") as out_f:
        for index, row in grouped.iterrows():
            out_f.write(f">{row['id']}\n")
            out_f.write(f"{row['seq']}\n")

        # Read the GFF3 file and extract the CDS sequences
    # Dictionary to store UTR lengths per gene
    utr_lengths = 0

    # Parse the GFF3 file
    with open(args.annotation) as in_handle:
        for rec in GFF.parse(in_handle, base_dict=seq_dict):
            for feature in rec.features:
                # ipdb.set_trace()
                if feature.type != "CDS":
                    utr_lengths += feature.location.end-feature.location.start

                if feature.type == "CDS":
                    utr_lengths -= feature.location.end-feature.location.start


    # Count how many IDs from the NT alignment exist in the AA translation
    matched_ids = 0
    length_mismatches = 0
    total_nt_seqs = 0

    # Make a dictionary for faster access to AA sequence lengths
    aa_lengths = {row["id"]: len(row["seq"]) for _, row in grouped.iterrows()}

    for record in SeqIO.parse(args.alignment, "fasta"):
        total_nt_seqs += 1
        nt_id = record.id
        nt_len = len(record.seq)

        if nt_id in aa_lengths:
            matched_ids += 1
            aa_len = aa_lengths[nt_id]

            # Each AA should come from 3 NTs (codon)
            if (nt_len-utr_lengths) // 3 != aa_len:
                length_mismatches += 1


    print(f"{length_mismatches} sequences had length mismatches between NT and AA.")
    print(f"{matched_ids} of {total_nt_seqs} NT sequences ({matched_ids / total_nt_seqs:.2%}) were found in the AA translation.")
