#!/usr/bin/env python3
import pandas as pd
import sys
import os

def summarize_lengths(lengths):
    mode_val = lengths.mode()
    mode_val = mode_val.iloc[0] if len(mode_val) > 0 else None

    return {
        "count": len(lengths),
        "min": lengths.min(),
        "q1": lengths.quantile(0.25),
        "median": lengths.median(),
        "mean": lengths.mean(),
        "mode": mode_val,
        "q3": lengths.quantile(0.75),
        "max": lengths.max(),
        "sd": lengths.std()
    }

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 gene_desert_stats.py gene_deserts/*.bed")
        sys.exit(1)

    files = sys.argv[1:]

    all_rows = []
    all_lengths = []

    for f in files:
        species = os.path.basename(f).replace(".gene_deserts.bed", "")

        df = pd.read_csv(f, sep="\t", header=None, names=["chr", "start", "end"])
        df["length"] = df["end"] - df["start"]

        stats_dict = summarize_lengths(df["length"])
        stats_dict["species"] = species
        all_rows.append(stats_dict)

        temp = pd.DataFrame({"species": species, "length": df["length"]})
        all_lengths.append(temp)

    # Save summary table
    summary_df = pd.DataFrame(all_rows)
    summary_df = summary_df[
        ["species", "count", "min", "q1", "median", "mean", "mode", "q3", "max", "sd"]
    ]
    summary_df.to_csv("desert_summary.tsv", sep="\t", index=False)

    # Save lengths for R
    lengths_df = pd.concat(all_lengths, ignore_index=True)
    lengths_df.to_csv("all_desert_lengths.tsv", sep="\t", index=False)

    print("✔ desert_summary.tsv generated")
    print("✔ all_desert_lengths.tsv generated")

if __name__ == "__main__":
    main()
