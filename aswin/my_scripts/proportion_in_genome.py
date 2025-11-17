import pandas as pd
import argparse

def read_genome_sizes(file):
    df = pd.read_csv(file, sep="\t", header=None, names=["chr", "size"])
    return df

def read_bed(file):
    df = pd.read_csv(file, sep="\t", header=None, names=["chr", "start", "end"])
    return df

def merge_intervals(df):
    merged = []
    for chrom, group in df.groupby("chr"):
        group = group.sort_values("start")
        cur_start, cur_end = None, None
        
        for _, row in group.iterrows():
            if cur_start is None:
                cur_start, cur_end = row["start"], row["end"]
            elif row["start"] <= cur_end:
                cur_end = max(cur_end, row["end"])
            else:
                merged.append((chrom, cur_start, cur_end))
                cur_start, cur_end = row["start"], row["end"]
        
        if cur_start is not None:
            merged.append((chrom, cur_start, cur_end))
    
    return pd.DataFrame(merged, columns=["chr", "start", "end"])


def calculate_coverage(genome_df, desert_df):
    result_rows = []
    total_genome = genome_df["size"].sum()

    # merge overlapping gene desert intervals
    merged_deserts = merge_intervals(desert_df)

    total_desert_bp = 0

    for _, row in genome_df.iterrows():
        chrom = row["chr"]
        chrom_size = row["size"]

        deserts = merged_deserts[merged_deserts["chr"] == chrom]

        chrom_desert_bp = (deserts["end"] - deserts["start"]).sum()
        total_desert_bp += chrom_desert_bp

        result_rows.append([
            chrom,
            chrom_size,
            chrom_desert_bp,
            chrom_desert_bp / chrom_size
        ])

    genome_covered = total_desert_bp / total_genome
    genome_not_covered = 1 - genome_covered

    summary = {
        "total_genome_bp": total_genome,
        "total_desert_bp": total_desert_bp,
        "percent_covered": genome_covered,
        "percent_not_covered": genome_not_covered
    }

    return pd.DataFrame(result_rows, columns=[
        "chr", "chrom_size", "desert_bp", "fraction_covered"
    ]), summary


def main(genome_file, desert_file):
    genome_df = read_genome_sizes(genome_file)
    desert_df = read_bed(desert_file)

    per_chr_df, summary = calculate_coverage(genome_df, desert_df)

    print("\n=== Per Chromosome Coverage ===")
    print(per_chr_df.to_string(index=False))

    print("\n=== Genome-wide Summary ===")
    print(f"Total genome size: {summary['total_genome_bp']:,} bp")
    print(f"Total desert bases: {summary['total_desert_bp']:,} bp")
    print(f"Genome covered by deserts: {summary['percent_covered']*100:.4f}%")
    print(f"Genome NOT covered: {summary['percent_not_covered']*100:.4f}%")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute genome coverage by gene deserts.")
    parser.add_argument("genome_file", help="Input genome sizes file")
    parser.add_argument("desert_file", help="Input gene desert BED file")
    args = parser.parse_args()

    main(args.genome_file, args.desert_file)
