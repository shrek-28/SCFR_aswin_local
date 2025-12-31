from pathlib import Path
import polars as pl 
from collections import defaultdict 
import glob
import re

## 1. DIRECTORY SETUP
def setup_species_workspace(base_dir, species:str) -> Path:
    '''
    creates a workspace directory inside exon_shadow/ 
    announces that it is being processed
    setup step
    '''
    out_dir = base_dir/"exon_shadow"/species
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"workspace created: {out_dir}")
    return out_dir

## 2: CONVERSION OF SCFR .out TO ZERO-BASED BED CSV
def scfr_to_bed_csv(scfr_file: Path, output_file:Path):
    ## reads the tab separated scfr as a csv file
    df = pl.read_csv(scfr_file, sep="\t", has_header=False) 

    ## filters where column 3 is not null 
    df = df.filter(pl.col(3).is_not_null())

    ## strand
    df = df.with_columns([
        pl.lit("1").alias("score"),
        pl.when(pl.col(3) == "-").then("-").otherwise("+").alias("strand")
    ])

    ## makes csv
    df.write_csv(output_file)

    ## saves csv and returns the new dataframe 
    print(f"SCFR BED CSV saved: {output_file}")
    return df

## 3: PARSE GTF ATTRIBUTES (WAS PRESENT IN METADATA)
def parse_attributes(attr_str: str) -> dict:
    attrs = {}
    for field in attr_str.strip().strip(";"):
        field = field.strip()
        if not field:
            continue
        if "=" in field:
            k, v = field.split("=", 1)
        else:
            parts = field.split(None, 1)
            if len(parts)!=2:
                continue
            k, v = parts
        attrs[k] = v.strip().strip('"')
    return attrs

## 4: getting gene name (present in metadata file)
def get_gene_name(attrs):
    for k in ["gene_name", "gene", "Name", "gene_id", "ID", "Parent"]:
        if k in attrs:
            return attrs[k].split(",")[0]
    return "NA"

# 5: load genome sizes (present in metadata file)
def get_chrom_sizes(species_pattern: str, base_dir: Path) -> dict:
    chrom_sizes = {}
    base_path = base_dir / "genome_reports" / "*.tsv"
    report_files = glob.glob(str(base_path))
    target_files = [f for f in report_files if species_pattern.lower() in f.lower()]
    if not target_files:
        print(f"Warning: No genome report found for {species_pattern}")
        return chrom_sizes
    target_file = target_files[0]
    with open(target_file, 'r') as f:
        for line in f:
            if line.startswith("#") or "Seq length" in line:
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 11:
                chrom = parts[8]
                try:
                    size = int(parts[10])
                    chrom_sizes[chrom] = size
                except ValueError:
                    continue
    return chrom_sizes    

# 7: FRAME CALCULATION (present in metadata file)
def calculate_genomic_frame(start: int, end: int, strand: str, phase: str, chrom_size: int) -> str:
    try:
        if phase == ".":
            return "NA"
        p = int(phase)
        if strand == "+":
            frame = ((start - (start - (start % 3)) + p - 1) % 3) + 1
            return str(frame)
        elif strand == "-":
            val = (chrom_size - end + 1)
            frame = abs((val % 3) + p)
            return str(frame if frame != 0 else 3)
    except:
        return "NA"
    return "NA"

# 8: PARSE GTF AND GET EXON METADATA (present in metadata file)
def parse_gtf_cds(gtf_file: Path, chrom_sizes: dict):
    cds_by_tx = defaultdict(list)
    gene_to_txs = defaultdict(set)
    exon_to_txs = defaultdict(set)

    with open(gtf_file, 'r') as fh:
        for line in fh:
            if line.startswith("#"): 
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9: 
                continue
            chrom, _, feature, start, end, _, strand, phase, attr_str = fields
            if feature != "CDS":
                continue

            start, end = int(start), int(end)
            chrom_size = chrom_sizes.get(chrom, 0)
            frame = calculate_genomic_frame(start, end, strand, phase, chrom_size)

            attrs = parse_attributes(attr_str)
            gene = get_gene_name(attrs)
            tx = attrs.get("transcript_id", "NA")

            exon_key = (chrom, start, end, gene, frame)
            cds_by_tx[tx].append({
                "chrom": chrom, "start": start, "end": end,
                "gene": gene, "tx": tx, "frame": frame, "strand": strand
            })
            gene_to_txs[gene].add(tx)
            exon_to_txs[exon_key].add(tx)

    return cds_by_tx, gene_to_txs, exon_to_txs

# 9: WRITE BED-LIKE CSV WITH EXON METADATA
def write_cds_bed_csv(output_file: Path, cds_by_tx, gene_to_txs, exon_to_txs):
    rows = []
    for tx, exons in cds_by_tx.items():
        strand = exons[0]["strand"]

        ## exon sorting in 5'->3' order
        if strand == "+":
            exons.sort(key=lambda x: x["start"])
        else:
            exons.sort(key=lambda x: x["start"], reverse=True)

        total_exons = len(exons)
        for exon_number, e in enumerate(exons, start=1):
            # getting exon gene
            gene = e["gene"]

            # getting chromosome number, start, end and frame
            exon_key = (e["chrom"], e["start"], e["end"], gene, e["frame"])

            ## from metadata metrics
            gene_tx_count = len(gene_to_txs[gene])
            exon_tx_count = len(exon_to_txs[exon_key])
            sharing = exon_tx_count / gene_tx_count

            ## EXON ORDER AND SPLICING
            exon_order = "first" if exon_number == 1 else ("last" if exon_number == total_exons else "middle")
            splicing = "constitutive" if sharing == 1.0 else ("unique" if exon_tx_count == 1 else "alternative")

            ## APPENDING RESULT
            rows.append({
                "chrom": e["chrom"], "start": e["start"], "end": e["end"],
                "gene": gene, "tx": e["tx"], "strand": e["strand"], "frame": e["frame"],
                "exon_number": exon_number, "gene_tx_count": gene_tx_count,
                "exon_tx_count": exon_tx_count, "exon_order": exon_order,
                "sharing": round(sharing, 3), "splicing": splicing
            })

    df = pl.DataFrame(rows)
    df.write_csv(output_file)
    print(f"CDS BED CSV saved: {output_file}")

## PIPELINE 
def gtf_to_bed_pipeline(gtf_file: Path, output_file: Path, species: str, base_dir: Path):
    chrom_sizes = get_chrom_sizes(species, base_dir)
    cds_by_tx, gene_to_txs, exon_to_txs = parse_gtf_cds(gtf_file, chrom_sizes)
    write_cds_bed_csv(output_file, cds_by_tx, gene_to_txs, exon_to_txs)

## overlap computation
def compute_scfr_cds_overlaps(scfr_csv: Path, cds_csv: Path, output_csv: Path):
    # Load SCFR and CDS CSVs
    scfr_df = pl.read_csv(scfr_csv)
    cds_df = pl.read_csv(cds_csv)

    # Make sure columns exist
    for col in ["start", "end", "strand", "frame"]:
        if col not in scfr_df.columns:
            raise ValueError(f"{col} missing in SCFR CSV")
        if col not in cds_df.columns:
            raise ValueError(f"{col} missing in CDS CSV")

    # Prepare for overlaps: merge all combinations of SCFR × CDS with same strand
    overlaps = []

    # Iterate over SCFR rows
    for scfr_row in scfr_df.iter_rows(named=True):
        strand = scfr_row["strand"]
        scfr_start, scfr_end = scfr_row["start"], scfr_row["end"]

        # Filter CDS by same strand
        cds_candidates = cds_df.filter(pl.col("strand") == strand)

        # Overlap condition: scfr.start <= cds.end and scfr.end >= cds.start
        cds_overlap = cds_candidates.filter(
            (pl.col("start") <= scfr_end) & (pl.col("end") >= scfr_start)
        )

        for cds_row in cds_overlap.iter_rows(named=True):
            combined = {**scfr_row, **cds_row}  # merge dictionaries
            overlaps.append(combined)

    overlap_df = pl.DataFrame(overlaps)
    overlap_df.write_csv(output_csv)
    print(f"SCFR × CDS overlaps saved: {output_csv}")
    return overlap_df

## Filter SCFRs containing exons with same frame and calculate shadow
def filter_containing_cds_with_shadow(overlap_csv: Path, output_csv: Path):
    df = pl.read_csv(overlap_csv)

    shadow_rows = []

    for row in df.iter_rows(named=True):
        # SCFR frame (remove any "-")
        scfr_frame = str(row["frame"]).replace("-", "")
        cds_frame = str(row.get("frame_right", row.get("frame")))  # depends on merge columns

        # Overlap coordinates
        scfr_start, scfr_end = row["start"], row["end"]
        cds_start, cds_end = row["start_right"] if "start_right" in row else row["start"], row["end_right"] if "end_right" in row else row["end"]

        # Check if SCFR contains CDS exon and frame matches
        if scfr_frame == cds_frame and scfr_start <= cds_start and scfr_end >= cds_end:
            strand = row["strand"]
            if strand == "+":
                upstream = cds_start - scfr_start
                downstream = scfr_end - cds_end
            else:
                upstream = scfr_end - cds_end
                downstream = cds_start - scfr_start

            row_dict = dict(row)
            row_dict["upstream_shadow"] = upstream
            row_dict["downstream_shadow"] = downstream
            shadow_rows.append(row_dict)

    shadow_df = pl.DataFrame(shadow_rows)
    shadow_df.write_csv(output_csv)
    print(f"SCFR containing CDS with shadow saved: {output_csv}")
    return shadow_df

## OVERLAP OUTPUTS
def initialize_overlap_outputs(out_dir: Path, species: str):
    single_exon_file = out_dir / f"{species}_single_exon.csv"
    multi_exon_file = out_dir / f"{species}_multi_exon.csv"
    exitron_file = out_dir / f"{species}_exitron_candidates.csv"

    single_cols = [
        "chr", "start", "end", "frame", "filler", "strand",
        "chr_exon", "exon_start", "exon_end", "gene", "transcript",
        "exon_strand", "exon_frame", "exon_number", "gene_tx_count",
        "exon_tx_count", "exon_order", "exon_sharing", "exon_splicing",
        "overlap_len", "upstream_len_in_scfr", "downstream_len_in_scfr", "merged_txs"
    ]
    multi_cols = [
        "chr", "start", "end", "frame", "filler",
        "strand", "first_exon_start", "last_exon_end", "gene", "transcript",
        "exon_strand", "exon_frame", "exon_count", "upstream_len_in_scfr",
        "downstream_len_in_scfr", "first_exon_number", "last_exon_number",
        "gene_tx_count", "first_exon_tx_count", "last_exon_tx_count",
        "first_exon_order", "last_exon_order", "first_exon_sharing",
        "last_exon_sharing", "first_exon_splicing", "last_exon_splicing",
        "first_exon_overlap", "last_exon_overlap", "merged_txs"
    ]
    exitron_cols = [
        "chr", "start", "end", "frame", "filler", "frame",
        "chr_exon_1", "exon_1_start", "exon_1_end",
        "exon_2_start", "exon_2_end", "gene", "transcript",
        "exon_strand", "exon_frame", "intron_start", "intron_end",
        "intron_length", "merged_txs"
    ]

    # Create empty DataFrames
    single_df = pl.DataFrame({c: [] for c in single_cols})
    multi_df = pl.DataFrame({c: [] for c in multi_cols})
    exitron_df = pl.DataFrame({c: [] for c in exitron_cols})

    return single_df, multi_df, exitron_df, single_exon_file, multi_exon_file, exitron_file

def get_scfr_strand(scfr_df: pl.DataFrame) -> str:
    """
    Returns the strand of a SCFR (first row's strand column).
    """
    if scfr_df.height == 0:
        return "NA"
    return scfr_df[0, "strand"]

def get_total_number_of_exons(scfr_df: pl.DataFrame) -> int:
    """
    Returns total number of exons overlapping a SCFR.
    """
    return scfr_df.height

def get_number_of_merged_exons(scfr_df: pl.DataFrame) -> int:
    """
    Returns number of merged overlapping exons for a SCFR (strand-specific),
    without using any lambda.
    """
    if scfr_df.height == 0:
        return 0

    # Get exon intervals from columns start_right, end_right
    intervals = []
    for row in scfr_df.iter_rows(named=True):
        intervals.append([row["start_right"], row["end_right"]])

    # Sort intervals by start manually
    for i in range(1, len(intervals)):
        key = intervals[i]
        j = i - 1
        while j >= 0 and intervals[j][0] > key[0]:
            intervals[j + 1] = intervals[j]
            j -= 1
        intervals[j + 1] = key

    # Merge overlapping intervals
    merged = []
    for start, end in intervals:
        if not merged:
            merged.append([start, end])
        else:
            last = merged[-1]
            if start <= last[1]:
                last[1] = end if end > last[1] else last[1]
            else:
                merged.append([start, end])

    return len(merged)

def deduplicate_single_exon(scfr_df: pl.DataFrame) -> pl.DataFrame:
    """
    Deduplicate identical single-exon SCFRs.
    Equivalent to shell awk script:
      awk '{ k=$7 FS $8 ...; c[k]++; if(!(k in f)){f[k]=$0;o[++n]=k} }'
    Returns a DataFrame with a 'count' column.
    """
    if scfr_df.height == 0:
        return scfr_df

    # Columns used to define uniqueness (from shell script)
    unique_cols = [
        "start_right", "end_right", "start", "end",
        "exon_strand", "exon_frame", "exon_order",
        "exon_sharing", "upstream_len_in_scfr", "downstream_len_in_scfr",
        "merged_txs"
    ]

    # Group by these columns and count occurrences
    dedup_df = scfr_df.groupby(unique_cols).agg([
        pl.count().alias("count")
    ])

    # Reattach first row values for other columns if needed
    # e.g., chrom, gene, tx, strand, frame
    other_cols = ["chrom", "gene", "tx", "strand", "frame"]
    for col in other_cols:
        dedup_df = dedup_df.with_columns(pl.col(col).first().over(unique_cols))

    return dedup_df

def process_multi_exon_scfr(scfr_df: pl.DataFrame) -> pl.DataFrame:
    """
    Process multi-exon SCFRs (multiple exons per SCFR or transcript).
    Sort exons by start/end based on strand, compute first/last exon info,
    add exon count and other metadata, then deduplicate identical transcripts.
    """
    if scfr_df.height == 0:
        return scfr_df

    # Group by transcript
    transcripts = scfr_df.groupby("tx").agg(pl.all())  # list of rows per transcript

    processed_rows = []

    for tx_group in transcripts.iter_rows(named=True):
        rows = tx_group["chrom"]  # all columns as lists
        # Extract start, end, etc.
        strand = rows[0]["strand"]

        # Sort exons by strand
        if strand == "+":
            sorted_exons = sorted(rows, key=lambda r: r["start"])
        else:
            sorted_exons = sorted(rows, key=lambda r: r["start"], reverse=True)

        # First exon
        first = sorted_exons[0]
        last = sorted_exons[-1]
        exon_count = len(sorted_exons)

        # Build a new row combining first and last exon info (like your paste+awk)
        new_row = {
            "chr": first["chrom"],
            "start": first["start"],
            "end": last["end"],
            "frame": first["frame"],
            "filler": first.get("filler", "1"),
            "strand": strand,
            "first_exon_start": first["start"],
            "last_exon_end": last["end"],
            "gene": first["gene"],
            "transcript": tx_group["tx"][0],
            "exon_strand": strand,
            "exon_frame": first["frame"],
            "exon_count": exon_count,
            "upstream_len_in_scfr": first.get("upstream_len_in_scfr", 0),
            "downstream_len_in_scfr": last.get("downstream_len_in_scfr", 0),
            "first_exon_number": first.get("exon_number", 1),
            "last_exon_number": last.get("exon_number", exon_count),
            "gene_tx_count": first.get("gene_tx_count", 1),
            "first_exon_tx_count": first.get("exon_tx_count", 1),
            "last_exon_tx_count": last.get("exon_tx_count", 1),
            "first_exon_order": first.get("exon_order", "NA"),
            "last_exon_order": last.get("exon_order", "NA"),
            "first_exon_sharing": first.get("exon_sharing", 1.0),
            "last_exon_sharing": last.get("exon_sharing", 1.0),
            "first_exon_splicing": first.get("exon_splicing", "NA"),
            "last_exon_splicing": last.get("exon_splicing", "NA"),
            "first_exon_overlap": first.get("overlap_len", 0),
            "last_exon_overlap": last.get("overlap_len", 0),
            "merged_txs": first.get("merged_txs", "")
        }
        processed_rows.append(new_row)

    multi_df = pl.DataFrame(processed_rows)

    # Deduplicate identical multi-exon transcripts
    unique_cols = [
        "first_exon_start", "last_exon_end", "first_exon_number", "last_exon_number",
        "gene_tx_count", "first_exon_tx_count", "last_exon_tx_count",
        "first_exon_order", "last_exon_order", "first_exon_sharing",
        "last_exon_sharing", "first_exon_splicing", "last_exon_splicing"
    ]

    dedup_df = multi_df.groupby(unique_cols).agg([
        pl.count().alias("count")
    ])

    return dedup_df



def process_scfr_shadow(shadow_csv: str, single_out: str, multi_out: str):
    """
    Python version of the SCFR shadow pipeline:
    - For each unique SCFR, extract overlaps
    - Compute strand, tnoe, mc
    - Classify single/multi exon SCFR
    - Deduplicate
    - Save single & multi exon outputs
    """
    df = pl.read_csv(shadow_csv)
    df = df.with_columns(
        (pl.col("chrom").cast(str) + "_" + pl.col("start").cast(str) + "_" +
         pl.col("end").cast(str) + "_" + pl.col("frame")).alias("scfr_id")
    )

    single_df = pl.DataFrame([])
    multi_df = pl.DataFrame([])

    for scfr_id in df["scfr_id"].unique().to_list():
        scfr_rows = df.filter(pl.col("scfr_id") == scfr_id)
        strand = get_scfr_strand(scfr_rows)
        tnoe = get_total_number_of_exons(scfr_rows)
        mc = get_number_of_merged_exons(scfr_rows)

        if tnoe == 1 or mc == 1:
            processed = deduplicate_single_exon(scfr_rows)
            single_df = pl.concat([single_df, processed])
        elif mc > 1:
            processed = process_multi_exon_scfr(scfr_rows)
            multi_df = pl.concat([multi_df, processed])

    single_df.write_csv(single_out)
    multi_df.write_csv(multi_out)
    print(f"Single-exon SCFRs saved to: {single_out}")
    print(f"Multi-exon SCFRs saved to: {multi_out}")

def calculate_exitrons(scfr_df: pl.DataFrame, strand: str) -> pl.DataFrame:
    """
    Calculate exitron candidates per transcript (like the shell pipeline)
    - Groups by transcript (column 'tx')
    - Sorts exons by start/end based on strand
    - Computes intron_start, intron_end, intron_length for consecutive exons
    - Deduplicates identical introns
    """
    if scfr_df.height == 0:
        return scfr_df

    exitron_rows = []

    # Process each unique transcript
    for tx in scfr_df["tx"].unique().to_list():
        tx_df = scfr_df.filter(pl.col("tx") == tx)
        
        # Sort exons according to strand
        if strand == "+":
            sorted_df = tx_df.sort([pl.col("start"), pl.col("end")])
        elif strand == "-":
            sorted_df = tx_df.sort([pl.col("start").reverse(), pl.col("end").reverse()])
        else:
            continue

        prev_end = None
        prev_start = None

        for row in sorted_df.iter_rows(named=True):
            if prev_end is None:
                prev_end = row["end"]
                prev_start = row["start"]
                continue
            
            # Calculate intron coordinates between consecutive exons
            intron_start = prev_end
            intron_end = row["start"]
            intron_length = abs(intron_end - intron_start)

            exitron_rows.append({
                "chr": row["chrom"],
                "start": row["start"],
                "end": row["end"],
                "frame": row["frame"],
                "filler": row.get("filler", "1"),
                "strand": strand,
                "exon_1_start": prev_start,
                "exon_1_end": prev_end,
                "exon_2_start": row["start"],
                "exon_2_end": row["end"],
                "gene": row["gene"],
                "transcript": tx,
                "exon_strand": row["strand"],
                "exon_frame": row["frame"],
                "intron_start": intron_start,
                "intron_end": intron_end,
                "intron_length": intron_length,
                "merged_txs": row.get("merged_txs", "")
            })

            prev_end = row["end"]
            prev_start = row["start"]

    exitron_df = pl.DataFrame(exitron_rows)

    # Deduplicate identical exitrons
    unique_cols = [
        "chr", "start", "end", "frame", "filler", "strand",
        "exon_1_start", "exon_1_end", "exon_2_start", "exon_2_end",
        "gene", "transcript", "exon_strand", "exon_frame", "intron_start",
        "intron_end", "intron_length"
    ]
    dedup_df = exitron_df.groupby(unique_cols).agg([pl.count().alias("count")])

    return dedup_df

def process_species(base_dir: Path, species: str):
    print(f"\nProcessing {species}")
    out_dir = setup_species_workspace(base_dir, species)

    # SCFR CSV
    scfr_file = base_dir / "SCFR_all" / f"{species}_SCFR_all.out"
    scfr_csv = out_dir / f"{species}_scfr_all.csv"
    scfr_to_bed_csv(scfr_file, scfr_csv)

    # CDS BED CSV
    gtf_dir = base_dir / "genes" / species
    gtf_files = list(gtf_dir.glob("GCF*.gtf"))
    if not gtf_files:
        print(f"No GTF file found for {species}")
        return
    gtf_file = gtf_files[0]

    cds_csv = out_dir / f"{species}_coding_exons.csv"
    gtf_to_bed_pipeline(gtf_file, cds_csv, species, base_dir)

    # Compute SCFR × CDS overlaps
    overlap_csv = out_dir / f"{species}_scfr_cds_all_overlaps.csv"
    compute_scfr_cds_overlaps(scfr_csv, cds_csv, overlap_csv)

    # Filter SCFRs containing exons with same frame and calculate shadow
    shadow_csv = out_dir / f"{species}_scfr_containing_cds.csv"
    filter_containing_cds_with_shadow(overlap_csv, shadow_csv)

     # -------------------------
    # 5. Classify exon shadow: single & multi-exon SCFRs
    # -------------------------
    single_out = out_dir / f"{species}_single_exon.csv"
    multi_out = out_dir / f"{species}_multi_exon.csv"
    process_scfr_shadow(shadow_csv, single_out, multi_out)  # single/multi-exon classification & deduplication

    # -------------------------
    # 6. Calculate exitron candidates
    # -------------------------
    scfr_shadow_df = pl.read_csv(shadow_csv)
    strand = get_scfr_strand(scfr_shadow_df)
    exitron_out = out_dir / f"{species}_exitron_candidates.csv"
    exitron_df = calculate_exitrons(scfr_shadow_df, strand)
    exitron_df.write_csv(exitron_out)

    print(f"Single-exon SCFRs: {single_out}")
    print(f"Multi-exon SCFRs: {multi_out}")
    print(f"Exitron candidates: {exitron_out}")
    print(f"Processing complete for {species}. All CSVs saved in {out_dir}")

# Example usage for all species
if __name__ == "__main__":
    import sys
    from pathlib import Path

    if len(sys.argv) < 3:
        print("Usage: python exon_shadow_pipeline.py <BASE_DIR> <SPECIES|all>")
        sys.exit(1)

    base_dir = Path(sys.argv[1])
    species_arg = sys.argv[2]

    all_species = [
        "human", "bonobo", "chimpanzee", "gorilla",
        "borangutan", "sorangutan", "gibbon"
    ]

    if species_arg == "all":
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor() as executor:
            for sp in all_species:
                executor.submit(process_species, base_dir, sp)
    else:
        process_species(base_dir, species_arg)