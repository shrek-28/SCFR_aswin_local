#!/usr/bin/env python3
import sys
import os
import csv
import re
import numpy as np
from collections import Counter
from Bio import SeqIO
from scipy.signal import find_peaks

# --------------------------------------------------------
# Robust SCFR extractor: parse basename, handle NC_..._1 -> NC_... .1
# --------------------------------------------------------
def extract_scfr_sequence(multifasta, csv_filename):
    """
    Given a multifasta and a CSV filename like:
      SCFR_NC_060925_1_110130430_110135794_frame2_fft.csv
    find the matching FASTA record whose description contains:
      NC_060925.1:110130430-110135794
    Returns the sequence (uppercase string).
    """
    basename = os.path.basename(csv_filename)

    # Find two reasonably large integers in filename → assume they are start and end
    nums = re.findall(r'(\d{5,12})', basename)
    if len(nums) < 2:
        raise ValueError(f"ERROR: Could not find start/end coordinates in filename: {basename}")
    # heuristics: take the last two numbers as start,end
    start, end = nums[-2], nums[-1]

    # Attempt to find accession token: prefer patterns like NC_060925_1 or NC_060925.1
    acc_match = re.search(r'(NC[_\.]\d+(?:[_\.]\d+)?)', basename, flags=re.IGNORECASE)
    if not acc_match:
        # fallback: any NC... token (less strict)
        acc_match = re.search(r'(NC[A-Za-z0-9_.-]+)', basename, flags=re.IGNORECASE)
    if not acc_match:
        raise ValueError(f"ERROR: Could not find accession token in filename: {basename}")

    acc_raw = acc_match.group(1)

    # Convert underscore-style accession NC_060925_1 -> NC_060925.1 to match fasta header
    m_acc = re.match(r'(NC)_(\d+)_(\d+)$', acc_raw, flags=re.IGNORECASE)
    if m_acc:
        acc_key = f"{m_acc.group(1)}_{m_acc.group(2)}.{m_acc.group(3)}"
    else:
        # if already in dot form (NC_060925.1) or NC060925.1 etc., keep as-is
        # also normalize any accidental underscores before digits that aren't the version number:
        acc_key = acc_raw.replace('__', '_')

    target_key = f"{acc_key}:{start}-{end}"
    # debug print
    # print("DEBUG: basename:", basename)
    # print("DEBUG: extracted acc_raw:", acc_raw, "-> acc_key:", acc_key)
    # print("DEBUG: start,end:", start, end)
    # print("DEBUG: target_key:", target_key)

    for record in SeqIO.parse(multifasta, "fasta"):
        if target_key in record.description:
            return str(record.seq).upper()

    # second attempt: some FASTA headers may not have the accession version (.1) — try with base accession
    acc_base = acc_key.split('.')[0]
    alt_key = f"{acc_base}:{start}-{end}"
    for record in SeqIO.parse(multifasta, "fasta"):
        if alt_key in record.description:
            return str(record.seq).upper()

    raise ValueError(f"ERROR: SCFR sequence {target_key} (or {alt_key}) not found in {multifasta}")

# --------------------------------------------------------
# Extract top motifs
# --------------------------------------------------------
def extract_motifs(seq, period, top_n=5):
    motifs = [seq[i:i+period] for i in range(0, len(seq) - period + 1)]
    return Counter(motifs).most_common(top_n)

# --------------------------------------------------------
# Load FFT CSV
# --------------------------------------------------------
def load_fft_csv(path):
    freqs, mags = [], []
    with open(path) as f:
        r = csv.reader(f)
        next(r)
        for row in r:
            freqs.append(float(row[0]))
            mags.append(float(row[1]))
    return np.array(freqs), np.array(mags)

# --------------------------------------------------------
# Table printer
# --------------------------------------------------------
def print_table(rows, headers):
    if not rows:
        print("(no rows to display)")
        return
    col_widths = [max(len(str(row[i])) for row in rows + [headers]) for i in range(len(headers))]
    fmt = "  ".join("{:<" + str(w) + "}" for w in col_widths)
    print(fmt.format(*headers))
    print("-" * (sum(col_widths) + 2 * (len(headers) - 1)))
    for row in rows:
        print(fmt.format(*row))

# --------------------------------------------------------
# MAIN
# --------------------------------------------------------
def main():

    if len(sys.argv) != 3:
        print("Usage: python3 scfr_fft_peak_motif_report.py FFT.csv multifasta.fasta")
        sys.exit(1)

    csv_file = sys.argv[1]
    multifasta = sys.argv[2]

    print(f"\nExtracting SCFR sequence for file: {os.path.basename(csv_file)} from {multifasta} ...")
    try:
        seq = extract_scfr_sequence(multifasta, csv_file)
    except Exception as e:
        print(str(e))
        sys.exit(2)

    N = len(seq)
    print(f"SCFR length = {N} bp")

    freqs, mags = load_fft_csv(csv_file)

    # Peak detection parameters (same as original pipeline)
    half = N // 2
    mean_mag = np.mean(mags[:half])
    std_mag = np.std(mags[:half])
    threshold = mean_mag + 2 * std_mag

    peaks, props = find_peaks(mags[:half], height = threshold)

    if "peak_heights" in props:
        heights = props["peak_heights"]
    else:
        heights = np.array([])

    if len(heights) == 0:
        print("\nNo peaks above threshold (mean + 2*std). Consider relaxing threshold.\n")
        sys.exit(0)

    top = np.argsort(heights)[::-1][:3]
    selected_indices = [peaks[i] for i in top if i < len(peaks)]

    rows = []
    for idx in selected_indices:
        if idx == 0:
            continue
        freq = freqs[idx]
        mag = mags[idx]
        period = int(round(N / idx))
        if 2 <= period <= 30:
            motif_list = extract_motifs(seq, period)
            motif_str = ", ".join(f"{m}:{c}" for m, c in motif_list)
        else:
            motif_str = "No motifs (period out of range)"
        rows.append([f"{freq:.6f}", f"{mag:.3f}", str(period), motif_str])

    print("\n=== Peak & Motif Report (SCFR) ===\n")
    print_table(rows, ["Frequency", "Magnitude", "Period", "Top Motifs"])
    print("\nDone.\n")

if __name__ == "__main__":
    main()

