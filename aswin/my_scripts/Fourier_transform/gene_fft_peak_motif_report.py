import sys
import csv
import numpy as np
from collections import Counter
from Bio import SeqIO
from scipy.signal import find_peaks

# --------------------------------------------------------
# Extract FASTA sequence from multipfast using protein_id
# --------------------------------------------------------
def extract_fasta_from_multifasta(multi_fasta, protein_id):
    """
    Extracts a sequence whose header contains protein_id=NP_xxxx
    protein_id: NP_001278296   (dot version also allowed)
    """
    base = protein_id.split(".")[0]   # NP_001278296

    for record in SeqIO.parse(multi_fasta, "fasta"):
        header = record.description
        if f"protein_id={protein_id}" in header or f"protein_id={base}" in header:
            return str(record.seq).upper()

    raise ValueError(f"ERROR: Protein ID {protein_id} not found in {multi_fasta}")

# --------------------------------------------------------
# Extract top 5 motifs for a given periodicity
# --------------------------------------------------------
def extract_motifs(seq, period, top_n=5):
    motifs = [seq[i:i+period] for i in range(0, len(seq) - period + 1)]
    return Counter(motifs).most_common(top_n)

# --------------------------------------------------------
# Load FFT CSV (Freq, Magnitude)
# --------------------------------------------------------
def load_fft_csv(path):
    freqs = []
    mags = []
    with open(path) as f:
        r = csv.reader(f)
        next(r)
        for row in r:
            freqs.append(float(row[0]))
            mags.append(float(row[1]))
    return np.array(freqs), np.array(mags)

# --------------------------------------------------------
# Pretty table printer
# --------------------------------------------------------
def print_table(rows, headers):
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
    if len(sys.argv) != 4:
        print("Usage: python3 analyze_fft_csv.py FFT.csv PROTEIN_ID multifasta.fasta")
        sys.exit(1)

    csv_file = sys.argv[1]
    protein_id = sys.argv[2]
    multifasta = sys.argv[3]

    # ----- Extract FASTA from multi-fasta -----
    print(f"\nExtracting protein '{protein_id}' from {multifasta} ...")
    seq = extract_fasta_from_multifasta(multifasta, protein_id)
    N = len(seq)
    print(f"Sequence length = {N} bp")

    # ----- read CSV FFT -----
    freqs, mags = load_fft_csv(csv_file)

    # ----- Identify significant peaks exactly like PDF pipeline -----
    half = N // 2
    mean_mag = np.mean(mags[:half])
    std_mag = np.std(mags[:half])

    peaks, props = find_peaks(mags[:half], height=mean_mag + 2 * std_mag)
    heights = props["peak_heights"]

    # Pick top 3 peaks
    top = np.argsort(heights)[::-1][:3]
    selected_indices = [peaks[i] for i in top]

    # --------------------------------------------------------
    # Compute motifs for each selected frequency
    # --------------------------------------------------------
    rows = []
    for idx in selected_indices:
        freq = freqs[idx]
        mag = mags[idx]

        if idx == 0:
            continue

        period = int(round(N / idx))

        # Filter for meaningful motif sizes
        if 2 <= period <= 30:
            motif_list = extract_motifs(seq, period)
        else:
            motif_list = []

        motif_str = ", ".join(f"{m}:{c}" for m, c in motif_list) if motif_list else "No motifs"

        rows.append([
            f"{freq:.6f}",
            f"{mag:.3f}",
            str(period),
            motif_str
        ])

    # --------------------------------------------------------
    # Print output as a neat table
    # --------------------------------------------------------
    print("\n=== Peak & Motif Report ===\n")
    print_table(rows, ["Frequency", "Magnitude", "Period", "Top Motifs"])

    print("\nDone.\n")

if __name__ == "__main__":
    main()

