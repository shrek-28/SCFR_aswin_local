import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
from collections import defaultdict

# Codon table (Standard)
codon_table = {
    'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L',
    'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
    'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M',
    'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
    'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',
    'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
    'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
    'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
    'TAT':'Y', 'TAC':'Y', 'TAA':'*', 'TAG':'*',
    'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q',
    'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
    'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
    'TGT':'C', 'TGC':'C', 'TGA':'*', 'TGG':'W',
    'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
    'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
    'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'
}

# Synonymous codons
amino_to_codons = defaultdict(list)
for codon, aa in codon_table.items():
    if aa != "*":
        amino_to_codons[aa].append(codon)

def is_valid_orf(seq):
    if len(seq) % 3 != 0:
        return False
    if not seq.startswith("ATG") or seq[-3:] not in ["TAA", "TAG", "TGA"]:
        return False
    protein = Seq(seq).translate(to_stop=True)
    return '*' not in protein[:-1]

def codon_usage(seq):
    usage = defaultdict(int)
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        if codon in codon_table:
            usage[codon] += 1
    return usage
def calculate_rscu(codon_counts):
    rscu = {}
    for aa, codons in amino_to_codons.items():
        total = sum(codon_counts[c] for c in codons)
        if total == 0:
            for c in codons:
                rscu[c] = 0.0
        else:
            expected = total / len(codons)
            for c in codons:
                rscu[c] = codon_counts[c] / expected if expected else 0.0
    return rscu

def sliding_metrics(seq, window=60, step=30):
    metrics = []
    for i in range(0, len(seq) - window + 1, step):
        sub = seq[i:i+window]
        gc = (sub.count('G') + sub.count('C')) / len(sub)
        gc_stretches = sum(1 for x in ['G', 'C'] if x*5 in sub)
        metrics.append((i, i+window, gc, gc_stretches))
    return metrics

def process_directory(fasta_dir, summary_out="global_metrics.tsv", region_out="regional_metrics.tsv", rscu_out="rscu.tsv"):
    global_data = []
    regional_data = []
    rscu_data = []

    for file in os.listdir(fasta_dir):
        if file.endswith(".fa") or file.endswith(".fasta"):
            clade = os.path.splitext(file)[0]
            path = os.path.join(fasta_dir, file)
            for record in SeqIO.parse(path, "fasta"):
                species = record.id
                seq = str(record.seq).upper().replace('-', '')
                orf_ok = is_valid_orf(seq)
                codons = codon_usage(seq)
                gc_content = (seq.count('G') + seq.count('C')) / len(seq)
                gc_stretch_count = sum(1 for x in ['G', 'C'] if x*5 in seq)

                # Global codon info
                global_data.append({
                    'clade': clade, 'species': species, 'length': len(seq),
                    'orf_valid': orf_ok, 'gc_content': round(gc_content, 4),
                    'gc_stretch_count': gc_stretch_count,
                    **{f'codon_{k}': v for k, v in codons.items()}
                })

                # Sliding window
                for start, end, gc, stretch in sliding_metrics(seq):
                    regional_data.append({
                        'clade': clade, 'species': species,
                        'start': start, 'end': end,
                        'gc_content': round(gc, 4),
                        'gc_stretch_count': stretch
                    })

                # RSCU calculation
                rscu = calculate_rscu(codons)
                for codon, val in rscu.items():
                    rscu_data.append({
                        'clade': clade,
                        'species': species,
                        'codon': codon,
                        'amino_acid': codon_table[codon],
                        'rscu': round(val, 4)
                    })

    pd.DataFrame(global_data).to_csv(summary_out, sep="\t", index=False)
    pd.DataFrame(regional_data).to_csv(region_out, sep="\t", index=False)
    pd.DataFrame(rscu_data).to_csv(rscu_out, sep="\t", index=False)

# ---------- MAIN ----------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute codon usage metrics from FASTA directory")
    parser.add_argument("input_dir", help="Path to folder containing FASTA files")
    parser.add_argument("output_dir", help="Directory to save TSV output files")
    args = parser.parse_args()

    input_dir = args.input_dir.rstrip("/")
    outdir = args.output_dir.rstrip("/")

    # Create output directory if missing
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    summary_out = os.path.join(outdir, "global_metrics.tsv")
    region_out = os.path.join(outdir, "regional_metrics.tsv")
    rscu_out = os.path.join(outdir, "rscu.tsv")

    process_directory(input_dir, summary_out, region_out, rscu_out)
