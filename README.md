import matplotlib.pyplot as plt
import pandas as pd
from collections import Counter

# -------- STEP 1: Read FASTA --------
def read_fasta(file_path, n_genes=20):
    genes = []
    seq = ""

    with open(file_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq:
                    genes.append(seq)
                    if len(genes) == n_genes:
                        break
                    seq = ""
            else:
                seq += line
    if len(genes) < n_genes and seq:
        genes.append(seq)

    return genes


# -------- STEP 2: Basic calculations --------
def gc_percent(seq):
    return (seq.count("G") + seq.count("C")) / len(seq) * 100


def kmer_count(seq, k):
    return Counter(seq[i:i+k] for i in range(len(seq) - k + 1))


def codon_usage(seq):
    return Counter(seq[i:i+3] for i in range(0, len(seq) - 2, 3))


# -------- STEP 3: Load data --------
file_path = r"C:\Users\Admin\Downloads\E-coli\ncbi_dataset\data\GCF_000005845.2\cds_from_genomic.fna"

genes = read_fasta(file_path, n_genes=20)

# -------- STEP 4: Compute metrics --------
lengths = []
gc_values = []
all_2mers = Counter()
all_3mers = Counter()
all_codons = Counter()

for seq in genes:
    lengths.append(len(seq))
    gc_values.append(gc_percent(seq))
    all_2mers.update(kmer_count(seq, 2))
    all_3mers.update(kmer_count(seq, 3))
    all_codons.update(codon_usage(seq))

# -------- STEP 5: Make table --------
df = pd.DataFrame({
    "Gene": [f"Gene_{i+1}" for i in range(len(genes))],
    "Length": lengths,
    "GC_percent": gc_values
})

print(df)

# -------- STEP 6: Plots --------
plt.bar(df["Gene"], df["Length"])
plt.xticks(rotation=90)
plt.title("Gene Lengths")
plt.ylabel("Base pairs")
plt.show()

plt.bar(df["Gene"], df["GC_percent"])
plt.xticks(rotation=90)
plt.title("GC Content (%)")
plt.ylabel("GC %")
plt.show()

# -------- STEP 7: Top k-mers & codons --------
print("\nTop 10 2-mers:")
print(all_2mers.most_common(10))

print("\nTop 10 3-mers:")
print(all_3mers.most_common(10))

print("\nTop 10 Codons:")
print(all_codons.most_common(10))

# -------- STEP 8: Genome-wide GC comparison --------

all_genes = read_fasta(file_path, n_genes=100000)  # read all genes

genome_gc = [gc_percent(seq) for seq in all_genes]
genome_gc_avg = sum(genome_gc) / len(genome_gc)

print("\nGenome-wide average GC%:", round(genome_gc_avg, 2))
print("Average GC% of selected 20 genes:", round(sum(gc_values)/len(gc_values), 2))

# -------- STEP 9: Heatmap of 2-mer frequencies --------

import numpy as np

bases = ["A", "T", "G", "C"]
all_2mer_keys = [a+b for a in bases for b in bases]

heatmap_data = []

for seq in genes:
    counts = kmer_count(seq, 2)
    total = sum(counts.values())
    heatmap_data.append([counts[k]/total for k in all_2mer_keys])

heatmap_data = np.array(heatmap_data)

plt.figure(figsize=(10,6))
plt.imshow(heatmap_data, aspect="auto")
plt.colorbar(label="Normalized frequency")
plt.xticks(range(len(all_2mer_keys)), all_2mer_keys, rotation=90)
plt.yticks(range(len(genes)), df["Gene"])
plt.title("2-mer frequency heatmap (20 E. coli genes)")
plt.show()

# -------- STEP 10: PCA on 2-mer composition --------

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

X = heatmap_data  # reuse normalized 2-mer matrix

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

pca = PCA(n_components=2)
X_pca = pca.fit_transform(X_scaled)

plt.figure()
plt.scatter(X_pca[:,0], X_pca[:,1])

for i, gene in enumerate(df["Gene"]):
    plt.text(X_pca[i,0], X_pca[i,1], gene)

plt.xlabel("PC1")
plt.ylabel("PC2")
plt.title("PCA of 2-mer composition (PC1 vs PC2)")
plt.show()

print("\nExplained variance by PC1 and PC2:")
print(pca.explained_variance_ratio_)
