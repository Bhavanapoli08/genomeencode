import os
from Bio import SeqIO
import matplotlib.pyplot as plt
from genomeencode.burros_wheeler import BurrosWheeler
from genomeencode.huffman import HuffmanTree

# ============================
# DATASET PATHS
# ============================


BASE_DIR = os.path.dirname(os.path.dirname(__file__))   # genomeencode/
DATA_DIR = os.path.join(BASE_DIR, "data")               # genomeencode/data

DATASET = {
    "Small":  os.path.join(DATA_DIR, "small.fasta"),
    "Medium": os.path.join(DATA_DIR, "medium.fasta"),
    "Large":  os.path.join(DATA_DIR, "large.fasta")
}

# ============================
# LOAD FASTA FILE
# ============================

def load_seq(path):
    record = SeqIO.read(path, "fasta")
    return str(record.seq)


# ============================
# COMPUTE BWT SIZE
# ============================

def compute_bwt_size(seq):
    t = BurrosWheeler.bwt_advanced(seq)
    return len(t.encode("utf-8"))


# ============================
# COMPUTE HUFFMAN SIZE
# ============================

def compute_huffman_size(seq):
    huff = HuffmanTree(seq)
    huff.get_codings(huff.root)

    binary = huff.seq_to_binstr()
    unicode_out = HuffmanTree.binstr_to_unicode(binary)

    return len(unicode_out.encode("utf-8"))


# ============================
# COMPUTE BWT + HUFFMAN SIZE
# ============================

def compute_combined_size(seq):
    t = BurrosWheeler.bwt_advanced(seq)

    huff = HuffmanTree(t)
    huff.get_codings(huff.root)

    binary = huff.seq_to_binstr()
    unicode_out = HuffmanTree.binstr_to_unicode(binary)

    return len(unicode_out.encode("utf-8"))


# ============================
# MAIN EXECUTION
# ============================

print("\n==== COMPRESSION RATIO METRICS ====\n")
print(f"{'Dataset':<10} {'Orig(Bytes)':<12} {'BWT(Bytes)':<12} {'Huff(Bytes)':<12} {'Comb(Bytes)':<12} {'BWT_Ratio':<12} {'Huff_Ratio':<12} {'Comb_Ratio':<12}")
print("-" * 110)

labels = []
orig_sizes = []
bwt_sizes = []
huff_sizes = []
comb_sizes = []

for label, file in DATASET.items():
    seq = load_seq(file)
    orig = len(seq.encode("utf-8"))

    bwt = compute_bwt_size(seq)
    huff = compute_huffman_size(seq)
    comb = compute_combined_size(seq)

    labels.append(label)
    orig_sizes.append(orig)
    bwt_sizes.append(bwt)
    huff_sizes.append(huff)
    comb_sizes.append(comb)

    print(f"{label:<10} {orig:<12} {bwt:<12} {huff:<12} {comb:<12} "
          f"{(bwt/orig):<12.3f} {(huff/orig):<12.3f} {(comb/orig):<12.3f}")


# ============================
# BAR GRAPH VISUALIZATION
# ============================

plt.figure(figsize=(12, 6))
bar_width = 0.25
x = range(len(labels))

plt.bar([n - bar_width for n in x], bwt_sizes, width=bar_width, label="BWT")
plt.bar(x, huff_sizes, width=bar_width, label="Huffman")
plt.bar([n + bar_width for n in x], comb_sizes, width=bar_width, label="BWT + Huffman")

plt.xticks(x, labels)
plt.ylabel("Compressed Size (Bytes)")
plt.title("Compression Size Comparison (BWT vs Huffman vs Combined)")
plt.legend()
plt.grid(axis="y")
plt.tight_layout()

plt.savefig("compression_ratio_graph.png", dpi=300)
plt.show()

