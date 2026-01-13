import time
import matplotlib.pyplot as plt
from Bio import SeqIO
from genomeencode.burros_wheeler import BurrosWheeler
from genomeencode.huffman import HuffmanTree
import os

# -----------------------------
# Dataset Paths
# -----------------------------
DATASET = {
    "Small":  "genomeencode/data/small.fasta",
    "Medium": "genomeencode/data/medium.fasta",
    "Large":  "genomeencode/data/large.fasta"
}


# -----------------------------
# Load FASTA
# -----------------------------
def load_seq(path):
    record = SeqIO.read(path, "fasta")
    return str(record.seq)


# -----------------------------
# Compression Functions
# -----------------------------
def time_bwt(seq):
    start = time.time()
    out = BurrosWheeler.bwt_advanced(seq)
    end = time.time()
    return end - start, out


def time_huffman(seq):
    start = time.time()
    huff = HuffmanTree(seq)
    huff.get_codings(huff.root)
    binstr = huff.seq_to_binstr()
    unicode_out = HuffmanTree.binstr_to_unicode(binstr)
    end = time.time()
    return end - start, unicode_out


def time_combined(seq):
    start = time.time()

    # BWT
    t = BurrosWheeler.bwt_advanced(seq)

    # Huffman
    huff = HuffmanTree(t)
    huff.get_codings(huff.root)
    binstr = huff.seq_to_binstr()
    unicode_out = HuffmanTree.binstr_to_unicode(binstr)

    end = time.time()
    return end - start, unicode_out


# -----------------------------
# Metrics Calculation
# -----------------------------
sizes = []
orig_sizes = []
bwt_speed = []
huff_speed = []
comb_speed = []


print("\n==== COMPRESSION SPEED METRICS (Bytes/sec) ====\n")

print(f"{'Dataset':<10} {'Orig(Bytes)':<12} {'BWT_Speed':<12} {'Huff_Speed':<12} {'Comb_Speed':<12}")
print("-" * 60)

for label, file in DATASET.items():

    seq = load_seq(file)
    orig_bytes = len(seq.encode("utf-8"))

    # BWT
    bwt_time, bwt_out = time_bwt(seq)
    bwt_spd = orig_bytes / bwt_time

    # Huffman
    huff_time, huff_out = time_huffman(seq)
    huff_spd = orig_bytes / huff_time

    # Combined
    comb_time, comb_out = time_combined(seq)
    comb_spd = orig_bytes / comb_time

    # Save for plotting
    sizes.append(label)
    orig_sizes.append(orig_bytes)
    bwt_speed.append(bwt_spd)
    huff_speed.append(huff_spd)
    comb_speed.append(comb_spd)

    print(f"{label:<10} {orig_bytes:<12} {bwt_spd:<12.2f} {huff_spd:<12.2f} {comb_spd:<12.2f}")


# -----------------------------
# Plot Compression Speed
# -----------------------------
plt.figure(figsize=(10, 6))
plt.plot(sizes, bwt_speed, marker='o', label="BWT Only")
plt.plot(sizes, huff_speed, marker='o', label="Huffman Only")
plt.plot(sizes, comb_speed, marker='o', label="BWT + Huffman")

plt.title("Compression Speed Comparison (Bytes/sec)")
plt.xlabel("Dataset Size")
plt.ylabel("Speed (Bytes/sec)")
plt.grid(True)
plt.legend()
plt.tight_layout()

plt.savefig("compression_speed_graph.png", dpi=300)
plt.show()
