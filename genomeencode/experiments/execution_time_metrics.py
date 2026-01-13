import time
import matplotlib.pyplot as plt
from genomeencode.burros_wheeler import BurrosWheeler
from genomeencode.huffman import HuffmanTree
from genomeencode.sequence import Sequence
from Bio import SeqIO
import os

# ----------------------------------------------------
# DATASET PATHS (Corrected)
# ----------------------------------------------------

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, "..", "data")

DATASET = {
    "Small":  os.path.join(DATA_DIR, "small.fasta"),
    "Medium": os.path.join(DATA_DIR, "medium.fasta"),
    "Large":  os.path.join(DATA_DIR, "large.fasta"),
}


# ----------------------------------------------------
# LOAD FASTA SEQUENCE
# ----------------------------------------------------
def load_seq(path):
    record = SeqIO.read(path, "fasta")
    return str(record.seq)


# ----------------------------------------------------
# MEASURE BWT TIME
# ----------------------------------------------------
def measure_bwt(seq):
    start = time.time()
    _ = BurrosWheeler.bwt_advanced(seq)
    return time.time() - start


# ----------------------------------------------------
# MEASURE HUFFMAN TIME
# ----------------------------------------------------
def measure_huffman(seq):
    start = time.time()

    huff = HuffmanTree(seq)
    huff.get_codings(huff.root)
    binary = huff.seq_to_binstr()
    _ = HuffmanTree.binstr_to_unicode(binary)

    return time.time() - start


# ----------------------------------------------------
# MEASURE BWT + HUFFMAN (COMBINED)
# ----------------------------------------------------
def measure_combined(seq):
    start = time.time()

    # BWT stage
    t = BurrosWheeler.bwt_advanced(seq)

    # Huffman stage
    huff = HuffmanTree(t)
    huff.get_codings(huff.root)
    binary = huff.seq_to_binstr()
    _ = HuffmanTree.binstr_to_unicode(binary)

    return time.time() - start


# ----------------------------------------------------
# MAIN EXECUTION
# ----------------------------------------------------
sizes = []
bwt_times = []
huff_times = []
combined_times = []

print("\n=== Running Execution Time Experiments ===\n")

for label, filepath in DATASET.items():

    seq = load_seq(filepath)
    sizes.append(label)

    print(f"\nProcessing {label} dataset... ({len(seq)} bases)")

    bwt_t = measure_bwt(seq)
    huff_t = measure_huffman(seq)
    comb_t = measure_combined(seq)

    bwt_times.append(bwt_t)
    huff_times.append(huff_t)
    combined_times.append(comb_t)

    print(f"  BWT Time        : {bwt_t:.6f} sec")
    print(f"  Huffman Time    : {huff_t:.6f} sec")
    print(f"  BWT + Huffman   : {comb_t:.6f} sec")


# ----------------------------------------------------
# PRINT TABLE
# ----------------------------------------------------
print("\n\n========== EXECUTION TIME TABLE ==========\n")
print(f"{'Dataset':<10} {'BWT':<12} {'Huffman':<12} {'BWT+Huffman':<12}")
print("-" * 50)

for i in range(len(sizes)):
    print(f"{sizes[i]:<10} {bwt_times[i]:<12.6f} {huff_times[i]:<12.6f} {combined_times[i]:<12.6f}")


# ----------------------------------------------------
# PLOT GRAPH
# ----------------------------------------------------
plt.figure(figsize=(10,6))
plt.plot(sizes, bwt_times, marker='o', linewidth=2, label="BWT Only")
plt.plot(sizes, huff_times, marker='o', linewidth=2, label="Huffman Only")
plt.plot(sizes, combined_times, marker='o', linewidth=2, label="BWT + Huffman")

plt.title("Execution Time Comparison", fontsize=16)
plt.xlabel("Dataset Size")
plt.ylabel("Execution Time (seconds)")
plt.grid(True)
plt.legend()
plt.tight_layout()

OUTPUT_PATH = os.path.join(BASE_DIR, "execution_time_graph.png")
plt.savefig(OUTPUT_PATH, dpi=300)

print(f"\nGraph saved to: {OUTPUT_PATH}\n")
plt.show()

