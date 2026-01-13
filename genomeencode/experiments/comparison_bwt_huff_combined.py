import time
import matplotlib.pyplot as plt
from genomeencode.burros_wheeler import BurrosWheeler
from genomeencode.huffman import HuffmanTree
from Bio import SeqIO
import os
import psutil


# -----------------------------
# Dataset paths (absolute inside package)
# -----------------------------
DATASET = {
    "Small": "genomeencode/data/small.fasta",
    "Medium": "genomeencode/data/medium.fasta",
    "Large": "genomeencode/data/large.fasta"
}


# -----------------------------
# Load DNA sequence
# -----------------------------
def load_seq(path):
    record = SeqIO.read(path, "fasta")
    return str(record.seq)


# -----------------------------
# Measure Memory Usage
# -----------------------------
def get_memory_kb():
    return psutil.Process().memory_info().rss / 1024


# -----------------------------
# BWT Only (Compression)
# -----------------------------
def run_bwt(seq):
    start_mem = get_memory_kb()
    t0 = time.time()

    bwt = BurrosWheeler.bwt_advanced(seq)

    t1 = time.time()
    end_mem = get_memory_kb()

    return {
        "time": t1 - t0,
        "mem": abs(end_mem - start_mem),
        "size": len(bwt.encode())
    }


# -----------------------------
# Huffman Only (Compression)
# -----------------------------
def run_huffman(seq):
    start_mem = get_memory_kb()
    t0 = time.time()

    h = HuffmanTree(seq)
    h.get_codings(h.root)
    binstr = h.seq_to_binstr()
    unicode_out = HuffmanTree.binstr_to_unicode(binstr)

    t1 = time.time()
    end_mem = get_memory_kb()

    return {
        "time": t1 - t0,
        "mem": abs(end_mem - start_mem),
        "size": len(unicode_out.encode())
    }


# -----------------------------------------
# Combined BWT + Huffman Compression
# -----------------------------------------
def run_combined(seq):
    start_mem = get_memory_kb()
    t0 = time.time()

    bwt = BurrosWheeler.bwt_advanced(seq)
    h = HuffmanTree(bwt)
    h.get_codings(h.root)
    binstr = h.seq_to_binstr()
    unicode_out = HuffmanTree.binstr_to_unicode(binstr)

    t1 = time.time()
    end_mem = get_memory_kb()

    return {
        "time": t1 - t0,
        "mem": abs(end_mem - start_mem),
        "size": len(unicode_out.encode())
    }



# ======================================================
# MAIN EXECUTION
# ======================================================
labels = []
bwt_sizes = []
huff_sizes = []
comb_sizes = []

bwt_speed = []
huff_speed = []
comb_speed = []

bwt_mem = []
huff_mem = []
comb_mem = []


print("\n==== BWT vs HUFFMAN vs COMBINED ====\n")
print(f"{'Dataset':<10} {'BWT_Size':<10} {'Huff_Size':<12} {'Comb_Size':<12} "
      f"{'BWT_Speed':<12} {'Huff_Speed':<12} {'Comb_Speed':<12} "
      f"{'BWT_Mem(KB)':<12} {'Huff_Mem(KB)':<14} {'Comb_Mem(KB)':<14}")
print("-" * 130)


for label, file in DATASET.items():

    seq = load_seq(file)
    labels.append(label)

    # --- Run methods ---
    b = run_bwt(seq)
    h = run_huffman(seq)
    c = run_combined(seq)

    # --- Compression Size ---
    bwt_sizes.append(b["size"])
    huff_sizes.append(h["size"])
    comb_sizes.append(c["size"])

    # --- Speed (bytes/sec) ---
    bwt_speed.append(len(seq) / b["time"])
    huff_speed.append(len(seq) / h["time"])
    comb_speed.append(len(seq) / c["time"])

    # --- Memory ---
    bwt_mem.append(b["mem"])
    huff_mem.append(h["mem"])
    comb_mem.append(c["mem"])

    # PRINT TABLE ROW
    print(f"{label:<10} {b['size']:<10} {h['size']:<12} {c['size']:<12} "
          f"{(len(seq)/b['time']):<12.2f} {(len(seq)/h['time']):<12.2f} {(len(seq)/c['time']):<12.2f} "
          f"{b['mem']:<12.2f} {h['mem']:<14.2f} {c['mem']:<14.2f}")


# =====================================================
# GRAPH: Combined Comparison
# =====================================================

plt.figure(figsize=(14, 8))

# Compression Size
plt.subplot(3, 1, 1)
plt.plot(labels, bwt_sizes, marker='o', label="BWT Size")
plt.plot(labels, huff_sizes, marker='o', label="Huffman Size")
plt.plot(labels, comb_sizes, marker='o', label="Combined Size")
plt.title("Compression Size Comparison")
plt.ylabel("Size (bytes)")
plt.legend()
plt.grid(True)

# Speed
plt.subplot(3, 1, 2)
plt.plot(labels, bwt_speed, marker='o', label="BWT Speed")
plt.plot(labels, huff_speed, marker='o', label="Huffman Speed")
plt.plot(labels, comb_speed, marker='o', label="Combined Speed")
plt.title("Compression Speed Comparison")
plt.ylabel("Bytes/sec")
plt.yscale('log')
plt.legend()
plt.grid(True)

# Memory
plt.subplot(3, 1, 3)
plt.plot(labels, bwt_mem, marker='o', label="BWT Memory")
plt.plot(labels, huff_mem, marker='o', label="Huffman Memory")
plt.plot(labels, comb_mem, marker='o', label="Combined Memory")
plt.title("Memory Usage Comparison (KB)")
plt.ylabel("Memory (KB)")
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.savefig("bwt_huffman_combined_comparison.png", dpi=300)
plt.show()
