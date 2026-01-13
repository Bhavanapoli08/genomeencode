import os
import time
import matplotlib.pyplot as plt
import psutil
from Bio import SeqIO
from genomeencode.burros_wheeler import BurrosWheeler
from genomeencode.huffman import HuffmanTree

# --------------------------------------------
# LOAD FASTA SEQUENCE
# --------------------------------------------
def load_seq(path):
    record = SeqIO.read(path, "fasta")
    return str(record.seq)

DATASET_PATH = os.path.join(
    os.path.dirname(__file__), "..", "data"
)

DATASET = {
    "Small":  os.path.join(DATASET_PATH, "small.fasta"),
    "Medium": os.path.join(DATASET_PATH, "medium.fasta"),
    "Large":  os.path.join(DATASET_PATH, "large.fasta")
}

# --------------------------------------------
# MEMORY MEASUREMENT WRAPPER
# --------------------------------------------
def memory_used_kb(func, *args):
    process = psutil.Process(os.getpid())
    mem_before = process.memory_info().rss
    func(*args)
    mem_after = process.memory_info().rss
    return (mem_after - mem_before) / 1024  # KB


# --------------------------------------------
# MEASURE BWT MEMORY
# --------------------------------------------
def measure_bwt_memory(seq):
    return memory_used_kb(BurrosWheeler.bwt_advanced, seq)


# --------------------------------------------
# MEASURE HUFFMAN MEMORY
# --------------------------------------------
def measure_huffman_memory(seq):
    def huff_task(s):
        h = HuffmanTree(s)
        h.get_codings(h.root)
        _ = h.seq_to_binstr()
    return memory_used_kb(huff_task, seq)


# --------------------------------------------
# MEASURE COMBINED MEMORY
# --------------------------------------------
def measure_combined_memory(seq):
    def combined_task(s):
        t = BurrosWheeler.bwt_advanced(s)
        h = HuffmanTree(t)
        h.get_codings(h.root)
        _ = h.seq_to_binstr()
    return memory_used_kb(combined_task, seq)


# --------------------------------------------
# MAIN
# --------------------------------------------
sizes = []
bwt_mem = []
huff_mem = []
comb_mem = []

print("\n==== MEMORY USAGE (KB) ====\n")
print(f"{'Dataset':<10} {'BWT(KB)':<12} {'Huff(KB)':<12} {'Combined(KB)':<15}")
print("-" * 60)

for label, path in DATASET.items():
    seq = load_seq(path)

    m1 = measure_bwt_memory(seq)
    m2 = measure_huffman_memory(seq)
    m3 = measure_combined_memory(seq)

    sizes.append(label)
    bwt_mem.append(m1)
    huff_mem.append(m2)
    comb_mem.append(m3)

    print(f"{label:<10} {m1:<12.2f} {m2:<12.2f} {m3:<15.2f}")


# --------------------------------------------
# DUAL AXIS PLOT (so orange is visible)
# --------------------------------------------
plt.figure(figsize=(12, 6))
fig, ax1 = plt.subplots(figsize=(14, 6))

# Left axis – BWT + Combined
ax1.set_xlabel("Dataset Size")
ax1.set_ylabel("Memory (KB) - BWT & Combined")
l1 = ax1.plot(sizes, bwt_mem, marker='o', label="BWT Memory", color='blue')
l3 = ax1.plot(sizes, comb_mem, marker='o', label="Combined Memory", color='green')
ax1.tick_params(axis='y')

# Right axis – Huffman
ax2 = ax1.twinx()
ax2.set_ylabel("Memory (KB) - Huffman")
l2 = ax2.plot(sizes, huff_mem, marker='o', label="Huffman Memory", color='orange')
ax2.tick_params(axis='y')

# Combine legends
lines = l1 + l2 + l3
labels = [l.get_label() for l in lines]
ax1.legend(lines, labels, loc="upper left")

plt.title("Memory Usage Comparison (KB)")
plt.grid(True)
plt.tight_layout()
plt.show()

