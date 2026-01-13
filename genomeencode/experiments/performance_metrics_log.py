import time
import matplotlib.pyplot as plt
from genomeencode.burros_wheeler import BurrosWheeler
from genomeencode.huffman import HuffmanTree
from genomeencode.sequence import Sequence
from genomeencode.experiments.load_sequences import load_fasta

# ---------------------------
# Load DNA Sequences
# ---------------------------

small = load_fasta("small.fasta")
medium = load_fasta("medium.fasta")
large = load_fasta("large.fasta")

sequences = [small, medium, large]
labels = ["Short DNA", "Medium DNA", "Large DNA"]

# ---------------------------
# Functions to measure time
# ---------------------------

def time_bwt(seq):
    start = time.time()
    BurrosWheeler.bwt_advanced(seq)
    return time.time() - start

def time_huffman(seq):
    start = time.time()
    tree = HuffmanTree(seq)
    tree.get_codings(tree.root)
    tree.seq_to_binstr()
    return time.time() - start

def time_bwt_huffman(seq):
    start = time.time()
    t = BurrosWheeler.bwt_advanced(seq)
    tree = HuffmanTree(t)
    tree.get_codings(tree.root)
    tree.seq_to_binstr()
    return time.time() - start

# ---------------------------
# Measure Times
# ---------------------------

bwt_times = []
huffman_times = []
combined_times = []

for seq in sequences:
    bwt_times.append(time_bwt(seq))
    huffman_times.append(time_huffman(seq))
    combined_times.append(time_bwt_huffman(seq))

# Add tiny offset so Huffman doesn't become flat
huffman_times = [t + 0.00001 for t in huffman_times]

# ---------------------------
# Plot (Logarithmic Scale)
# ---------------------------

plt.figure(figsize=(10, 6))

x = range(len(labels))

plt.plot(x, bwt_times, marker='o', label="BWT", linewidth=2)
plt.plot(x, huffman_times, marker='o', label="Huffman", linewidth=2)
plt.plot(x, combined_times, marker='o', label="BWT + Huffman", linewidth=2)

plt.xticks(x, labels)
plt.yscale("log")   # <<< key change (Huffman visible now)

plt.xlabel("DNA Size")
plt.ylabel("Time (seconds, log scale)")
plt.title("Compression Time Comparison (Logarithmic Scale)")
plt.legend()
plt.grid(True, which="both", ls="--", linewidth=0.5)
plt.tight_layout()

plt.show()
