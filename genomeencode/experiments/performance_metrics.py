import time
import os
import matplotlib.pyplot as plt

from genomeencode.burros_wheeler import BurrosWheeler
from genomeencode.huffman import HuffmanTree


# -------------------------------------------------------
# Helper functions
# -------------------------------------------------------
def measure_time(fn, *args):
    """Measure execution time of any function."""
    start = time.time()
    result = fn(*args)
    end = time.time()
    return result, end - start


def compression_ratio(original, compressed):
    """Calculate compression ratio = original_size / compressed_size."""
    return len(original) / max(len(compressed), 1)


# -------------------------------------------------------
# Test DNA sequences
# -------------------------------------------------------
short_dna  = "ACGTACGTACGTACGT"
medium_dna = "ACGTTGCA" * 2000        # ~16,000 bases
large_dna  = "ACGT" * 20000           # ~80,000 bases

datasets = {
    "Short DNA": short_dna,
    "Medium DNA": medium_dna,
    "Large DNA": large_dna
}


# -------------------------------------------------------
# Storage for plotting
# -------------------------------------------------------
bwt_times = []
huff_times = []
combined_times = []

decode_times = []
ratios = []
labels = []


# -------------------------------------------------------
# Actual Experiment Loop
# -------------------------------------------------------
for label, dna in datasets.items():

    print(f"\n--- Testing on {label} ({len(dna)} bases) ---")

    labels.append(label)

    # ------------------ BWT Only ------------------------
    (_, t_bwt) = measure_time(lambda seq: BurrosWheeler.bwt_advanced(seq), dna)
    bwt_times.append(t_bwt)

    # ------------------ Huffman Only --------------------
    huff_tree = HuffmanTree(dna)
    huff_tree.get_codings(huff_tree.root)

    (_, t_huff) = measure_time(lambda: huff_tree.seq_to_binstr())
    huff_times.append(t_huff)

    # ------------------ Combined (BWT + Huffman) --------
    bwt_seq = BurrosWheeler.bwt_advanced(dna)
    ct = HuffmanTree(bwt_seq)
    ct.get_codings(ct.root)

    (_, t_comb) = measure_time(lambda: ct.seq_to_binstr())
    combined_times.append(t_comb)

    # ------------------ Decompression -------------------
    bin_data = ct.seq_to_binstr()
    _, t_decode = measure_time(lambda: HuffmanTree.binstr_to_seq(bin_data, ct.codes))
    decode_times.append(t_decode)

    # ------------------ Compression Ratio ---------------
    ratio = compression_ratio(dna, bin_data)
    ratios.append(ratio)


# -------------------------------------------------------
# Plot 1 — Compression Time
# -------------------------------------------------------
plt.figure(figsize=(10, 6))
plt.plot(labels, bwt_times, marker='o', label="BWT")
plt.plot(labels, huff_times, marker='o', label="Huffman")
plt.plot(labels, combined_times, marker='o', label="BWT + Huffman")
plt.title("Compression Time Comparison")
plt.xlabel("DNA Size")
plt.ylabel("Time (seconds)")
plt.grid(True)
plt.legend()
plt.show()


# -------------------------------------------------------
# Plot 2 — Decompression Time
# -------------------------------------------------------
plt.figure(figsize=(10, 6))
plt.plot(labels, decode_times, marker='o', label="Full Decompression")
plt.title("Decompression Time Comparison")
plt.xlabel("DNA Size")
plt.ylabel("Time (seconds)")
plt.grid(True)
plt.legend()
plt.show()


# -------------------------------------------------------
# Plot 3 — Compression Ratio
# -------------------------------------------------------
plt.figure(figsize=(10, 6))
plt.bar(labels, ratios, color="green")
plt.title("Compression Ratio (Higher is Better)")
plt.xlabel("DNA Size")
plt.ylabel("Compression Ratio")
plt.grid(axis='y')
plt.show()

print("\nDONE — All metrics calculated and plots generated.")
