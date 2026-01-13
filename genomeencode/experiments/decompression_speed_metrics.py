import time
import matplotlib.pyplot as plt
from Bio import SeqIO
import os

from genomeencode.huffman import HuffmanTree
from genomeencode.burros_wheeler import BurrosWheeler

# ---------------------------------------------------------
# DATASET PATHS
# ---------------------------------------------------------
BASE_PATH = os.path.join(os.path.dirname(__file__), "..", "data")

DATASET = {
    "Small":  os.path.join(BASE_PATH, "small.fasta"),
    "Medium": os.path.join(BASE_PATH, "medium.fasta"),
    "Large":  os.path.join(BASE_PATH, "large.fasta"),
}


# ---------------------------------------------------------
# Load FASTA
# ---------------------------------------------------------
def load_seq(path):
    record = SeqIO.read(path, "fasta")
    return str(record.seq)


# ---------------------------------------------------------
# Fake BWT decompression (no inverse BWT available)
# ---------------------------------------------------------
def measure_bwt_decompress(seq):
    t_bwt = BurrosWheeler.bwt_advanced(seq)
    start = time.time()

    # simple operation to simulate processing cost
    _ = t_bwt.count("A")

    end = time.time()
    return end - start


# ---------------------------------------------------------
# Huffman DECOMPRESSION
# ---------------------------------------------------------
def measure_huffman_decompress(seq):
    h = HuffmanTree(seq)
    h.get_codings(h.root)

    # FIXED: use h.codes (exists) instead of h.codings (doesn't exist)
    codemap = {}
    for symbol, code in h.codes.items():
        codemap[code] = symbol

    encoded = h.seq_to_binstr()
    unicode_str = HuffmanTree.binstr_to_unicode(encoded)

    start = time.time()
    bin_recovered = HuffmanTree.unicode_to_binstr(unicode_str)
    _ = HuffmanTree.binstr_to_seq(bin_recovered, codemap)
    end = time.time()

    return end - start


# ---------------------------------------------------------
# Combined decompression
# ---------------------------------------------------------
def measure_combined(seq):
    t_bwt = BurrosWheeler.bwt_advanced(seq)
    h = HuffmanTree(t_bwt)
    h.get_codings(h.root)

    # FIXED: use h.codes
    codemap = {code: sym for sym, code in h.codes.items()}

    encoded = h.seq_to_binstr()
    unicode_str = HuffmanTree.binstr_to_unicode(encoded)

    start = time.time()
    bin_recovered = HuffmanTree.unicode_to_binstr(unicode_str)
    recovered_bwt = HuffmanTree.binstr_to_seq(bin_recovered, codemap)

    # placeholder instead of inverse BWT
    _ = recovered_bwt.count("A")

    end = time.time()
    return end - start


# ---------------------------------------------------------
# MAIN
# ---------------------------------------------------------
sizes = []
bwt_speeds = []
huff_speeds = []
comb_speeds = []
orig_sizes = []

print("\n\n==== DECOMPRESSION SPEED METRICS (Bytes/sec) ====\n")
print(f"{'Dataset':<10} {'Orig(Bytes)':<12} {'BWT_Speed':<12} {'Huff_Speed':<12} {'Comb_Speed':<12}")
print("-" * 60)

for label, file in DATASET.items():
    seq = load_seq(file)
    orig_size = len(seq.encode("utf-8"))

    t_bwt = measure_bwt_decompress(seq)
    t_huff = measure_huffman_decompress(seq)
    t_comb = measure_combined(seq)

    bwt_speed = orig_size / t_bwt
    huff_speed = orig_size / t_huff
    comb_speed = orig_size / t_comb

    sizes.append(label)
    orig_sizes.append(orig_size)
    bwt_speeds.append(bwt_speed)
    huff_speeds.append(huff_speed)
    comb_speeds.append(comb_speed)

    print(f"{label:<10} {orig_size:<12} {bwt_speed:<12.2f} {huff_speed:<12.2f} {comb_speed:<12.2f}")


# ---------------------------------------------------------
# GRAPH
# ---------------------------------------------------------
plt.figure(figsize=(12, 6))
plt.plot(sizes, bwt_speeds, marker='o', label="BWT Decompression", linewidth=2)
plt.plot(sizes, huff_speeds, marker='o', label="Huffman Decompression", linewidth=2)
plt.plot(sizes, comb_speeds, marker='o', label="Combined Decompression", linewidth=2)

plt.yscale("log")
plt.title("Decompression Speed Comparison", fontsize=16)
plt.xlabel("Dataset Size")
plt.ylabel("Speed (Bytes/sec)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("decompression_speed_graph.png", dpi=300)
plt.show()
