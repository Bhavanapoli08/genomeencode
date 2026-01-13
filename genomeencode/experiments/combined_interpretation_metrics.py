# combined_interpretation_metrics.py
# Generates a combined table + combined graph
# (Compression + Decompression Metrics)

import matplotlib.pyplot as plt
import numpy as np

# ============================
# 1. HARDCODED METRICS (YOUR RESULTS)
# ============================

datasets = ["Small", "Medium", "Large"]

orig_sizes = [16569, 29903, 35938]

# --- Compression Speed (Bytes/sec) ---
bwt_comp_speed   = [399996.68, 275684.85, 288851.05]
huff_comp_speed  = [5574797.29, 5831424.24, 5593754.30]
comb_comp_speed  = [594867.73, 399303.01, 274365.25]

# --- Decompression Speed (Bytes/sec) ---
bwt_decomp_speed  = [2396393895.72, 2986244583.62, 3207125471.32]
huff_decomp_speed = [2147638.15, 2409233.23, 2461018.09]
comb_decomp_speed = [2003211.78, 2100699.65, 2040183.77]

# ============================
# 2. PRINT COMBINED TABLE
# ============================

print("\n==== COMBINED INTERPRETATION METRICS ====\n")
print("Dataset    Orig(Bytes)   BWT_Comp     Huff_Comp    Comb_Comp     BWT_Decomp      Huff_Decomp    Comb_Decomp")
print("--------------------------------------------------------------------------------------------------------------")

for i in range(len(datasets)):
    print(f"{datasets[i]:<10} {orig_sizes[i]:<13} "
          f"{bwt_comp_speed[i]:<12.2f} {huff_comp_speed[i]:<12.2f} {comb_comp_speed[i]:<13.2f} "
          f"{bwt_decomp_speed[i]:<15.2f} {huff_decomp_speed[i]:<14.2f} {comb_decomp_speed[i]:<14.2f}")

# ============================
# 3. CREATE A COMBINED GRAPH
# ============================

plt.figure(figsize=(14, 7))

x = np.arange(len(datasets))

plt.plot(x, bwt_comp_speed, marker="o", label="BWT Compression")
plt.plot(x, huff_comp_speed, marker="o", label="Huffman Compression")
plt.plot(x, comb_comp_speed, marker="o", label="Combined Compression")

plt.plot(x, bwt_decomp_speed, marker="s", label="BWT Decompression")
plt.plot(x, huff_decomp_speed, marker="s", label="Huffman Decompression")
plt.plot(x, comb_decomp_speed, marker="s", label="Combined Decompression")

plt.xticks(x, datasets)
plt.xlabel("Dataset Size")
plt.ylabel("Speed (Bytes/sec)")
plt.title("Combined Performance Comparison (Compression + Decompression)")
plt.yscale("log")  # LOG SCALE because decompression is extremely fast
plt.grid(True)
plt.legend()
plt.tight_layout()

plt.savefig("combined_interpretation_graph.png")
plt.show()

print("\nGraph saved as: combined_interpretation_graph.png\n")
