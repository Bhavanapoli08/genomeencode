import matplotlib.pyplot as plt
import pandas as pd

# -------------------------------
# METRICS (from your output)
# -------------------------------
data = {
    "Dataset": ["Small", "Medium", "Large"],
    "Original (Bytes)": [16569, 29903, 35938],
    "BWT (Bytes)": [16570, 29904, 35939],
    "Huffman (Bytes)": [7146, 12089, 13920],
    "BWT+Huffman (Bytes)": [7144, 13093, 15676],
    "BWT Ratio": [1.000, 0.0 + 1.000, 1.000],
    "Huffman Ratio": [0.431, 0.404, 0.387],
    "Combined Ratio": [0.431, 0.438, 0.436]
}

df = pd.DataFrame(data)

# -------------------------------
# CREATE TABLE IMAGE
# -------------------------------
fig, ax = plt.subplots(figsize=(10, 3))
ax.axis('off')

table = ax.table(
    cellText=df.values,
    colLabels=df.columns,
    cellLoc='center',
    loc='center'
)

table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1.2, 1.2)

plt.title("Compression Ratio Metrics Table", fontsize=14, pad=20)

plt.savefig("compression_ratio_table.png", dpi=300, bbox_inches='tight')
plt.show()

print("\nTable Image saved as: compression_ratio_table.png")
