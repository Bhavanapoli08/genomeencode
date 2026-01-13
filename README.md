# 🧬 GenomeCode

**GenomeCode** is a graphical and algorithmic application designed for efficient **DNA sequence compression and decompression** using a combination of **Burrows–Wheeler Transform (BWT)** and **Huffman Coding**.  
It aims to minimize storage requirements for large genomic datasets while preserving complete sequence integrity.

---

## 🚀 Features

- **Interactive Interface**
  - Input DNA sequences manually, from a file, or generate them randomly.
  - Perform compression and decompression operations in real time.

- **Efficient Compression**
  - **Burrows–Wheeler Transform (BWT):** Groups similar patterns together to enhance redundancy.
  - **Huffman Coding:** Assigns shorter binary codes to frequently occurring nucleotides for compact storage.

- **Accurate Decompression**
  - Restores the original DNA sequence perfectly through **Huffman Decoding** and **Inverse BWT**.

- **Visualization**
  - Demonstrates how BWT and Huffman algorithms work step-by-step in both compression and decompression.

---

## 🧠 Core Concepts

| Technique | Purpose | Description |
|------------|----------|-------------|
| **Burrows–Wheeler Transform (BWT)** | Pattern Grouping | Rearranges the DNA sequence to cluster similar symbols, improving compressibility. |
| **Huffman Coding** | Entropy Reduction | Assigns variable-length binary codes based on nucleotide frequency for lossless compression. |
| **Inverse BWT & Huffman Decoding** | Decompression | Reconstructs the original DNA sequence from compressed data. |

---

## 🧩 System Overview

**Compression Flow**
DNA Sequence → Burrows–Wheeler Transform → Huffman Coding → Compressed Output


**Decompression Flow**

Compressed Data → Huffman Decoding → Inverse BWT → Original DNA Sequence


---

## ⚙️ Installation

### Prerequisites
- Python 3.8 or higher  
- Required packages:
  ```bash
  pip install numpy matplotlib
Run the Application
python genomecode.py

🧪 Example Usage

Input Sequence
ACGTACGTGGGTTT

Compression
→ Burrows–Wheeler Transform
→ Huffman Encoding
→ Compressed Output (Binary Data)

Decompression
→ Huffman Decoding
→ Inverse BWT
→ Original Sequence: ACGTACGTGGGTTT

📊 Visualization

GenomeCode includes a visual module that:

Shows how BWT reorganizes DNA symbols.

Displays Huffman trees for binary encoding.

Animates the decompression pipeline for educational use.

🧬 Applications
Genomic data compression and archiving

Educational tool for algorithm visualization

Research on biological data compression methods

Integration into larger bioinformatics workflows

