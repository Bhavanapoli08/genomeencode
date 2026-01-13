import os
from Bio import SeqIO

def load_fasta(filename):
    # Correct path: genomeencode/data/
    base_path = os.path.join(os.path.dirname(__file__), "..", "data", filename)
    base_path = os.path.abspath(base_path)
    return SeqIO.read(base_path, "fasta").seq

if __name__ == "__main__":
    print("Small:", len(load_fasta("small.fasta")))
    print("Medium:", len(load_fasta("medium.fasta")))
    print("Large:", len(load_fasta("large.fasta")))
