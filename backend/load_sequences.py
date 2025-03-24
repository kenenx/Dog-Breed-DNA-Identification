from Bio import SeqIO
from Bio.Seq import Seq

def load_sequences(file_path):
    """Load sequences from a FASTA file."""
    sequences = []
    with open(file_path, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequences.append(str(record.seq))
    return sequences

if __name__ == "__main__":
    fasta_file_path = "data/dog_breeds.fa"  # Ensure this path is correct
    sequences = load_sequences(fasta_file_path)
    print(sequences)