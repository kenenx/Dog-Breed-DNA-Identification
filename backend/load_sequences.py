from Bio import SeqIO
from Bio.Seq import Seq

def get_breed(header):
    """Get dog breed from FASTA header."""
    parts = header.split("[breed=")
    if len(parts) > 1:
        return parts[1].split("]")[0]
    return "Unknown"

def load_sequences(file_path):
    """Load sequences and store them in a dictionary by breed."""
    breed_sequences = {}
    with open(file_path, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            breed = get_breed(record.description)
            breed_sequences.setdefault(breed, []).append(str(record.seq))
    return breed_sequences

if __name__ == "__main__":
    fasta_file_path = "data/dog_breeds.fa"  # Ensure this path is correct
    sequences = load_sequences(fasta_file_path)
    print(sequences.keys())