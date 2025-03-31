from typing import Dict, List
from Bio import SeqIO

def get_breed(header: str) -> str:
    """Get dog breed from FASTA header."""
    parts = header.split("[breed=")
    if len(parts) > 1:
        return parts[1].split("]")[0]
    return "Unknown"

def load_sequences(file_path: str) -> Dict[str, List[str]]:
    """Load sequences and store them in a dictionary by breed."""
    sequences = {}
    with open(file_path, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            breed = get_breed(record.description)
            sequences.setdefault(breed, []).append(str(record.seq))
    return sequences