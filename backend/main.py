from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner

def get_breed(header):
    """Get dog breed from FASTA header."""
    parts = header.split("[breed=")
    if len(parts) > 1:
        return parts[1].split("]")[0]
    return "Unknown"

def load_sequences(file_path):
    """Load sequences and store them in a dictionary by breed."""
    sequences = {}
    with open(file_path, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            breed = get_breed(record.description)
            sequences.setdefault(breed, []).append(str(record.seq))
    return sequences

def find_closest_match(reference_sequences, query_file_path):
    """Align query sequences to the reference sequences and return the closest match and difference."""
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    best_match = None
    best_score = float('-inf')

    # Load query sequences
    query_sequences = []
    with open(query_file_path, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            query_sequences.append(str(record.seq))

    # Align query sequences to reference sequences
    for query_seq in query_sequences:
        for breed, sequences in reference_sequences.items():
            for ref_seq in sequences:
                alignment = aligner.align(ref_seq, query_seq)
                score = alignment[0].score
                print(f"Aligning Mystery to {breed}, score: {score}") # Function crashes if this debug line is removed as the dataset is too big

                if score > best_score:
                    best_score = score
                    best_match = (breed, ref_seq)


    return best_match, best_score

if __name__ == "__main__":
    reference_file_path = "data/dog_breeds.fa"
    query_file_path = "data/mystery.fa"

    reference_sequences = load_sequences(reference_file_path)
    # print(reference_sequences.keys())

    best_match, best_score = find_closest_match(reference_sequences, query_file_path)
    print(f"Best match: {best_match[0]}")
    print(f"Score: {best_score}")
