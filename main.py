from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
from typing import Dict, List, Tuple
import numpy as np
import random

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

def calculate_probabilities(scores: Dict[str, float]) -> Dict[str, float]:
    """Calculate probabilities for each breed based on alignment scores using a numerically stable softmax."""
    max_score = max(scores.values())  # Find the maximum score

    # Scale scores to prevent large differences
    scaled_scores = {breed: score / 1000 for breed, score in scores.items()}  # Scale down scores
    max_scaled_score = max(scaled_scores.values())

    # Apply numerically stable softmax
    exp_scores = {breed: np.exp(score - max_scaled_score) for breed, score in scaled_scores.items()}
    total_score = sum(exp_scores.values())
    probabilities = {breed: exp_score / total_score for breed, exp_score in exp_scores.items()}

    # Sort probabilities by size and print the top 5
    sorted_probabilities = sorted(probabilities.items(), key=lambda x: x[1], reverse=True)
    print("Top 5 Probabilities:")
    for breed, prob in sorted_probabilities[:5]:
        print(f"{breed}: {prob:.4f}")

    return probabilities


def compute_p_value(best_score: float, reference_sequences: Dict[str, List[str]], query_length: int, num_simulations: int = 1000) -> float:
    """Compute p-value by comparing the best score to random sequence scores."""
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    random_scores = []

    for _ in range(num_simulations):
        random_query = ''.join(random.choices("ACGT", k=query_length))  # Generate random sequence
        for breed, sequences in reference_sequences.items():
            for ref_seq in sequences:
                alignment = aligner.align(ref_seq, random_query)
                random_scores.append(alignment[0].score)

    # Calculate p-value as the proportion of random scores >= best_score
    p_value = sum(1 for score in random_scores if score >= best_score) / len(random_scores)
    return p_value

def find_closest_match(reference_sequences: Dict[str, List[str]], query_file_path: str) -> Tuple[str, float, Dict[str, float], float]:
    """Align query sequences to the reference sequences and return the closest match, probabilities, and p-value."""
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    best_match = None
    best_score = float('-inf')
    scores = {}

    # Load query sequences
    query_sequences = load_sequences(query_file_path)

    # Align query sequences to reference sequences
    for query_seqs in query_sequences.values():
        for query_seq in query_seqs:
            for breed, sequences in reference_sequences.items():
                for ref_seq in sequences:
                    alignment = aligner.align(ref_seq, query_seq)
                    score = alignment[0].score
                    print(f"Aligning Mystery to {breed}, score: {score}")
                    
                    scores[breed] = max(scores.get(breed, float('-inf')), score)  # Track max score per breed

                    if score > best_score:
                        best_score = score
                        best_match = (breed, ref_seq)

    # Calculate probabilities
    probabilities = calculate_probabilities(scores)

    # Compute p-value
    query_length = len(query_sequences[list(query_sequences.keys())[0]][0])  # Length of the first query sequence
    p_value = compute_p_value(best_score, reference_sequences, query_length)

    return best_match, best_score, probabilities, p_value



if __name__ == "__main__":
    reference_file_path = "data/dog_breeds.fa"
    query_file_path = "data/mystery.fa"

    reference_sequences = load_sequences(reference_file_path)
    # print(reference_sequences.keys())

    best_match, best_score = find_closest_match(reference_sequences, query_file_path)
    print(f"Best match: {best_match[0]}")
    print(f"Score: {best_score}")
