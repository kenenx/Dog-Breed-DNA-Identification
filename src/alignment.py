from typing import Dict, List, Tuple
import numpy as np
import random
from Bio.Align import PairwiseAligner

def calculate_probabilities(scores: Dict[str, float]) -> Dict[str, float]:
    """
    Calculate probabilities for each breed based on alignment scores using a numerically stable softmax.

    Args:
        scores (Dict[str, float]): A dictionary where keys are breed names and values are alignment scores.

    Returns:
        Dict[str, float]: A dictionary where keys are breed names and values are probabilities.
    """
    max_score = max(scores.values())  # Find the maximum score to normalize
    normalized_scores = {breed: score - max_score for breed, score in scores.items()}  # Prevent overflow

    # Apply softmax to calculate probabilities
    exp_scores = {breed: np.exp(score) for breed, score in normalized_scores.items()}
    total_score = sum(exp_scores.values())
    probabilities = {breed: exp_score / total_score for breed, exp_score in exp_scores.items()}

    # Print the top 5 probabilities for debugging
    sorted_probabilities = sorted(probabilities.items(), key=lambda x: x[1], reverse=True)
    for breed, prob in sorted_probabilities[:5]:
        print(f"{breed}: {prob:.4f}")

    return probabilities

def find_closest_match(reference_sequences: Dict[str, List[str]], query_sequences: Dict[str, List[str]]) -> Tuple[str, float, Dict[str, float], Dict[str, float]]:
    """
    Align query sequences to the reference sequences and return the closest match, probabilities, and scores.

    Args:
        reference_sequences (Dict[str, List[str]]): A dictionary where keys are breed names and values are lists of DNA sequences.
        query_sequences (Dict[str, List[str]]): A dictionary where keys are query names and values are lists of DNA sequences.

    Returns:
        Tuple[str, float, Dict[str, float], Dict[str, float]]:
            - The best matching breed and its reference sequence.
            - The best alignment score.
            - A dictionary of probabilities for each breed.
            - A dictionary of alignment scores for each breed.
    """
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    scores = {}

    def align_query_to_reference(query_seq: str):
        """
        Align a single query sequence to all reference sequences.

        Args:
            query_seq (str): The query DNA sequence.

        Returns:
            Tuple[str, float, Dict[str, float]]:
                - The best matching breed and its reference sequence.
                - The best alignment score for this query.
                - A dictionary of alignment scores for all breeds.
        """
        local_best_score = float('-inf')
        local_best_match = None
        local_scores = {}

        for breed, sequences in reference_sequences.items():
            for ref_seq in sequences:
                alignment = aligner.align(ref_seq, query_seq)
                score = alignment[0].score
                print(f"Aligning Mystery to {breed}, score: {score}")

                local_scores[breed] = max(local_scores.get(breed, float('-inf')), score)
                if score > local_best_score:
                    local_best_score = score
                    local_best_match = (breed, ref_seq)

        return local_best_match, local_best_score, local_scores

    # Sequentially align all query sequences
    best_match = None
    best_score = float('-inf')
    for query_seqs in query_sequences.values():
        for query_seq in query_seqs:
            local_best_match, local_best_score, local_scores = align_query_to_reference(query_seq)
            scores.update(local_scores)
            if local_best_score > best_score:
                best_score = local_best_score
                best_match = local_best_match

    # Calculate probabilities based on scores
    probabilities = calculate_probabilities(scores)

    return best_match, best_score, probabilities, scores