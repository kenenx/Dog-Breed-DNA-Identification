from utils import load_sequences
from alignment import find_closest_match
from tabulate import tabulate 

if __name__ == "__main__":
    # File paths for reference and query DNA sequences
    reference_file_path = "data/dog_breeds.fa"
    query_file_path = "data/mystery.fa"

    # Load reference and query sequences from FASTA files
    reference_sequences = load_sequences(reference_file_path)
    query_sequences = load_sequences(query_file_path)

    # Perform alignment to find the closest match
    best_match, best_score, probabilities, scores = find_closest_match(reference_sequences, query_sequences)

    # Prepare data for the output table
    sorted_probabilities = sorted(probabilities.items(), key=lambda x: x[1], reverse=True)
    table_data = []
    for breed, prob in sorted_probabilities[:6]:  # Include top 6 entries
        table_data.append([breed, probabilities[breed], scores[breed]])

    # Print the results in a formatted table
    print(tabulate(table_data, headers=["Breed", "Probability", "Alignment Score"], tablefmt="pretty"))