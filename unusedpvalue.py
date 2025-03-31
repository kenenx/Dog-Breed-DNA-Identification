 # def compute_p_value(best_score: float, reference_sequences: Dict[str, List[str]], query_length: int, num_simulations: int = 1000, p_value_threshold: float = 0.01) -> float:
   """Compute p-value by comparing the best score to random sequence scores."""
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    threshold_count = 0  # Count of random scores >= best_score
    total_scores = 0  # Total number of scores processed

    def align_random_query(random_query: str) -> int:
        """Align a random query to reference sequences and count scores >= best_score."""
        local_count = 0
        for breed, sequences in reference_sequences.items():
            for ref_seq in sequences:
                alignment = aligner.align(ref_seq, random_query)
                score = alignment[0].score
                if score >= best_score:
                    local_count += 1
        return local_count
        futures = []
        for _ in range(num_simulations):
            random_query = ''.join(random.choices("ACGT", k=query_length))  # Generate random sequence
            futures.append(executor.submit(align_random_query, random_query))

        for future in futures:
            threshold_count += future.result()
            total_scores += 1

            # Early stopping if p-value is below the threshold
            if threshold_count / total_scores <= p_value_threshold:
                break

    # Calculate p-value as the proportion of random scores >= best_score
    p_value = threshold_count / total_scores
    return p_value