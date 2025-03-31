import unittest
from unittest.mock import patch
from src.alignment import calculate_probabilities, find_closest_match

class TestAlignment(unittest.TestCase):
    def test_calculate_probabilities(self):
        """Test the calculate_probabilities function."""
        try:
            scores = {"Beagle": 100, "Golden Retriever": 200, "Bulldog": 50}
            probabilities = calculate_probabilities(scores)
            self.assertIn("Beagle", probabilities)
            self.assertIn("Golden Retriever", probabilities)
            self.assertIn("Bulldog", probabilities)
            self.assertAlmostEqual(sum(probabilities.values()), 1.0, places=2)
            print("test_calculate_probabilities: PASS")
        except AssertionError as e:
            print("test_calculate_probabilities: FAIL")
            raise e

    @patch("Bio.Align.PairwiseAligner.align")
    def test_find_closest_match(self, mock_align):
        """Test the find_closest_match function."""
        try:
            mock_align.return_value = [type("MockAlignment", (object,), {"score": 100})()]
            reference_sequences = {"Beagle": ["ATCG"], "Golden Retriever": ["GCTA"]}
            query_sequences = {"Mystery": ["ATCG"]}
            best_match, best_score, probabilities, scores = find_closest_match(reference_sequences, query_sequences)
            self.assertEqual(best_match[0], "Beagle")
            self.assertEqual(best_score, 100)
            self.assertIn("Beagle", probabilities)
            self.assertIn("Golden Retriever", probabilities)
            self.assertAlmostEqual(sum(probabilities.values()), 1.0, places=2)
            self.assertEqual(scores["Beagle"], 100)
            print("test_find_closest_match: PASS")
        except AssertionError as e:
            print("test_find_closest_match: FAIL")
            raise e

if __name__ == "__main__":
    unittest.main()
