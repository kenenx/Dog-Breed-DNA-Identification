import unittest
from unittest.mock import patch, mock_open
from src.alignment import find_closest_match  # Corrected import path
from src.utils import load_sequences  # Corrected import path

class TestIntegration(unittest.TestCase):
    @patch("builtins.open", new_callable=mock_open, read_data=">dog1 [breed=Beagle]\nATCG\n>dog2 [breed=Golden Retriever]\nGCTA\n")
    @patch("Bio.Align.PairwiseAligner.align")
    def test_integration(self, mock_align, mock_file):
        """Test the integration of loading sequences and finding the closest match."""
        try:
            mock_align.return_value = [type("MockAlignment", (object,), {"score": 100})()]
            reference_sequences = load_sequences("data/dog_breeds.fa")
            query_sequences = load_sequences("data/mystery.fa")
            best_match, best_score, probabilities, scores = find_closest_match(reference_sequences, query_sequences)

            self.assertEqual(best_match[0], "Beagle")
            self.assertEqual(best_score, 100)
            self.assertIn("Beagle", probabilities)
            self.assertIn("Golden Retriever", probabilities)
            self.assertAlmostEqual(sum(probabilities.values()), 1.0, places=2)
            self.assertEqual(scores["Beagle"], 100)
            print("test_integration: PASS")
        except AssertionError as e:
            print("test_integration: FAIL")
            raise e

if __name__ == "__main__":
    unittest.main()
