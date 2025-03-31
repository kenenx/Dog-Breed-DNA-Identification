import unittest
from unittest.mock import patch, mock_open
from src.utils import get_breed, load_sequences
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

class TestUtils(unittest.TestCase):
    def test_get_breed(self):
        try:
            self.assertEqual(get_breed("dog1 [breed=Beagle]"), "Beagle")
            self.assertEqual(get_breed("dog2 [breed=Golden Retriever]"), "Golden Retriever")
            self.assertEqual(get_breed("dog3"), "Unknown")
            print("test_get_breed: PASS")
        except AssertionError as e:
            print("test_get_breed: FAIL")
            raise e

    @patch("builtins.open", new_callable=mock_open, read_data=">dog1 [breed=Beagle]\nATCG\n>dog2 [breed=Golden Retriever]\nGCTA\n")
    @patch("Bio.SeqIO.parse")
    def test_load_sequences(self, mock_seqio_parse, mock_file):
        try:
            mock_seqio_parse.return_value = [
                SeqRecord(Seq("ATCG"), id="dog1", description="dog1 [breed=Beagle]"),
                SeqRecord(Seq("GCTA"), id="dog2", description="dog2 [breed=Golden Retriever]")
            ]

            sequences = load_sequences("mock_file_path")
            self.assertIn("Beagle", sequences)
            self.assertIn("Golden Retriever", sequences)
            self.assertEqual(sequences["Beagle"], ["ATCG"])
            self.assertEqual(sequences["Golden Retriever"], ["GCTA"])
            print("test_load_sequences: PASS")
        except AssertionError as e:
            print("test_load_sequences: FAIL")
            raise e

if __name__ == "__main__":
    unittest.main()
