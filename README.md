# Dog-Breed DNA Identification Tool

This project is a Python-based tool for identifying the closest DNA sequence match from a database of dog breeds. It uses pairwise sequence alignment to compare query DNA sequences against a reference database.

## Features

- Parses DNA sequences from FASTA files.
- Identifies dog breeds from FASTA headers.
- Performs global pairwise alignment to find the closest match.
- Outputs the top 6 matching breeds, along with their alignment scores and probabilities in a table format.

## Requirements

- Python 3.7 or higher

Install the required dependencies using:

```sh
pip install -r requirements.txt
```

## Usage

1. Place the reference DNA sequences in `data/dog_breeds.fa`.
2. Place the query DNA sequences in `data/mystery.fa`.
3. Run the script: 
   ```sh
   python src/main.py
   ```
4. The output will display the top 6 matching breeds in a table format.

## Example Output

```
+-------------------------------------------------+-----------------------+------------------+
| Breed                                           | Probability           | Alignment Score  |
+-------------------------------------------------+-----------------------+------------------+
| English Springer Spaniel                        | 0.8389                | 16711.0          |
| Portuguese Warren dog, small size, weired hair | 0.1135                | 16699.0          |
| Mixed breed                                    | 0.0154                | 16689.0          |
| Sloughi                                        | 0.0154                | 16689.0          |
| Portuguese Sheepdog                            | 0.0057                | 16689.0          |
| Portuguese Warren dog, medium size, weired hair| 0.0021                | 16689.0          |
+-------------------------------------------------+-----------------------+------------------+
```

## Tests

Run the tests using the following command:

```sh
python -m unittest discover -s tests
```