# Dog-Breed DNA Identification Tool

This project is a Python-based tool for identifying the closest DNA sequence match from a database of dog breeds. It uses pairwise sequence alignment to compare query DNA sequences against a reference database.

## Features

- Parses DNA sequences from FASTA files.
- Identifies dog breeds from FASTA headers.
- Performs global pairwise alignment to find the closest match.
- Outputs the best matching breed and alignment score.

## Requirements

- Python 3.7 or higher
- Biopython library

Install the required dependencies using:

```sh
pip install biopython
```

## Usage

1. Place the reference DNA sequences in data/dog_breeds.fa.
2. Place the query DNA sequences in data/mystery.fa.
3. Run the script: 
   ```sh
   python main.py
   ```
4. The output will display the closest matching breed and alignment score

# Example Output

```
Aligning Mystery to Aidi, score: 1234.56
Aligning Mystery to Sloughi, score: 987.65
Best match: Aidi
Score: 1234.56
```