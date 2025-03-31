# Dog-Breed DNA Identification Tool

This project is a Python-based tool for identifying the closest DNA sequence match from a database of dog breeds. It uses pairwise sequence alignment to compare query DNA sequences against a reference database.

## Features

- Parses DNA sequences from FASTA files.
- Identifies dog breeds from FASTA headers.
- Performs global pairwise alignment to find the closest match.
- Outputs the top 6 matching breeds, along with their alignment scores and probabilities.
- Provides a web-based interface to view results dynamically.

## Requirements

- Python 3.7 or higher

Install the required dependencies using:

```sh
pip install -r requirements.txt
```

## Usage

### Command-Line Interface
1. Place the reference DNA sequences in `data/dog_breeds.fa`.
2. Place the query DNA sequences in `data/mystery.fa`.
3. Run the script: 
   ```sh
   python src/main.py
   ```
4. The output will display the top 6 matching breeds in a table format.

### Web Interface
1. Place the reference DNA sequences in `data/dog_breeds.fa`.
2. Place the query DNA sequences in `data/mystery.fa`.
3. Run the Flask web server:
   ```sh
   python src/app.py
   ```
4. Open your browser and navigate to `http://127.0.0.1:5000/`.
5. Enter the file paths for the reference and query DNA sequences (default paths are pre-filled).
6. Click "Submit" to start the alignment process.
7. A loading page will display the progress of the alignment, including console messages such as alignment scores.
8. Once the alignment is complete, the results page will display the top 6 matching breeds, their probabilities, and alignment scores.

## Example Output (Web Interface)

### Loading Page
```
Processing...
Loading sequences...
Performing alignment...
Aligning Mystery to Beagle, score: 100
Aligning Mystery to Golden Retriever, score: 95
Alignment complete.
```

### Results Page
```
+-------------------+-------------+------------------+
| Breed             | Probability | Alignment Score  |
+-------------------+-------------+------------------+
| Beagle            | 0.8389      | 100              |
| Golden Retriever  | 0.1611      | 95               |
+-------------------+-------------+------------------+
```

## Tests

Run the tests using the following command:

```sh
python -m unittest discover -s tests
```