# Protein Engineering Toolkit

`protein_eng` is a Python package for protein sequence analysis, conversions, and embedding generation.

## Installation

```bash
pip install .
```

## Available Functions

### Sequence Analysis
- `calculate_sequence_length(sequence)`: Returns the length of a protein sequence.
- `is_valid_sequence(sequence)`: Checks if a sequence contains only valid amino acids.
- `residue_to_single_letter(residue)`: Converts 3-letter amino acid code (e.g., 'ALA') to 1-letter code ('A').
- `single_letter_to_residue(letter)`: Converts 1-letter amino acid code (e.g., 'A') to 3-letter code ('ALA').
- `extract_cdr_regions(sequence)`: Extracts CDR1, CDR2, and CDR3 regions from an antibody sequence.

### File Conversions
- `fasta_to_csv(fasta_path, csv_path)`: Converts a FASTA file to a CSV file.
- `csv_to_fasta(csv_path, target_column, fasta_path)`: Converts a CSV file to a FASTA file.
- `pdb_to_cif(pdb_path, cif_path)`: Converts a PDB structure file to MMCIF format.
- `cif_to_pdb(cif_path, pdb_path)`: Converts an MMCIF structure file to PDB format.

### Structure Analysis
- `extract_chains_from_pdb(pdb_path)`: Extracts sequences for all chains in a PDB file.

### Data Cleaning
- `clean_csv_sequences(csv_path, target_column)`: Removes rows with invalid protein sequences from a CSV file.

### Embeddings
- `get_esm_embedding(sequence)`: Generates embeddings using the ESM-2 model.
- `get_protbert_embedding(sequence)`: Generates embeddings using the ProtBERT model.

## Usage Example

```python
import protein_eng

# Calculate sequence length
length = protein_eng.calculate_sequence_length("ACDEF")

# Convert FASTA to CSV
protein_eng.fasta_to_csv("input.fasta", "output.csv")
```
