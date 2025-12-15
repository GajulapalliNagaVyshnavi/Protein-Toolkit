
"""
Sequence Analysis Package
A comprehensive toolkit for protein sequence analysis and manipulation
"""

__version__ = "1.0.0"
__author__ = "JoVy"

# Import all functions to make them available at package level
from .sequence import (
    calculate_sequence_length,
    residue_to_single_letter,
    single_letter_to_residue,
    extract_chains_from_pdb,
    is_valid_sequence,
    clean_csv_sequences,
    extract_cdr_regions
)

from .conversions import (
    fasta_to_csv,
    csv_to_fasta,
    pdb_to_cif,
    cif_to_pdb,
)

from .embeddings import (  
    get_esm_embedding,
    get_protbert_embedding
)

# Define what gets imported with "from sequence_analysis import *"
__all__ = [
    'calculate_sequence_length',
    'residue_to_single_letter',
    'single_letter_to_residue',
    'extract_chains_from_pdb',
    'is_valid_sequence',
    'clean_csv_sequences',
    'extract_cdr_regions',
    'fasta_to_csv',
    'csv_to_fasta',
    'pdb_to_cif',
    'cif_to_pdb',
    'get_esm_embedding',
    'get_protbert_embedding'
]