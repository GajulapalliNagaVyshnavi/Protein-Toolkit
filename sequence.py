"""
Sequence Analysis Package
A comprehensive toolkit for protein sequence analysis and manipulation
"""

import pandas as pd
from typing import Optional, Dict, List, Tuple
import re
from pathlib import Path


# Amino acid conversion dictionaries
THREE_TO_ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

ONE_TO_THREE = {v: k for k, v in THREE_TO_ONE.items()}

# Valid amino acids
VALID_AA = set('ACDEFGHIKLMNPQRSTVWY')


def calculate_sequence_length(sequence: str) -> int:
    """
    Calculate the length of a protein sequence.
    
    Args:
        sequence: Protein sequence string
        
    Returns:
        Length of the sequence
        
    Example:
        >>> calculate_sequence_length("ACDEFGHIKLM")
        11
    """
    return len(sequence.strip())


def residue_to_single_letter(residue: str) -> str:
    """
    Convert three-letter amino acid code to single-letter code.
    
    Args:
        residue: Three-letter amino acid code (e.g., 'HIS', 'ALA')
        
    Returns:
        Single-letter amino acid code (e.g., 'H', 'A')
        
    Raises:
        ValueError: If residue code is not recognized
        
    Example:
        >>> residue_to_single_letter("HIS")
        'H'
    """
    residue_upper = residue.strip().upper()
    if residue_upper not in THREE_TO_ONE:
        raise ValueError(f"Unknown residue code: {residue}")
    return THREE_TO_ONE[residue_upper]


def single_letter_to_residue(letter: str) -> str:
    """
    Convert single-letter amino acid code to three-letter code.
    
    Args:
        letter: Single-letter amino acid code (e.g., 'H', 'A')
        
    Returns:
        Three-letter amino acid code (e.g., 'HIS', 'ALA')
        
    Raises:
        ValueError: If letter code is not recognized
        
    Example:
        >>> single_letter_to_residue("H")
        'HIS'
    """
    letter_upper = letter.strip().upper()
    if letter_upper not in ONE_TO_THREE:
        raise ValueError(f"Unknown amino acid letter: {letter}")
    return ONE_TO_THREE[letter_upper]


def extract_chains_from_pdb(pdb_path: str) -> Dict[str, str]:
    """
    Extract all chains from a PDB file.
    
    Args:
        pdb_path: Path to PDB file
        
    Returns:
        Dictionary with chain IDs as keys and sequences as values
        
    Example:
        >>> chains = extract_chains_from_pdb("protein.pdb")
        >>> for chain_id, sequence in chains.items():
        ...     print(f"Chain {chain_id}: {sequence}")
    """
    chains = {}
    current_chain = None
    residues = {}
    
    try:
        with open(pdb_path, 'r') as f:
            for line in f:
                if line.startswith('ATOM'):
                    chain_id = line[21].strip()
                    if not chain_id:
                        chain_id = 'A'
                    
                    residue_name = line[17:20].strip()
                    residue_num = line[22:26].strip()
                    
                    if chain_id not in residues:
                        residues[chain_id] = {}
                    
                    if residue_num not in residues[chain_id]:
                        try:
                            aa_code = THREE_TO_ONE.get(residue_name, 'X')
                            residues[chain_id][residue_num] = aa_code
                        except:
                            residues[chain_id][residue_num] = 'X'
        
        for chain_id in residues:
            sorted_residues = sorted(residues[chain_id].items(), 
                                   key=lambda x: int(re.sub('[^0-9]', '', x[0]) or '0'))
            chains[chain_id] = ''.join([res[1] for res in sorted_residues])
        
        print("Extracted Chains:")
        for chain_id, sequence in chains.items():
            print(f"Chain {chain_id}: {sequence[:50]}{'...' if len(sequence) > 50 else ''} (Length: {len(sequence)})")
        
        return chains
        
    except FileNotFoundError:
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")
    except Exception as e:
        raise Exception(f"Error reading PDB file: {str(e)}")


def is_valid_sequence(sequence: str, allow_gaps: bool = False) -> bool:
    """
    Check if a protein sequence is valid.
    
    Args:
        sequence: Protein sequence string
        allow_gaps: Whether to allow gap characters ('-')
        
    Returns:
        True if sequence is valid, False otherwise
        
    Example:
        >>> is_valid_sequence("ACDEFGHIKLM")
        True
        >>> is_valid_sequence("ACDEFGHIKLMXYZ")
        False
    """
    if not sequence or not sequence.strip():
        return False
    
    sequence_upper = sequence.strip().upper()
    valid_chars = VALID_AA.copy()
    
    if allow_gaps:
        valid_chars.add('-')
    
    return all(char in valid_chars for char in sequence_upper)


def clean_csv_sequences(csv_path: str, 
                       target_column: str, 
                       output_path: Optional[str] = None,
                       allow_gaps: bool = False) -> pd.DataFrame:
    """
    Remove invalid sequences from a CSV file.
    
    Args:
        csv_path: Path to input CSV file
        target_column: Name of the column containing sequences
        output_path: Path to save cleaned CSV (optional)
        allow_gaps: Whether to allow gap characters in sequences
        
    Returns:
        DataFrame with invalid sequences removed
        
    Example:
        >>> df = clean_csv_sequences("sequences.csv", "sequence", "cleaned.csv")
        >>> print(f"Kept {len(df)} valid sequences")
    """
    try:
        df = pd.read_csv(csv_path)
    except FileNotFoundError:
        raise FileNotFoundError(f"CSV file not found: {csv_path}")
    except Exception as e:
        raise Exception(f"Error reading CSV file: {str(e)}")
    
    if target_column not in df.columns:
        raise ValueError(f"Column '{target_column}' not found in CSV. Available columns: {list(df.columns)}")
    
    initial_count = len(df)
    df_cleaned = df[df[target_column].apply(lambda x: is_valid_sequence(str(x), allow_gaps))].copy()
    removed_count = initial_count - len(df_cleaned)
    
    print(f"Total sequences: {initial_count}")
    print(f"Valid sequences: {len(df_cleaned)}")
    print(f"Invalid sequences removed: {removed_count}")
    
    if output_path:
        df_cleaned.to_csv(output_path, index=False)
        print(f"Cleaned CSV saved to: {output_path}")
    
    return df_cleaned



def extract_cdr_regions(sequence: str, scheme: str = "kabat") -> Dict[str, str]:
    """
    Extract CDR regions from an antibody sequence.
    Uses simplified position-based extraction (Kabat numbering approximation).
    
    Args:
        sequence: Antibody sequence string
        scheme: Numbering scheme ('kabat', 'chothia', 'imgt') - default 'kabat'
        
    Returns:
        Dictionary containing CDR1, CDR2, CDR3 sequences
        
    Note:
        This is a simplified implementation. For production use,
        consider using specialized tools like ANARCI or AbRSA.
        
    Example:
        >>> cdrs = extract_cdr_regions("EVQLVESGGGLVQPGG...")
        >>> print(f"CDR1: {cdrs['CDR1']}")
        >>> print(f"CDR2: {cdrs['CDR2']}")
        >>> print(f"CDR3: {cdrs['CDR3']}")
    """
    seq_len = len(sequence)
    cdrs = {}
    
    if scheme.lower() == "kabat":
        if seq_len >= 35:
            cdrs['CDR1'] = sequence[26:35]
        if seq_len >= 66:
            cdrs['CDR2'] = sequence[50:66]
        if seq_len >= 102:
            cdrs['CDR3'] = sequence[95:102]
    elif scheme.lower() == "chothia":
        if seq_len >= 34:
            cdrs['CDR1'] = sequence[26:34]
        if seq_len >= 56:
            cdrs['CDR2'] = sequence[50:56]
        if seq_len >= 102:
            cdrs['CDR3'] = sequence[95:102]
    elif scheme.lower() == "imgt":
        if seq_len >= 39:
            cdrs['CDR1'] = sequence[27:39]
        if seq_len >= 63:
            cdrs['CDR2'] = sequence[56:63]
        if seq_len >= 118:
            cdrs['CDR3'] = sequence[105:118]
    else:
        raise ValueError(f"Unknown scheme: {scheme}. Use 'kabat', 'chothia', or 'imgt'")
    
    if not cdrs:
        print(f"Warning: Sequence too short ({seq_len} aa) to extract CDRs using {scheme} scheme")
    else:
        print(f"CDR Regions ({scheme.upper()} scheme):")
        for cdr_name, cdr_seq in cdrs.items():
            print(f"{cdr_name}: {cdr_seq}")
    
    return cdrs


# # Example usage
# if __name__ == "__main__":
#     print("=== Sequence Analysis Package ===\n")
    
#     # 1. Sequence length
#     seq = "MFSKLAHLQRFAVLSRGVHSSVASATSVATKKTVQGPPTSDDIFEREYKYGAHNYHPLPVALERGKGIYLWDVEGRKYFDFLSSYSAVNQGHCHPKIVNALKSQVDKLTLTSRAFYNNVLGEYEEYITKLFNYHKVLPMNTGVEAGETACKLARKWGYTVKGIQKYKAKIVFAAGNFWGRTLSAISSSTDPTSYDGFGPFMPGFDIIPYNDLPALERALQDPNVAAFMVEPIQGEAGVVVPDPGYLMGVRELCTRHQVLFIADEIQTGLARTGRWLAVDYENVRPDIVLLGKALSGGLYPVSAVLCDDDIMLTIKPGEHGSTYGGNPLGCRVAIAALEVLEEENLAENADKLGIILRNELMKLPSDVVTAVRGKGLLNAIVIKETKDWDAWKVCLRLRDNGLLAKPTHGDIIRFAPPLVIKEDELRESIEIINKTILSF"
#     print(f"1. Sequence Length: {calculate_sequence_length(seq)}\n")
    
#     # 2. Residue to single letter
#     print(f"2. HIS -> {residue_to_single_letter('HIS')}\n")
    
#     # 3. Single letter to residue
#     print(f"3. H -> {single_letter_to_residue('H')}\n")
    
#     # 4. Extract chains (example - requires PDB file)
#     chains = extract_chains_from_pdb("/home/boltzmann6/hackathon_dec/sample_inputs/WT_aminotransferase.pdb")
    
#     # 5. Valid sequence check
#     print(f"4. Is valid? {is_valid_sequence('MFSKLAHLQRFAVLSRGVHSSVASATSVATKKTVQGPPTSDDIFEREYKYGAHNYHPLPVALERGKGIYLWDVEGRKYFDFLSSYSAVNQGHCHPKIVNALKSQVDKLTLTSRAFYNNVLGEYEEYITKLFNYHKVLPMNTGVEAGETACKLARKWGYTVKGIQKYKAKIVFAAGNFWGRTLSAISSSTDPTSYDGFGPFMPGFDIIPYNDLPALERALQDPNVAAFMVEPIQGEAGVVVPDPGYLMGVRELCTRHQVLFIADEIQTGLARTGRWLAVDYENVRPDIVLLGKALSGGLYPVSAVLCDDDIMLTIKPGEHGSTYGGNPLGCRVAIAALEVLEEENLAENADKLGIILRNELMKLPSDVVTAVRGKGLLNAIVIKETKDWDAWKVCLRLRDNGLLAKPTHGDIIRFAPPLVIKEDELRESIEIINKTILSF')}")
#     # print(f"   Is 'ACDEFXYZ' valid? {is_valid_sequence('ACDEFXYZ')}\n")
    
#     # 6. Clean CSV (example - requires CSV file)
#     df = clean_csv_sequences("/home/boltzmann6/hackathon_dec/sample_inputs/sample.csv", "Sequence", "/home/boltzmann6/hackathon_dec/sample_inputs/cleaned.csv")
    
#     # 7. CDR extraction
#     antibody_seq = "MFSKLAHLQRFAVLSRGVHSSVASATSVATKKTVQGPPTSDDIFEREYKYGAHNYHPLPVALERGKGIYLWDVEGRKYFDFLSSYSAVNQGHCHPKIVNALKSQVDKLTLTSRAFYNNVLGEYEEYITKLFNYHKVLPMNTGVEAGETACKLARKWGYTVKGIQKYKAKIVFAAGNFWGRTLSAISSSTDPTSYDGFGPFMPGFDIIPYNDLPALERALQDPNVAAFMVEPIQGEAGVVVPDPGYLMGVRELCTRHQVLFIADEIQTGLARTGRWLAVDYENVRPDIVLLGKALSGGLYPVSAVLCDDDIMLTIKPGEHGSTYGGNPLGCRVAIAALEVLEEENLAENADKLGIILRNELMKLPSDVVTAVRGKGLLNAIVIKETKDWDAWKVCLRLRDNGLLAKPTHGDIIRFAPPLVIKEDELRESIEIINKTILSF"
#     print("5. CDR Extraction:")
#     cdrs = extract_cdr_regions(antibody_seq)
#     print()
    