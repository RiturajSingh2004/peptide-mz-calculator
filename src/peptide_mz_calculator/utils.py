"""
Utilities module for peptide m/z calculator.

This module provides helper functions that don't fit into the other specialized modules.
"""

import numpy as np
from typing import Dict, List, Any, Tuple, Optional, Union
import re


def get_aa_composition(sequence: str) -> Dict[str, int]:
    """
    Calculate the amino acid composition of a peptide.
    
    Args:
        sequence: Amino acid sequence
        
    Returns:
        Dictionary of amino acids and their counts
    """
    composition = {}
    for aa in sequence:
        if aa in composition:
            composition[aa] += 1
        else:
            composition[aa] = 1
    
    return composition


def calculate_molecular_weight(composition: Dict[str, int]) -> float:
    """
    Calculate molecular weight from amino acid composition.
    
    Args:
        composition: Dictionary of amino acids and their counts
        
    Returns:
        Molecular weight in Daltons
    """
    # Average masses of amino acids (Da)
    aa_masses = {
        'A': 71.0788,
        'C': 103.1388,
        'D': 115.0886,
        'E': 129.1155,
        'F': 147.1766,
        'G': 57.0519,
        'H': 137.1411,
        'I': 113.1594,
        'K': 128.1741,
        'L': 113.1594,
        'M': 131.1926,
        'N': 114.1038,
        'P': 97.1167,
        'Q': 128.1307,
        'R': 156.1875,
        'S': 87.0782,
        'T': 101.1051,
        'V': 99.1326,
        'W': 186.2132,
        'Y': 163.1760
    }
    
    mw = 0.0
    for aa, count in composition.items():
        if aa in aa_masses:
            mw += aa_masses[aa] * count
    
    # Add mass of water (terminal OH and H)
    mw += 18.01528
    
    return mw


def validate_sequence(sequence: str) -> Tuple[bool, List[str]]:
    """
    Validate if sequence contains only valid amino acid codes.
    
    Args:
        sequence: Amino acid sequence
        
    Returns:
        Tuple of (is_valid, invalid_characters)
    """
    valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
    invalid_chars = [aa for aa in sequence if aa not in valid_aa]
    
    is_valid = len(invalid_chars) == 0
    
    return is_valid, invalid_chars


def format_formula(formula: str) -> str:
    """
    Format a chemical formula with HTML subscripts.
    
    Args:
        formula: Chemical formula string
        
    Returns:
        HTML formatted formula
    """
    # Replace digits with subscript HTML
    formatted = ""
    i = 0
    while i < len(formula):
        if formula[i].isalpha():
            formatted += formula[i]
            i += 1
        else:
            # Collect all consecutive digits
            digits = ""
            while i < len(formula) and formula[i].isdigit():
                digits += formula[i]
                i += 1
            formatted += f"<sub>{digits}</sub>"
    
    return formatted


def parse_modification_string(mod_string: str) -> Dict[str, Any]:
    """
    Parse a modification string into a structured format.
    
    Args:
        mod_string: Modification string (e.g., "Phospho(S12)")
        
    Returns:
        Dictionary with modification details
    """
    # Pattern: [ModName]([Residue][Position])
    pattern = r"([A-Za-z]+)\(([A-Z])(\d+|N-term|C-term)\)"
    match = re.match(pattern, mod_string)
    
    if match:
        mod_name, residue, position = match.groups()
        
        return {
            "name": mod_name,
            "residue": residue,
            "position": position
        }
    else:
        # Simple modification without position
        return {
            "name": mod_string,
            "residue": None,
            "position": None
        }


def get_aa_molecular_formula(aa: str) -> Dict[str, int]:
    """
    Get the molecular formula of an amino acid.
    
    Args:
        aa: Single-letter amino acid code
        
    Returns:
        Dictionary with elements as keys and counts as values
    """
    # Molecular formulas of amino acids
    aa_formulas = {
        'A': {'C': 3, 'H': 7, 'N': 1, 'O': 2},
        'C': {'C': 3, 'H': 7, 'N': 1, 'O': 2, 'S': 1},
        'D': {'C': 4, 'H': 7, 'N': 1, 'O': 4},
        'E': {'C': 5, 'H': 9, 'N': 1, 'O': 4},
        'F': {'C': 9, 'H': 11, 'N': 1, 'O': 2},
        'G': {'C': 2, 'H': 5, 'N': 1, 'O': 2},
        'H': {'C': 6, 'H': 9, 'N': 3, 'O': 2},
        'I': {'C': 6, 'H': 13, 'N': 1, 'O': 2},
        'K': {'C': 6, 'H': 14, 'N': 2, 'O': 2},
        'L': {'C': 6, 'H': 13, 'N': 1, 'O': 2},
        'M': {'C': 5, 'H': 11, 'N': 1, 'O': 2, 'S': 1},
        'N': {'C': 4, 'H': 8, 'N': 2, 'O': 3},
        'P': {'C': 5, 'H': 9, 'N': 1, 'O': 2},
        'Q': {'C': 5, 'H': 10, 'N': 2, 'O': 3},
        'R': {'C': 6, 'H': 14, 'N': 4, 'O': 2},
        'S': {'C': 3, 'H': 7, 'N': 1, 'O': 3},
        'T': {'C': 4, 'H': 9, 'N': 1, 'O': 3},
        'V': {'C': 5, 'H': 11, 'N': 1, 'O': 2},
        'W': {'C': 11, 'H': 12, 'N': 2, 'O': 2},
        'Y': {'C': 9, 'H': 11, 'N': 1, 'O': 3}
    }
    
    if aa in aa_formulas:
        return aa_formulas[aa]
    else:
        return {}


def calculate_mass_error(
    theoretical_mz: float, 
    experimental_mz: float, 
    unit: str = "ppm"
) -> float:
    """
    Calculate mass error between theoretical and experimental m/z values.
    
    Args:
        theoretical_mz: Theoretical m/z value
        experimental_mz: Experimental m/z value
        unit: Unit of error ("ppm" or "da")
        
    Returns:
        Mass error in the specified unit
    """
    if unit.lower() == "ppm":
        # Parts per million
        error = (experimental_mz - theoretical_mz) / theoretical_mz * 1e6
    else:
        # Absolute error in Daltons
        error = experimental_mz - theoretical_mz
    
    return error


def is_valid_fasta(fasta_string: str) -> bool:
    """
    Check if a string is valid FASTA format.
    
    Args:
        fasta_string: String to check
        
    Returns:
        True if valid FASTA format, False otherwise
    """
    # Basic validation
    if not fasta_string.strip():
        return False
    
    lines = fasta_string.strip().split('\n')
    
    # Must start with '>'
    if not lines[0].startswith('>'):
        return False
    
    # Must have at least one sequence line
    if len(lines) < 2:
        return False
    
    current_is_header = True
    
    for line in lines:
        line = line.strip()
        if not line:
            continue
            
        if line.startswith('>'):
            if not current_is_header:
                current_is_header = True
            else:
                # Two headers in a row (empty sequence)
                if len(lines) <= 1:
                    return False
        else:
            current_is_header = False
            
            # Check for non-amino acid characters in sequence
            valid_chars = set("ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy*- \t")
            if any(c not in valid_chars for c in line):
                return False
    
    return True


def get_common_ptms() -> Dict[str, Dict[str, Any]]:
    """
    Get a dictionary of common post-translational modifications.
    
    Returns:
        Dictionary of modification names and their properties
    """
    return {
        "Phosphorylation": {
            "targets": ["S", "T", "Y"],
            "mass_shift": 79.9663,
            "formula_add": {"P": 1, "O": 3},
            "formula_remove": {"H": 1},
            "pyopenms_name": "Phospho"
        },
        "Oxidation": {
            "targets": ["M"],
            "mass_shift": 15.9949,
            "formula_add": {"O": 1},
            "formula_remove": {},
            "pyopenms_name": "Oxidation"
        },
        "Carbamidomethyl": {
            "targets": ["C"],
            "mass_shift": 57.0215,
            "formula_add": {"C": 2, "H": 3, "N": 1, "O": 1},
            "formula_remove": {"H": 1},
            "pyopenms_name": "Carbamidomethyl"
        },
        "Acetylation": {
            "targets": ["K", "N-term"],
            "mass_shift": 42.0106,
            "formula_add": {"C": 2, "H": 2, "O": 1},
            "formula_remove": {"H": 1},
            "pyopenms_name": "Acetyl"
        },
        "Methylation": {
            "targets": ["K", "R"],
            "mass_shift": 14.0157,
            "formula_add": {"C": 1, "H": 2},
            "formula_remove": {},
            "pyopenms_name": "Methyl"
        },
        "Dimethylation": {
            "targets": ["K", "R"],
            "mass_shift": 28.0313,
            "formula_add": {"C": 2, "H": 4},
            "formula_remove": {},
            "pyopenms_name": "Dimethyl"
        }
    }


if __name__ == "__main__":
    # Example usage
    sequence = "PEPTIDE"
    
    # Get amino acid composition
    composition = get_aa_composition(sequence)
    print(f"Composition of {sequence}: {composition}")
    
    # Calculate approximate molecular weight
    mw = calculate_molecular_weight(composition)
    print(f"Approximate molecular weight: {mw:.2f} Da")
    
    # Validate sequence
    is_valid, invalid_chars = validate_sequence(sequence)
    print(f"Is valid sequence: {is_valid}")
    
    if not is_valid:
        print(f"Invalid characters: {invalid_chars}")