"""
Parser module for handling different file formats.

This module provides functions for parsing FASTA, CSV, and other file formats
commonly used in proteomics research.
"""

import pandas as pd
import re
import io
from typing import List, Dict, Any, Optional, Union, TextIO


def parse_fasta(fasta_content: str) -> List[Dict[str, str]]:
    """
    Parse a FASTA file and extract sequences.
    
    Args:
        fasta_content: Content of the FASTA file
        
    Returns:
        List of dictionaries with header and sequence
    
    Examples:
        >>> content = ">Protein1\\nACDEFG\\n>Protein2\\nPEPTIDE"
        >>> parse_fasta(content)
        [{'header': 'Protein1', 'sequence': 'ACDEFG'}, 
         {'header': 'Protein2', 'sequence': 'PEPTIDE'}]
    """
    sequences = []
    current_header = None
    current_sequence = ""
    
    for line in fasta_content.split('\n'):
        line = line.strip()
        if not line:
            continue
        
        if line.startswith('>'):
            if current_header is not None:
                sequences.append({
                    'header': current_header,
                    'sequence': current_sequence
                })
            current_header = line[1:]
            current_sequence = ""
        else:
            current_sequence += line
    
    if current_header is not None:
        sequences.append({
            'header': current_header,
            'sequence': current_sequence
        })
    
    return sequences


def validate_fasta_format(fasta_content: str) -> bool:
    """
    Validate whether a string is in proper FASTA format.
    
    Args:
        fasta_content: Content to validate
        
    Returns:
        True if valid FASTA format, False otherwise
    """
    # Check if content starts with '>'
    if not fasta_content.strip().startswith('>'):
        return False
    
    # Count headers
    headers = re.findall(r'^>', fasta_content, re.MULTILINE)
    if not headers:
        return False
    
    # Parse and check each entry
    entries = parse_fasta(fasta_content)
    for entry in entries:
        # Check for empty headers or sequences
        if not entry['header'] or not entry['sequence']:
            return False
        
        # Check for valid amino acid characters
        valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
        if not all(aa in valid_aa for aa in entry['sequence']):
            return False
    
    return True


def parse_csv(
    csv_content: str,
    delimiter: str = ','
) -> List[Dict[str, Any]]:
    """
    Parse a CSV file containing peptide sequences.
    
    Args:
        csv_content: Content of the CSV file
        delimiter: Delimiter character in the CSV file
        
    Returns:
        List of dictionaries containing the column headers and values from the CSV
    
    Examples:
        >>> content = "Sequence,Mass\\nPEPTIDE,799.4\\nACDEFG,668.2"
        >>> parse_csv(content)
        [{'Sequence': 'PEPTIDE', 'Mass': 799.4}, 
         {'Sequence': 'ACDEFG', 'Mass': 668.2}]
    """
    try:
        # Create a buffer from the content
        buffer = io.StringIO(csv_content)
        
        # Read CSV into DataFrame
        df = pd.read_csv(buffer, sep=delimiter)
        
        # Convert DataFrame to list of dictionaries
        records = df.to_dict('records')
        
        return records
    
    except Exception as e:
        # Re-raise with context
        raise ValueError(f"Error parsing CSV: {str(e)}") from e


def parse_txt(txt_content: str) -> List[str]:
    """
    Parse a plain text file with one sequence per line.
    
    Args:
        txt_content: Content of the text file
        
    Returns:
        List of sequences
    
    Examples:
        >>> content = "PEPTIDE\\nACDEFG\\nKLMNPQ"
        >>> parse_txt(content)
        ['PEPTIDE', 'ACDEFG', 'KLMNPQ']
    """
    sequences = []
    
    for line in txt_content.split('\n'):
        line = line.strip()
        if line:
            sequences.append(line)
    
    return sequences


def detect_file_format(file_content: str) -> str:
    """
    Detect the format of a file based on its content.
    
    Args:
        file_content: Content of the file
        
    Returns:
        Format string ("fasta", "csv", "tsv", or "txt")
    """
    # Check for FASTA format
    if file_content.strip().startswith('>'):
        if validate_fasta_format(file_content):
            return "fasta"
    
    # Check for CSV/TSV
    first_line = file_content.split('\n')[0].strip()
    if ',' in first_line:
        return "csv"
    elif '\t' in first_line:
        return "tsv"
    
    # Default to plain text
    return "txt"


def find_sequence_column(csv_records: List[Dict[str, Any]]) -> Optional[str]:
    """
    Find the column that likely contains peptide sequences.
    
    Args:
        csv_records: List of dictionaries from CSV parsing
        
    Returns:
        Column name, or None if no sequence column found
    """
    if not csv_records:
        return None
    
    # Potential column names for sequences
    possible_names = ['Sequence', 'sequence', 'SEQUENCE', 'Peptide', 'peptide', 'PEPTIDE']
    
    # Check for exact column name match
    columns = csv_records[0].keys()
    for name in possible_names:
        if name in columns:
            return name
    
    # Check for partial match
    for column in columns:
        for name in possible_names:
            if name.lower() in column.lower():
                return column
    
    # If nothing matches, look for a column with valid sequences
    valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
    
    for column in columns:
        # Check the first few records
        sample_values = [record[column] for record in csv_records[:5] if column in record]
        
        # Convert to strings
        sample_values = [str(value).upper() for value in sample_values]
        
        # Check if values look like peptide sequences
        if all(isinstance(val, str) and all(aa in valid_aa for aa in val) for val in sample_values):
            return column
    
    # No good candidates found
    return None


def write_fasta(sequences: List[Dict[str, str]]) -> str:
    """
    Generate FASTA format content from a list of sequences.
    
    Args:
        sequences: List of dictionaries with header and sequence
        
    Returns:
        FASTA formatted string
    """
    fasta_content = ""
    
    for seq in sequences:
        header = seq.get('header', 'Sequence')
        sequence = seq.get('sequence', '')
        
        # Ensure header starts with '>'
        if not header.startswith('>'):
            header = '>' + header
        
        # Add the entry
        fasta_content += f"{header}\n{sequence}\n"
    
    return fasta_content


def parse_file(
    file_content: str,
    file_format: Optional[str] = None
) -> List[Dict[str, Any]]:
    """
    Parse a file based on its format.
    
    Args:
        file_content: Content of the file
        file_format: Format of the file (optional, auto-detects if None)
        
    Returns:
        List of parsed records
    """
    # Auto-detect format if not provided
    if file_format is None:
        file_format = detect_file_format(file_content)
    
    # Parse based on format
    if file_format == "fasta":
        return parse_fasta(file_content)
    
    elif file_format == "csv":
        records = parse_csv(file_content, delimiter=',')
        sequence_column = find_sequence_column(records)
        
        if sequence_column:
            # Extract just the sequences
            return [{'sequence': record[sequence_column]} for record in records if sequence_column in record]
        return records
    
    elif file_format == "tsv":
        records = parse_csv(file_content, delimiter='\t')
        sequence_column = find_sequence_column(records)
        
        if sequence_column:
            # Extract just the sequences
            return [{'sequence': record[sequence_column]} for record in records if sequence_column in record]
        return records
    
    elif file_format == "txt":
        sequences = parse_txt(file_content)
        return [{'sequence': seq} for seq in sequences]
    
    else:
        raise ValueError(f"Unsupported file format: {file_format}")


if __name__ == "__main__":
    # Example usage
    fasta_example = ">Protein1\nACDEFG\n>Protein2\nPEPTIDE"
    parsed_fasta = parse_fasta(fasta_example)
    print("Parsed FASTA:", parsed_fasta)
    
    csv_example = "Sequence,Mass\nPEPTIDE,799.4\nACDEFG,668.2"
    parsed_csv = parse_csv(csv_example)
    print("Parsed CSV:", parsed_csv)
    
    txt_example = "PEPTIDE\nACDEFG\nKLMNPQ"
    parsed_txt = parse_txt(txt_example)
    print("Parsed TXT:", parsed_txt)