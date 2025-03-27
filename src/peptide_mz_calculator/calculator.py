"""
Calculator module for peptide m/z calculation.

This module provides functions for calculating mass-to-charge ratios,
generating isotope patterns, and other mass spectrometry-related calculations.
"""

import pyopenms
from pyopenms import AASequence, EmpiricalFormula
import numpy as np
from typing import Dict, List, Tuple, Optional, Union, Any


def calculate_mz(
    sequence: str,
    charge_state: int = 1,
    modifications: Optional[List[Dict[str, str]]] = None,
    isotope_type: str = "monoisotopic"
) -> Dict[str, Any]:
    """
    Calculate the m/z ratio of a peptide sequence using PyOpenMS.
    
    Args:
        sequence: Amino acid sequence in single-letter code
        charge_state: Charge state of the peptide
        modifications: List of modifications to apply
        isotope_type: Type of isotope mass to use
        
    Returns:
        Dictionary containing m/z value, mass, and other calculated properties
    
    Raises:
        ValueError: If sequence contains invalid amino acids or modifications
        RuntimeError: If calculation fails
    """
    try:
        # Validate sequence
        valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
        invalid_aa = [aa for aa in sequence if aa not in valid_aa]
        if invalid_aa:
            raise ValueError(f"Invalid amino acids in sequence: {', '.join(invalid_aa)}")
        
        # Create AASequence object
        peptide = AASequence.fromString(sequence)
        
        # Apply modifications if any
        if modifications:
            for mod in modifications:
                if mod['position'] == 'N-term':
                    peptide.setNTerminalModification(mod['name'])
                elif mod['position'] == 'C-term':
                    peptide.setCTerminalModification(mod['name'])
                else:
                    pos = int(mod['position']) - 1  # 0-based indexing
                    peptide.setModification(pos, mod['name'])
        
        # Calculate properties
        if isotope_type.lower() == "monoisotopic":
            mass = peptide.getMonoWeight()
        else:
            mass = peptide.getAverageWeight()
        
        # Calculate m/z by adding proton mass for each charge
        proton_mass = 1.007276466
        mz = (mass + charge_state * proton_mass) / charge_state
        
        # Get formula
        formula = peptide.getFormula()
        formula_string = str(formula)
        
        # Get fragment ions
        fragment_ions = generate_fragment_ions(sequence, charge_state)
        
        # Prepare result
        result = {
            "sequence": sequence,
            "charge": charge_state,
            "mass": mass,
            "mz": mz,
            "formula": formula_string,
            "fragment_ions": fragment_ions,
            "isotope_type": isotope_type
        }
        
        return result
    
    except Exception as e:
        # Re-raise as RuntimeError with more context
        raise RuntimeError(f"Error calculating m/z: {str(e)}") from e


def generate_isotope_pattern(
    peptide_sequence: str,
    charge_state: int = 1
) -> Tuple[List[float], List[float]]:
    """
    Generate theoretical isotope pattern for a peptide sequence.
    
    Args:
        peptide_sequence: Amino acid sequence
        charge_state: Charge state of the peptide
        
    Returns:
        Tuple of (mz_values, intensities) arrays for plotting
    
    Raises:
        RuntimeError: If pattern generation fails
    """
    try:
        peptide = AASequence.fromString(peptide_sequence)
        peptide_weight = peptide.getMonoWeight()
        
        # Create isotope distribution
        isotope_dist = pyopenms.CoarseIsotopePatternGenerator(5)
        iso_pattern = pyopenms.IsotopeDistribution()
        isotope_dist.estimateFromPeptideWeight(peptide_weight)
        iso_pattern = isotope_dist.getIsotopeDistribution(iso_pattern)
        
        masses = []
        intensities = []
        
        # Calculate m/z values and get intensities
        for i in range(iso_pattern.size()):
            # Convert mass to m/z by adding proton mass and dividing by charge
            proton_mass = 1.007276466
            mz = (peptide_weight + (i * 1.003) + (charge_state * proton_mass)) / charge_state
            masses.append(mz)
            intensities.append(iso_pattern.getIntensity(i))
        
        # If the pattern is empty (rare but possible), use theoretical approximation
        if not masses:
            # Base mz from the monoisotopic mass
            mono_mass = peptide.getMonoWeight(pyopenms.Residue.ResidueType.Full, charge_state)
            mz = mono_mass / charge_state
            
            # Approximate isotope pattern based on peptide size
            seq_length = len(peptide_sequence)
            
            # Pattern intensities scale with peptide size
            if seq_length < 5:
                # Small peptides have simpler patterns
                intensity_pattern = [100, 50, 15, 3, 0.5]
            elif seq_length < 10:
                # Medium peptides
                intensity_pattern = [100, 65, 30, 10, 2]
            else:
                # Larger peptides
                intensity_pattern = [100, 80, 50, 25, 10]
            
            # Scale intensities based on charge state
            intensity_pattern = [i * (1 - (0.1 * (charge_state - 1))) for i in intensity_pattern]
            intensity_pattern = [max(i, 0) for i in intensity_pattern]
            
            # Generate m/z values for each isotope peak
            masses = [mz + (i * 1.003) / charge_state for i in range(len(intensity_pattern))]
            intensities = intensity_pattern
        
        # Normalize intensities to 100%
        max_intensity = max(intensities)
        normalized_intensities = [i * 100 / max_intensity for i in intensities]
        
        return masses, normalized_intensities
    
    except Exception as e:
        raise RuntimeError(f"Error generating isotope pattern: {str(e)}") from e


def generate_fragment_ions(
    peptide_sequence: str,
    charge_state: int = 1
) -> Dict[str, List[Dict[str, Union[int, float]]]]:
    """
    Generate theoretical fragment ions (b and y) for a peptide sequence.
    
    Args:
        peptide_sequence: Amino acid sequence
        charge_state: Charge state
        
    Returns:
        Dictionary containing b_ions and y_ions lists
    
    Raises:
        RuntimeError: If fragment ion generation fails
    """
    try:
        peptide = AASequence.fromString(peptide_sequence)
        
        fragment_ions = {}
        
        # b-ions
        b_ions = []
        for i in range(1, len(peptide_sequence)):
            b_fragment = peptide.getPrefix(i)
            b_mass = b_fragment.getMonoWeight(pyopenms.Residue.ResidueType.BIon, charge_state)
            b_mz = b_mass / charge_state
            b_ions.append({
                "position": i,
                "mass": b_mass,
                "mz": b_mz,
                "sequence": str(b_fragment)
            })
        fragment_ions["b_ions"] = b_ions
        
        # y-ions
        y_ions = []
        for i in range(1, len(peptide_sequence)):
            y_fragment = peptide.getSuffix(i)
            y_mass = y_fragment.getMonoWeight(pyopenms.Residue.ResidueType.YIon, charge_state)
            y_mz = y_mass / charge_state
            y_ions.append({
                "position": i,
                "mass": y_mass,
                "mz": y_mz,
                "sequence": str(y_fragment)
            })
        fragment_ions["y_ions"] = y_ions
        
        return fragment_ions
    
    except Exception as e:
        raise RuntimeError(f"Error generating fragment ions: {str(e)}") from e


def get_theoretical_spectrum(
    peptide_sequence: str,
    min_charge: int = 1,
    max_charge: int = 3
) -> Dict[str, List[Dict[str, Any]]]:
    """
    Generate a comprehensive theoretical spectrum for a peptide.
    
    Args:
        peptide_sequence: Amino acid sequence
        min_charge: Minimum charge state to consider
        max_charge: Maximum charge state to consider
        
    Returns:
        Dictionary containing various ion types and their m/z values
    """
    try:
        peptide = AASequence.fromString(peptide_sequence)
        spectrum = {}
        
        # Precursor ions at different charge states
        precursors = []
        for z in range(min_charge, max_charge + 1):
            mass = peptide.getMonoWeight()
            mz = (mass + z * 1.007276466) / z
            precursors.append({
                "charge": z,
                "mz": mz,
                "type": "precursor"
            })
        spectrum["precursor_ions"] = precursors
        
        # Fragment ions at different charge states
        all_b_ions = []
        all_y_ions = []
        
        for z in range(min_charge, max_charge + 1):
            # Get the ions for this charge state
            fragment_ions = generate_fragment_ions(peptide_sequence, z)
            
            # Add charge information
            for ion in fragment_ions["b_ions"]:
                ion["charge"] = z
                all_b_ions.append(ion)
            
            for ion in fragment_ions["y_ions"]:
                ion["charge"] = z
                all_y_ions.append(ion)
        
        spectrum["b_ions"] = all_b_ions
        spectrum["y_ions"] = all_y_ions
        
        # Also add other ion types (a-ions, c-ions, etc.) if needed
        
        return spectrum
    
    except Exception as e:
        raise RuntimeError(f"Error generating theoretical spectrum: {str(e)}") from e


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
        # Error in parts per million
        error = (experimental_mz - theoretical_mz) / theoretical_mz * 1e6
    else:
        # Error in Daltons
        error = experimental_mz - theoretical_mz
    
    return error


def match_peaks(
    theoretical_peaks: List[float],
    experimental_peaks: List[float],
    tolerance: float = 10,
    tolerance_unit: str = "ppm"
) -> List[Dict[str, Any]]:
    """
    Match theoretical peaks to experimental peaks within a tolerance.
    
    Args:
        theoretical_peaks: List of theoretical m/z values
        experimental_peaks: List of experimental m/z values
        tolerance: Tolerance for matching
        tolerance_unit: Unit of tolerance ("ppm" or "da")
        
    Returns:
        List of matched peaks with error information
    """
    matches = []
    
    for theo_mz in theoretical_peaks:
        best_match = None
        min_error = float('inf')
        
        for exp_mz in experimental_peaks:
            error = calculate_mass_error(theo_mz, exp_mz, tolerance_unit)
            
            # Check if the error is within tolerance
            if abs(error) <= tolerance and abs(error) < abs(min_error):
                min_error = error
                best_match = exp_mz
        
        if best_match:
            matches.append({
                "theoretical_mz": theo_mz,
                "experimental_mz": best_match,
                "error": min_error,
                "error_unit": tolerance_unit
            })
    
    return matches


def get_element_composition(peptide_sequence: str) -> Dict[str, int]:
    """
    Get elemental composition of a peptide sequence.
    
    Args:
        peptide_sequence: Amino acid sequence
        
    Returns:
        Dictionary of elements and their counts
    """
    try:
        peptide = AASequence.fromString(peptide_sequence)
        formula = peptide.getFormula()
        
        composition = {}
        
        # Get the string representation of the formula
        formula_str = str(formula)
        
        # Parse the formula string to get element counts
        element = ""
        count_str = ""
        reading_element = True
        
        for char in formula_str:
            if char.isalpha():
                if not reading_element and element:
                    # Finished reading a count, save the element and count
                    count = int(count_str) if count_str else 1
                    composition[element] = count
                    element = char
                    count_str = ""
                    reading_element = True
                else:
                    # Continue reading element
                    element += char
                    reading_element = True
            elif char.isdigit():
                count_str += char
                reading_element = False
        
        # Process the last element
        if element:
            count = int(count_str) if count_str else 1
            composition[element] = count
        
        return composition
    
    except Exception as e:
        raise RuntimeError(f"Error getting element composition: {str(e)}") from e


if __name__ == "__main__":
    # Example usage
    test_peptide = "PEPTIDE"
    result = calculate_mz(test_peptide, charge_state=2)
    
    print(f"Peptide: {test_peptide}")
    print(f"m/z: {result['mz']:.4f}")
    print(f"Mass: {result['mass']:.4f}")
    print(f"Formula: {result['formula']}")
    
    # Generate isotope pattern
    masses, intensities = generate_isotope_pattern(test_peptide, charge_state=2)
    
    print("\nIsotope Pattern:")
    for i in range(len(masses)):
        print(f"  Peak {i+1}: m/z = {masses[i]:.4f}, Intensity = {intensities[i]:.1f}%")