"""
Module for handling peptide modifications.

This module provides functions for visualizing, selecting, and applying
modifications to peptide sequences.
"""

import streamlit as st
import pandas as pd
import json
import os
from typing import Dict, List, Any, Optional, Union
import pyopenms

# Default modification data
DEFAULT_MODIFICATIONS = {
    "common_modifications": [
        {
            "name": "Phosphorylation",
            "targets": ["S", "T", "Y"],
            "mass_shift": 79.9663,
            "formula_add": {"P": 1, "O": 3},
            "formula_remove": {"H": 1},
            "pyopenms_name": "Phospho"
        },
        {
            "name": "Oxidation",
            "targets": ["M"],
            "mass_shift": 15.9949,
            "formula_add": {"O": 1},
            "formula_remove": {},
            "pyopenms_name": "Oxidation"
        },
        {
            "name": "Carbamidomethyl",
            "targets": ["C"],
            "mass_shift": 57.0215,
            "formula_add": {"C": 2, "H": 3, "N": 1, "O": 1},
            "formula_remove": {"H": 1},
            "pyopenms_name": "Carbamidomethyl"
        },
        {
            "name": "Acetylation",
            "targets": ["K", "N-term"],
            "mass_shift": 42.0106,
            "formula_add": {"C": 2, "H": 2, "O": 1},
            "formula_remove": {"H": 1},
            "pyopenms_name": "Acetyl"
        },
        {
            "name": "Methylation",
            "targets": ["K", "R"],
            "mass_shift": 14.0157,
            "formula_add": {"C": 1, "H": 2},
            "formula_remove": {},
            "pyopenms_name": "Methyl"
        }
    ],
    "custom_modifications": [],
    "n_terminal_modifications": [
        {
            "name": "Acetylation",
            "mass_shift": 42.0106,
            "formula_add": {"C": 2, "H": 2, "O": 1},
            "formula_remove": {"H": 1},
            "pyopenms_name": "Acetyl"
        },
        {
            "name": "Formylation",
            "mass_shift": 27.9949,
            "formula_add": {"C": 1, "O": 1},
            "formula_remove": {"H": 1},
            "pyopenms_name": "Formyl"
        }
    ],
    "c_terminal_modifications": [
        {
            "name": "Amidation",
            "mass_shift": -0.9840,
            "formula_add": {"N": 1, "H": 1},
            "formula_remove": {"O": 1},
            "pyopenms_name": "Amidated"
        }
    ]
}


# Predefined modification sets for common experiments
MODIFICATION_SETS = {
    "Standard Proteomics": [
        {"name": "Carbamidomethyl", "targets": ["C"], "fixed": True},
        {"name": "Oxidation", "targets": ["M"], "fixed": False}
    ],
    "Phosphoproteomics": [
        {"name": "Carbamidomethyl", "targets": ["C"], "fixed": True},
        {"name": "Oxidation", "targets": ["M"], "fixed": False},
        {"name": "Phosphorylation", "targets": ["S", "T", "Y"], "fixed": False}
    ],
    "Acetylome Analysis": [
        {"name": "Carbamidomethyl", "targets": ["C"], "fixed": True},
        {"name": "Oxidation", "targets": ["M"], "fixed": False},
        {"name": "Acetylation", "targets": ["K"], "fixed": False},
        {"name": "Acetylation", "targets": ["N-term"], "fixed": False}
    ],
    "Ubiquitinome": [
        {"name": "Carbamidomethyl", "targets": ["C"], "fixed": True},
        {"name": "Oxidation", "targets": ["M"], "fixed": False},
        {"name": "GlyGly", "targets": ["K"], "fixed": False}
    ],
    "TMT Labeling": [
        {"name": "Carbamidomethyl", "targets": ["C"], "fixed": True},
        {"name": "Oxidation", "targets": ["M"], "fixed": False},
        {"name": "TMT6plex", "targets": ["K"], "fixed": True},
        {"name": "TMT6plex", "targets": ["N-term"], "fixed": True}
    ],
    "iTRAQ Labeling": [
        {"name": "Carbamidomethyl", "targets": ["C"], "fixed": True},
        {"name": "Oxidation", "targets": ["M"], "fixed": False},
        {"name": "iTRAQ4plex", "targets": ["K"], "fixed": True},
        {"name": "iTRAQ4plex", "targets": ["N-term"], "fixed": True}
    ],
    "SILAC Light/Heavy": [
        {"name": "Carbamidomethyl", "targets": ["C"], "fixed": True},
        {"name": "Oxidation", "targets": ["M"], "fixed": False},
        {"name": "Label:13C(6)15N(4)", "targets": ["R"], "fixed": False},
        {"name": "Label:13C(6)15N(2)", "targets": ["K"], "fixed": False}
    ]
}


def load_modifications(file_path: Optional[str] = None) -> Dict[str, List[Dict[str, Any]]]:
    """
    Load modifications from a JSON file or use defaults.
    
    Args:
        file_path: Path to the modifications JSON file
        
    Returns:
        Dictionary of modification lists
    """
    if file_path and os.path.exists(file_path):
        try:
            with open(file_path, 'r') as f:
                return json.load(f)
        except Exception as e:
            st.error(f"Error loading modifications: {str(e)}")
    
    return DEFAULT_MODIFICATIONS


def save_modifications(modifications: Dict[str, List[Dict[str, Any]]], file_path: str) -> bool:
    """
    Save modifications to a JSON file.
    
    Args:
        modifications: Dictionary of modification lists
        file_path: Path to save the JSON file
        
    Returns:
        True if successful, False otherwise
    """
    try:
        with open(file_path, 'w') as f:
            json.dump(modifications, f, indent=4)
        return True
    except Exception as e:
        st.error(f"Error saving modifications: {str(e)}")
        return False


def get_applicable_modifications(sequence: str, modifications: Dict[str, List[Dict[str, Any]]]) -> Dict[str, List[Dict[str, Any]]]:
    """
    Get modifications applicable to a peptide sequence.
    
    Args:
        sequence: Peptide sequence
        modifications: Dictionary of modification lists
        
    Returns:
        Dictionary of applicable modifications
    """
    applicable = {
        "residue_specific": [],
        "n_terminal": modifications.get("n_terminal_modifications", []),
        "c_terminal": modifications.get("c_terminal_modifications", [])
    }
    
    # Check common and custom modifications
    all_mods = modifications.get("common_modifications", []) + modifications.get("custom_modifications", [])
    
    for mod in all_mods:
        targets = mod.get("targets", [])
        
        # Check if any target residue is in the sequence
        if any(aa in sequence for aa in targets):
            applicable["residue_specific"].append(mod)
    
    return applicable


def formula_to_string(formula: Dict[str, int]) -> str:
    """
    Convert a formula dictionary to a string.
    
    Args:
        formula: Dictionary of elements and counts
        
    Returns:
        Formula string (e.g., "C2H3NO")
    """
    if not formula:
        return ""
    
    # Sort elements in a standard order
    standard_order = ["C", "H", "N", "O", "P", "S"]
    
    # Sort first by standard order, then alphabetically for other elements
    formula_items = sorted(formula.items(), key=lambda x: (
        standard_order.index(x[0]) if x[0] in standard_order else 100,
        x[0]
    ))
    
    # Build the formula string
    formula_str = ""
    for element, count in formula_items:
        if count == 1:
            formula_str += element
        elif count > 1:
            formula_str += f"{element}{count}"
    
    return formula_str


def render_modification_selector(
    sequence: str,
    modifications: Dict[str, List[Dict[str, Any]]],
    current_mods: List[Dict[str, Any]] = None
) -> List[Dict[str, Any]]:
    """
    Render a UI for selecting modifications to apply to a peptide.
    
    Args:
        sequence: Peptide sequence
        modifications: Dictionary of modification lists
        current_mods: Currently applied modifications
        
    Returns:
        List of selected modifications
    """
    if current_mods is None:
        current_mods = []
    
    st.markdown("### Add Modifications")
    
    # Get applicable modifications
    applicable = get_applicable_modifications(sequence, modifications)
    
    # Create tabs for different modification types
    mod_tabs = st.tabs(["Residue Specific", "N-Terminal", "C-Terminal", "Custom"])
    
    selected_mods = current_mods.copy()
    
    # Residue specific modifications
    with mod_tabs[0]:
        st.markdown("#### Residue Specific Modifications")
        
        if not applicable["residue_specific"]:
            st.info("No applicable residue-specific modifications for this sequence.")
        else:
            # Create a visual representation of the sequence with positions
            cols = st.columns(len(sequence) + 1)
            
            with cols[0]:
                st.markdown("**Pos**")
                st.markdown("**AA**")
                st.markdown("**Mod**")
            
            for i, aa in enumerate(sequence):
                with cols[i + 1]:
                    st.markdown(f"**{i + 1}**")
                    st.markdown(f"**{aa}**")
                    
                    # Find applicable modifications for this amino acid
                    aa_mods = [mod for mod in applicable["residue_specific"] 
                              if aa in mod.get("targets", [])]
                    
                    if aa_mods:
                        # Check if position already has a modification
                        existing_mod = next((m for m in selected_mods 
                                          if m.get("position") == str(i + 1)), None)
                        
                        mod_options = ["None"] + [m["name"] for m in aa_mods]
                        default_idx = 0
                        
                        if existing_mod:
                            # Find the index of the existing modification
                            try:
                                default_idx = mod_options.index(existing_mod["name"])
                            except ValueError:
                                default_idx = 0
                        
                        selected_mod = st.selectbox(
                            f"",
                            options=mod_options,
                            index=default_idx,
                            key=f"mod_aa_{i}"
                        )
                        
                        # Update modifications based on selection
                        if selected_mod != "None":
                            mod_details = next(m for m in aa_mods if m["name"] == selected_mod)
                            
                            # Remove existing modification at this position if any
                            selected_mods = [m for m in selected_mods 
                                           if m.get("position") != str(i + 1)]
                            
                            # Add the new modification
                            selected_mods.append({
                                "name": mod_details["pyopenms_name"],
                                "position": str(i + 1),
                                "residue": aa,
                                "display_name": mod_details["name"]
                            })
                        else:
                            # Remove modification at this position if "None" selected
                            selected_mods = [m for m in selected_mods 
                                           if m.get("position") != str(i + 1)]
                    else:
                        st.text("N/A")
    
    # N-Terminal modifications
    with mod_tabs[1]:
        st.markdown("#### N-Terminal Modifications")
        
        if not applicable["n_terminal"]:
            st.info("No N-terminal modifications available.")
        else:
            # Check if N-terminus already has a modification
            existing_mod = next((m for m in selected_mods 
                              if m.get("position") == "N-term"), None)
            
            mod_options = ["None"] + [m["name"] for m in applicable["n_terminal"]]
            default_idx = 0
            
            if existing_mod:
                # Find the index of the existing modification
                try:
                    # Convert from pyopenms_name back to display name
                    existing_name = next((m["name"] for m in applicable["n_terminal"] 
                                        if m["pyopenms_name"] == existing_mod["name"]), None)
                    if existing_name:
                        default_idx = mod_options.index(existing_name)
                except ValueError:
                    default_idx = 0
            
            selected_mod = st.selectbox(
                "Select N-terminal modification:",
                options=mod_options,
                index=default_idx,
                key="mod_nterm"
            )
            
            # Update modifications based on selection
            if selected_mod != "None":
                mod_details = next(m for m in applicable["n_terminal"] if m["name"] == selected_mod)
                
                # Remove existing N-terminal modification if any
                selected_mods = [m for m in selected_mods 
                               if m.get("position") != "N-term"]
                
                # Add the new modification
                selected_mods.append({
                    "name": mod_details["pyopenms_name"],
                    "position": "N-term",
                    "display_name": mod_details["name"]
                })
            else:
                # Remove N-terminal modification if "None" selected
                selected_mods = [m for m in selected_mods 
                               if m.get("position") != "N-term"]
    
    # C-Terminal modifications
    with mod_tabs[2]:
        st.markdown("#### C-Terminal Modifications")
        
        if not applicable["c_terminal"]:
            st.info("No C-terminal modifications available.")
        else:
            # Check if C-terminus already has a modification
            existing_mod = next((m for m in selected_mods 
                              if m.get("position") == "C-term"), None)
            
            mod_options = ["None"] + [m["name"] for m in applicable["c_terminal"]]
            default_idx = 0
            
            if existing_mod:
                # Find the index of the existing modification
                try:
                    # Convert from pyopenms_name back to display name
                    existing_name = next((m["name"] for m in applicable["c_terminal"] 
                                        if m["pyopenms_name"] == existing_mod["name"]), None)
                    if existing_name:
                        default_idx = mod_options.index(existing_name)
                except ValueError:
                    default_idx = 0
            
            selected_mod = st.selectbox(
                "Select C-terminal modification:",
                options=mod_options,
                index=default_idx,
                key="mod_cterm"
            )
            
            # Update modifications based on selection
            if selected_mod != "None":
                mod_details = next(m for m in applicable["c_terminal"] if m["name"] == selected_mod)
                
                # Remove existing C-terminal modification if any
                selected_mods = [m for m in selected_mods 
                               if m.get("position") != "C-term"]
                
                # Add the new modification
                selected_mods.append({
                    "name": mod_details["pyopenms_name"],
                    "position": "C-term",
                    "display_name": mod_details["name"]
                })
            else:
                # Remove C-terminal modification if "None" selected
                selected_mods = [m for m in selected_mods 
                               if m.get("position") != "C-term"]
    
    # Custom modification editor
    with mod_tabs[3]:
        st.markdown("#### Custom Modification")
        
        with st.form("custom_mod_form"):
            st.markdown("Define a new custom modification")
            
            custom_name = st.text_input("Modification Name", "")
            
            # Target residues (multi-select)
            all_aas = list("ACDEFGHIKLMNPQRSTVWY")
            target_residues = st.multiselect(
                "Target Residues",
                options=all_aas + ["N-term", "C-term"],
                default=[]
            )
            
            # Mass shift
            mass_shift = st.number_input("Mass Shift (Da)", value=0.0, format="%.4f")
            
            # Formula addition
            st.markdown("##### Formula Addition")
            formula_add_cols = st.columns(6)
            formula_add = {}
            
            for i, elem in enumerate(["C", "H", "N", "O", "P", "S"]):
                with formula_add_cols[i % 6]:
                    count = st.number_input(f"{elem}+", value=0, min_value=0, key=f"add_{elem}")
                    if count > 0:
                        formula_add[elem] = count
            
            # Formula removal
            st.markdown("##### Formula Removal")
            formula_rem_cols = st.columns(6)
            formula_remove = {}
            
            for i, elem in enumerate(["C", "H", "N", "O", "P", "S"]):
                with formula_rem_cols[i % 6]:
                    count = st.number_input(f"{elem}-", value=0, min_value=0, key=f"rem_{elem}")
                    if count > 0:
                        formula_remove[elem] = count
            
            # Submit button
            submit_custom = st.form_submit_button("Add Custom Modification")
            
            if submit_custom and custom_name and target_residues:
                # Create custom modification
                custom_mod = {
                    "name": custom_name,
                    "targets": [r for r in target_residues if r not in ["N-term", "C-term"]],
                    "mass_shift": mass_shift,
                    "formula_add": formula_add,
                    "formula_remove": formula_remove,
                    "pyopenms_name": custom_name.replace(" ", "_")
                }
                
                # Add to custom modifications
                if "custom_modifications" not in modifications:
                    modifications["custom_modifications"] = []
                
                modifications["custom_modifications"].append(custom_mod)
                
                # Add to terminal modifications if applicable
                if "N-term" in target_residues:
                    if "n_terminal_modifications" not in modifications:
                        modifications["n_terminal_modifications"] = []
                    
                    n_term_mod = custom_mod.copy()
                    n_term_mod.pop("targets", None)
                    modifications["n_terminal_modifications"].append(n_term_mod)
                
                if "C-term" in target_residues:
                    if "c_terminal_modifications" not in modifications:
                        modifications["c_terminal_modifications"] = []
                    
                    c_term_mod = custom_mod.copy()
                    c_term_mod.pop("targets", None)
                    modifications["c_terminal_modifications"].append(c_term_mod)
                
                st.success(f"Added custom modification: {custom_name}")
        
        # Display current custom modifications
        st.markdown("#### Current Custom Modifications")
        
        custom_mods = modifications.get("custom_modifications", [])
        
        if not custom_mods:
            st.info("No custom modifications defined.")
        else:
            # Create a DataFrame for custom mods
            custom_df = []
            
            for mod in custom_mods:
                custom_df.append({
                    "Name": mod["name"],
                    "Targets": ", ".join(mod.get("targets", [])),
                    "Mass Shift": f"{mod['mass_shift']:.4f}",
                    "Formula Add": formula_to_string(mod.get("formula_add", {})),
                    "Formula Remove": formula_to_string(mod.get("formula_remove", {}))
                })
            
            st.dataframe(pd.DataFrame(custom_df))
            
            # Option to remove custom modifications
            if st.button("Remove All Custom Modifications"):
                modifications["custom_modifications"] = []
                st.success("All custom modifications removed.")
    
    # Display summary of selected modifications
    st.markdown("### Applied Modifications")
    
    if not selected_mods:
        st.info("No modifications applied.")
    else:
        mod_df = []
        
        for mod in selected_mods:
            position = mod.get("position", "")
            residue = mod.get("residue", "")
            
            if position == "N-term":
                position_display = "N-terminus"
            elif position == "C-term":
                position_display = "C-terminus"
            else:
                position_display = f"{residue}{position}"
            
            mod_df.append({
                "Modification": mod.get("display_name", mod.get("name", "")),
                "Position": position_display
            })
        
        st.dataframe(pd.DataFrame(mod_df))
    
    return selected_mods


def register_custom_modification_with_pyopenms(mod_def: Dict[str, Any]) -> bool:
    """
    Register a custom modification with PyOpenMS.
    
    Args:
        mod_def: Modification definition
        
    Returns:
        True if successful, False otherwise
    """
    try:
        # Get modification name
        mod_name = mod_def.get("pyopenms_name", mod_def.get("name", "")).replace(" ", "_")
        
        # Create formula strings
        formula_add = formula_to_string(mod_def.get("formula_add", {}))
        formula_remove = formula_to_string(mod_def.get("formula_remove", {}))
        
        # Get modification sites
        sites = "".join(mod_def.get("targets", []))
        
        # Create a modification
        mod = pyopenms.ResidueModification()
        mod.setId(mod_name)
        mod.setFullName(mod_def.get("name", ""))
        mod.setDiffFormula(pyopenms.EmpiricalFormula(formula_add))
        mod.setDiffMonoMass(mod_def.get("mass_shift", 0.0))
        
        if sites:
            mod.setOrigin(sites[0])  # Use the first target as the primary one
        
        # Add the modification to the database
        mod_db = pyopenms.ModificationsDB.getInstance()
        mod_db.addModification(mod)
        
        return True
    
    except Exception as e:
        st.error(f"Error registering modification with PyOpenMS: {str(e)}")
        return False


def visualize_modified_peptide(sequence: str, modifications: List[Dict[str, Any]]):
    """
    Visualize a peptide with its modifications.
    
    Args:
        sequence: Peptide sequence
        modifications: List of applied modifications
    """
    st.markdown("### Modified Peptide Visualization")
    
    # Create a colored sequence display
    html = """
    <style>
        .peptide-container {
            font-family: monospace;
            font-size: 20px;
            letter-spacing: 2px;
            margin: 20px 0;
        }
        .aa {
            padding: 5px;
            margin: 2px;
            border-radius: 4px;
            display: inline-block;
            width: 30px;
            height: 30px;
            text-align: center;
            line-height: 30px;
        }
        .modified {
            position: relative;
        }
        .modified::after {
            content: "*";
            position: absolute;
            top: -10px;
            right: -5px;
            color: red;
            font-weight: bold;
        }
        .nterm {
            border-left: 3px solid blue;
        }
        .cterm {
            border-right: 3px solid green;
        }
        .mod-label {
            font-size: 12px;
            margin-top: 5px;
            color: #666;
        }
        .mod-tooltip {
            position: relative;
            display: inline-block;
        }
        .mod-tooltip .tooltiptext {
            visibility: hidden;
            width: 120px;
            background-color: black;
            color: #fff;
            text-align: center;
            padding: 5px;
            border-radius: 6px;
            position: absolute;
            z-index: 1;
            bottom: 125%;
            left: 50%;
            margin-left: -60px;
            opacity: 0;
            transition: opacity 0.3s;
        }
        .mod-tooltip:hover .tooltiptext {
            visibility: visible;
            opacity: 1;
        }
    </style>
    """
    
    # Start HTML for the peptide container
    html += '<div class="peptide-container">'
    
    # Process each amino acid
    for i, aa in enumerate(sequence):
        # Check if this position has a modification
        pos_mods = [m for m in modifications if m.get("position") == str(i + 1)]
        
        # Determine CSS classes
        css_classes = ["aa"]
        
        if pos_mods:
            css_classes.append("modified")
        
        if i == 0 and any(m.get("position") == "N-term" for m in modifications):
            css_classes.append("nterm")
        
        if i == len(sequence) - 1 and any(m.get("position") == "C-term" for m in modifications):
            css_classes.append("cterm")
        
        # Set background color based on amino acid properties
        aa_colors = {
            # Hydrophobic
            'A': '#FFCC99', 'V': '#FFCC99', 'L': '#FFCC99', 'I': '#FFCC99',
            'M': '#FFCC99', 'F': '#FFCC99', 'W': '#FFCC99', 'P': '#FFCC99',
            # Polar
            'G': '#99CCFF', 'S': '#99CCFF', 'T': '#99CCFF', 'C': '#99CCFF',
            'Y': '#99CCFF', 'N': '#99CCFF', 'Q': '#99CCFF',
            # Acidic
            'D': '#FF9999', 'E': '#FF9999',
            # Basic
            'K': '#CC99FF', 'R': '#CC99FF', 'H': '#CC99FF'
        }
        
        bg_color = aa_colors.get(aa, '#FFFFFF')
        
        # Create tooltip content for modifications
        tooltip = ""
        if pos_mods:
            mod_names = [m.get("display_name", m.get("name", "Unknown")) for m in pos_mods]
            tooltip = f'<span class="tooltiptext">{", ".join(mod_names)}</span>'
        
        # Add the amino acid with position number
        html += f'<div class="mod-tooltip"><div class="{" ".join(css_classes)}" style="background-color: {bg_color}">{aa}</div>'
        html += f'<div class="mod-label">{i + 1}</div>{tooltip}</div>'
    
    html += '</div>'
    
    # N-terminal modification indicator
    nterm_mods = [m for m in modifications if m.get("position") == "N-term"]
    if nterm_mods:
        mod_names = [m.get("display_name", m.get("name", "Unknown")) for m in nterm_mods]
        html += f'<div><strong>N-terminal:</strong> {", ".join(mod_names)}</div>'
    
    # C-terminal modification indicator
    cterm_mods = [m for m in modifications if m.get("position") == "C-term"]
    if cterm_mods:
        mod_names = [m.get("display_name", m.get("name", "Unknown")) for m in cterm_mods]
        html += f'<div><strong>C-terminal:</strong> {", ".join(mod_names)}</div>'
    
    # Render the HTML
    st.markdown(html, unsafe_allow_html=True)
    
    # Show modified sequence as text
    if modifications:
        # Create the sequence with modification syntax
        modified_seq = sequence
        
        # Function to add modifications to the sequence string
        def add_mod_to_sequence(seq, position, mod_name):
            if position == "N-term":
                return f"{mod_name}-{seq}"
            elif position == "C-term":
                return f"{seq}-{mod_name}"
            else:
                pos = int(position)
                return f"{seq[:pos]}{seq[pos-1]}({mod_name}){seq[pos:]}"
        
        # Start with the unmodified sequence
        modified_seq_str = sequence
        
        # Add each modification
        for mod in modifications:
            pos = mod.get("position", "")
            name = mod.get("name", "")
            
            if pos and name:
                modified_seq_str = add_mod_to_sequence(modified_seq_str, pos, name)
        
        st.text_area("Modified Sequence (OpenMS format)", modified_seq_str, height=100)
        
        # Calculate mass and m/z with modifications
        try:
            mod_peptide = pyopenms.AASequence.fromString(modified_seq_str)
            mono_mass = mod_peptide.getMonoWeight()
            
            st.markdown(f"**Monoisotopic Mass with Modifications:** {mono_mass:.4f} Da")
        except Exception as e:
            st.warning(f"Unable to calculate modified mass: {str(e)}")
            st.markdown("Note: Some custom modifications may require registration with PyOpenMS.")


def apply_modification_set(sequence: str, mod_set_name: str, all_modifications: Dict[str, List[Dict[str, Any]]]) -> List[Dict[str, Any]]:
    """
    Apply a predefined modification set to a peptide sequence.
    
    Args:
        sequence: Peptide sequence
        mod_set_name: Name of the modification set
        all_modifications: Dictionary of all available modifications
        
    Returns:
        List of applied modifications
    """
    if mod_set_name not in MODIFICATION_SETS:
        return []
    
    # Get the modification set
    mod_set = MODIFICATION_SETS[mod_set_name]
    
    # Create applied modifications list
    applied_mods = []
    
    # Get common modifications lookup
    common_mods = {m["name"]: m for m in all_modifications.get("common_modifications", [])}
    
    # Process each modification in the set
    for mod_template in mod_set:
        mod_name = mod_template["name"]
        targets = mod_template["targets"]
        
        # Get the modification details
        mod_details = common_mods.get(mod_name)
        if not mod_details:
            continue
        
        # Apply to specified targets
        for target in targets:
            if target == "N-term":
                # Apply to N-terminus
                applied_mods.append({
                    "name": mod_details["pyopenms_name"],
                    "position": "N-term",
                    "display_name": mod_details["name"]
                })
            elif target == "C-term":
                # Apply to C-terminus
                applied_mods.append({
                    "name": mod_details["pyopenms_name"],
                    "position": "C-term",
                    "display_name": mod_details["name"]
                })
            else:
                # Apply to specific residues
                for i, aa in enumerate(sequence):
                    if aa == target:
                        applied_mods.append({
                            "name": mod_details["pyopenms_name"],
                            "position": str(i + 1),
                            "residue": aa,
                            "display_name": mod_details["name"]
                        })
    
    return applied_mods


def render_modification_set_selector(sequence: str, all_modifications: Dict[str, List[Dict[str, Any]]]) -> List[Dict[str, Any]]:
    """
    Render a UI for selecting a predefined modification set.
    
    Args:
        sequence: Peptide sequence
        all_modifications: Dictionary of all available modifications
        
    Returns:
        List of selected modifications based on the chosen set
    """
    st.markdown("### Apply Predefined Modification Set")
    
    mod_set_names = ["None"] + list(MODIFICATION_SETS.keys())
    selected_set = st.selectbox(
        "Select a predefined modification set:",
        options=mod_set_names,
        index=0
    )
    
    if selected_set != "None":
        # Apply the selected modification set
        applied_mods = apply_modification_set(sequence, selected_set, all_modifications)
        
        # Show description of the modification set
        if selected_set in MODIFICATION_SETS:
            st.markdown("#### Set Contents")
            
            # Create a DataFrame of modifications in this set
            set_df = []
            for mod in MODIFICATION_SETS[selected_set]:
                set_df.append({
                    "Modification": mod["name"],
                    "Targets": ", ".join(mod["targets"]),
                    "Type": "Fixed" if mod.get("fixed", False) else "Variable"
                })
            
            st.dataframe(pd.DataFrame(set_df))
            
            return applied_mods
    
    return []


def get_modification_mass_shift(modifications: List[Dict[str, Any]], all_modifications: Dict[str, List[Dict[str, Any]]]) -> float:
    """
    Calculate the total mass shift from all applied modifications.
    
    Args:
        modifications: List of applied modifications
        all_modifications: Dictionary of all available modifications
        
    Returns:
        Total mass shift in Daltons
    """
    total_shift = 0.0
    
    # Create a lookup dictionary for modification details
    mod_lookup = {}
    
    # Add common modifications
    for mod in all_modifications.get("common_modifications", []):
        mod_lookup[mod["name"]] = mod
        mod_lookup[mod["pyopenms_name"]] = mod
    
    # Add custom modifications
    for mod in all_modifications.get("custom_modifications", []):
        mod_lookup[mod["name"]] = mod
        mod_lookup[mod["pyopenms_name"]] = mod
    
    # Add N-terminal modifications
    for mod in all_modifications.get("n_terminal_modifications", []):
        mod_lookup[mod["name"]] = mod
        mod_lookup[mod["pyopenms_name"]] = mod
    
    # Add C-terminal modifications
    for mod in all_modifications.get("c_terminal_modifications", []):
        mod_lookup[mod["name"]] = mod
        mod_lookup[mod["pyopenms_name"]] = mod
    
    # Calculate total mass shift
    for mod in modifications:
        mod_name = mod.get("name", "")
        mod_details = mod_lookup.get(mod_name)
        
        if mod_details:
            total_shift += mod_details.get("mass_shift", 0.0)
    
    return total_shift


def format_modified_sequence(sequence: str, modifications: List[Dict[str, Any]]) -> str:
    """
    Format a sequence with modifications in OpenMS format.
    
    Args:
        sequence: Peptide sequence
        modifications: List of applied modifications
        
    Returns:
        Formatted sequence string
    """
    if not modifications:
        return sequence
    
    # Function to add modifications to the sequence string
    def add_mod_to_sequence(seq, position, mod_name):
        if position == "N-term":
            return f"{mod_name}-{seq}"
        elif position == "C-term":
            return f"{seq}-{mod_name}"
        else:
            pos = int(position)
            return f"{seq[:pos-1]}{seq[pos-1]}({mod_name}){seq[pos:]}"
    
    # Sort modifications by position (N-term first, then numbered positions, C-term last)
    def mod_sort_key(mod):
        pos = mod.get("position", "")
        if pos == "N-term":
            return -1
        elif pos == "C-term":
            return 1000
        else:
            try:
                return int(pos)
            except:
                return 500
    
    sorted_mods = sorted(modifications, key=mod_sort_key)
    
    # Start with the unmodified sequence
    result = sequence
    
    # Add each modification
    for mod in sorted_mods:
        pos = mod.get("position", "")
        name = mod.get("name", "")
        
        if pos and name:
            result = add_mod_to_sequence(result, pos, name)
    
    return result