"""
Visualization module for peptide m/z calculator.

This module provides functions for creating interactive visualizations
of mass spectra and other related plots.
"""

import streamlit as st
from typing import Dict, List, Any, Optional, Union
import plotly.graph_objects as go
import plotly.express as px
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from typing import List, Dict, Tuple, Optional, Any, Union


def plot_spectrum(
    mz_values: List[float],
    intensities: List[float],
    title: str = "Theoretical Mass Spectrum",
    plot_type: str = "lines",
    color: str = "#1E88E5",
    show_peaks: bool = True,
    highlight_monoisotopic: bool = True
) -> go.Figure:
    """
    Create an interactive plot of a mass spectrum.
    
    Args:
        mz_values: m/z values for the x-axis
        intensities: Corresponding intensity values for the y-axis
        title: Plot title
        plot_type: Type of plot ("lines", "stems", or "both")
        color: Color for the plot
        show_peaks: Whether to show peak annotations
        highlight_monoisotopic: Whether to highlight the monoisotopic peak
        
    Returns:
        Plotly figure object that can be displayed or saved
    """
    fig = go.Figure()
    
    if plot_type in ["lines", "both"]:
        fig.add_trace(go.Scatter(
            x=mz_values,
            y=intensities,
            mode='lines',
            name='Isotope Pattern',
            line=dict(color=color, width=2)
        ))
    
    if plot_type in ["stems", "both"]:
        for i, (mz, intensity) in enumerate(zip(mz_values, intensities)):
            fig.add_trace(go.Scatter(
                x=[mz, mz],
                y=[0, intensity],
                mode='lines',
                line=dict(color=color, width=1),
                showlegend=False,
                hoverinfo='skip'
            ))
    
    # Add peak annotations if requested
    if show_peaks:
        for i, (mz, intensity) in enumerate(zip(mz_values, intensities)):
            # Only label peaks with significant intensity
            if intensity > max(intensities) * 0.10:
                peak_color = "red" if i == 0 and highlight_monoisotopic else color
                fig.add_trace(go.Scatter(
                    x=[mz],
                    y=[intensity],
                    mode='markers+text',
                    marker=dict(color=peak_color, size=8),
                    text=[f"{mz:.4f}"],
                    textposition="top center",
                    showlegend=False,
                    hoverinfo='text',
                    hovertext=f"m/z: {mz:.4f}<br>Intensity: {intensity:.1f}%"
                ))
    
    fig.update_layout(
        title=title,
        xaxis_title='m/z',
        yaxis_title='Relative Intensity (%)',
        template='plotly_white',
        height=500,
        margin=dict(l=40, r=40, t=60, b=40),
        showlegend=False,
        hovermode='closest'
    )
    
    # Add grid lines for better readability
    fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='lightgray')
    fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='lightgray')
    
    return fig


def plot_fragment_ions(
    peptide_sequence: str,
    fragment_ions: Dict[str, List[Dict[str, Any]]],
    title: str = "Fragment Ions"
) -> go.Figure:
    """
    Create a plot of fragment ions.
    
    Args:
        peptide_sequence: Amino acid sequence
        fragment_ions: Dictionary of fragment ions
        title: Plot title
        
    Returns:
        Plotly figure object
    """
    fig = go.Figure()
    
    # Extract b-ions and y-ions
    b_ions = fragment_ions.get("b_ions", [])
    y_ions = fragment_ions.get("y_ions", [])
    
    # Add b-ions to the plot
    b_mz_values = [ion["mz"] for ion in b_ions]
    b_intensities = [80] * len(b_mz_values)  # Arbitrary intensity for visualization
    
    fig.add_trace(go.Scatter(
        x=b_mz_values,
        y=b_intensities,
        mode='markers+text',
        marker=dict(color="#1E88E5", size=10, symbol="circle"),
        text=[f"b{ion['position']}" for ion in b_ions],
        textposition="top center",
        name='b-ions'
    ))
    
    # Add y-ions to the plot
    y_mz_values = [ion["mz"] for ion in y_ions]
    y_intensities = [40] * len(y_mz_values)  # Lower for visual separation
    
    fig.add_trace(go.Scatter(
        x=y_mz_values,
        y=y_intensities,
        mode='markers+text',
        marker=dict(color="#FFA000", size=10, symbol="square"),
        text=[f"y{ion['position']}" for ion in y_ions],
        textposition="bottom center",
        name='y-ions'
    ))
    
    fig.update_layout(
        title=title,
        xaxis_title='m/z',
        yaxis_title='',
        template='plotly_white',
        height=300,
        showlegend=True,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        hovermode='closest'
    )
    
    # Hide y-axis ticks since intensities are arbitrary
    fig.update_yaxes(showticklabels=False)
    
    return fig


def plot_aa_composition(
    composition: Dict[str, int]
) -> go.Figure:
    """
    Create a bar chart of amino acid composition.
    
    Args:
        composition: Dictionary of amino acids and their counts
        
    Returns:
        Plotly figure object
    """
    # Convert composition to lists for plotting
    aa_codes = list(composition.keys())
    counts = list(composition.values())
    
    # Create bar chart
    fig = go.Figure(data=[
        go.Bar(
            x=aa_codes,
            y=counts,
            marker_color="#1E88E5",
            text=counts,
            textposition='auto'
        )
    ])
    
    fig.update_layout(
        title="Amino Acid Composition",
        xaxis_title="Amino Acid",
        yaxis_title="Count",
        template="plotly_white",
        height=400
    )
    
    return fig


def plot_mass_distribution(
    masses: List[float],
    title: str = "Mass Distribution"
) -> go.Figure:
    """
    Create a histogram of peptide masses.
    
    Args:
        masses: List of peptide masses
        title: Plot title
        
    Returns:
        Plotly figure object
    """
    fig = go.Figure()
    
    fig.add_trace(go.Histogram(
        x=masses,
        nbinsx=20,
        marker_color='#1E88E5',
        opacity=0.7
    ))
    
    fig.update_layout(
        title=title,
        xaxis_title="Mass (Da)",
        yaxis_title="Count",
        template="plotly_white",
        height=400
    )
    
    return fig


def plot_comparison(
    spectra: List[Dict[str, Union[str, List[float], List[float]]]],
    title: str = "Spectrum Comparison"
) -> go.Figure:
    """
    Create a comparison plot of multiple spectra.
    
    Args:
        spectra: List of dictionaries with 'name', 'mz_values', 'intensities'
        title: Plot title
        
    Returns:
        Plotly figure object
    """
    fig = go.Figure()
    
    colors = ["#1E88E5", "#FFA000", "#43A047", "#7B1FA2", "#E53935"]
    
    for i, spectrum in enumerate(spectra):
        name = spectrum.get("name", f"Spectrum {i+1}")
        mz_values = spectrum.get("mz_values", [])
        intensities = spectrum.get("intensities", [])
        
        color = colors[i % len(colors)]
        
        fig.add_trace(go.Scatter(
            x=mz_values,
            y=intensities,
            mode='lines',
            name=name,
            line=dict(color=color, width=2)
        ))
    
    fig.update_layout(
        title=title,
        xaxis_title='m/z',
        yaxis_title='Relative Intensity (%)',
        template='plotly_white',
        height=500,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        hovermode='closest'
    )
    
    return fig


def plot_sequence_coverage(
    peptide_sequence: str,
    covered_positions: List[int],
    title: str = "Sequence Coverage"
) -> go.Figure:
    """
    Create a visual representation of sequence coverage.
    
    Args:
        peptide_sequence: Amino acid sequence
        covered_positions: List of positions (0-based) covered
        title: Plot title
        
    Returns:
        Plotly figure object
    """
    fig = go.Figure()
    
    # Create a list of amino acids with position
    aa_positions = [(i, aa) for i, aa in enumerate(peptide_sequence)]
    
    # Calculate coverage percentage
    coverage_percent = len(covered_positions) / len(peptide_sequence) * 100
    
    # Add amino acids as text labels
    for i, aa in aa_positions:
        is_covered = i in covered_positions
        color = "#43A047" if is_covered else "#E53935"
        
        fig.add_trace(go.Scatter(
            x=[i],
            y=[0],
            mode='markers+text',
            marker=dict(color=color, size=30, symbol="square"),
            text=[aa],
            textposition="middle center",
            textfont=dict(color="white", size=14),
            showlegend=False,
            hoverinfo='text',
            hovertext=f"Position: {i+1}<br>Amino Acid: {aa}<br>{'Covered' if is_covered else 'Not Covered'}"
        ))
    
    # Add coverage percentage as annotation
    fig.add_annotation(
        x=len(peptide_sequence) / 2,
        y=1,
        text=f"Coverage: {coverage_percent:.1f}%",
        showarrow=False,
        font=dict(size=14)
    )
    
    fig.update_layout(
        title=title,
        xaxis_title='Position',
        template='plotly_white',
        height=200,
        margin=dict(l=40, r=40, t=60, b=40),
        xaxis=dict(tickmode='linear', tick0=0, dtick=1),
        yaxis=dict(showticklabels=False, range=[-1, 2])
    )
    
    return fig


def generate_fragment_map(
    peptide_sequence: str,
    fragment_ions: Dict[str, List[Dict[str, Any]]]
) -> go.Figure:
    """
    Generate a fragment map visualization.
    
    Args:
        peptide_sequence: Amino acid sequence
        fragment_ions: Dictionary of fragment ions
        
    Returns:
        Plotly figure object
    """
    fig = go.Figure()
    
    # Extract b-ions and y-ions information
    b_ions = fragment_ions.get("b_ions", [])
    y_ions = fragment_ions.get("y_ions", [])
    
    # Create lists to track which positions have ions
    b_positions = [ion["position"] for ion in b_ions]
    y_positions = [len(peptide_sequence) - ion["position"] for ion in y_ions]
    
    # Draw the peptide sequence
    for i, aa in enumerate(peptide_sequence):
        # Add amino acid letter
        fig.add_annotation(
            x=i,
            y=0,
            text=aa,
            showarrow=False,
            font=dict(size=16)
        )
        
        # Add b-ion marker if present
        if i+1 in b_positions:
            ion_info = next((ion for ion in b_ions if ion["position"] == i+1), None)
            mz_text = f"{ion_info['mz']:.2f}" if ion_info else ""
            
            fig.add_trace(go.Scatter(
                x=[i+0.5],
                y=[0.3],
                mode='markers+text',
                marker=dict(color="#1E88E5", size=10, symbol="triangle-down"),
                text=[f"b{i+1}"],
                textposition="top center",
                showlegend=i==0,
                name="b-ions",
                hoverinfo='text',
                hovertext=f"b{i+1}<br>m/z: {mz_text}"
            ))
        
        # Add y-ion marker if present
        if i+1 in y_positions:
            ion_index = len(peptide_sequence) - (i+1)
            ion_info = next((ion for ion in y_ions if ion["position"] == ion_index), None)
            mz_text = f"{ion_info['mz']:.2f}" if ion_info else ""
            
            fig.add_trace(go.Scatter(
                x=[i+0.5],
                y=[-0.3],
                mode='markers+text',
                marker=dict(color="#FFA000", size=10, symbol="triangle-up"),
                text=[f"y{ion_index}"],
                textposition="bottom center",
                showlegend=i==0,
                name="y-ions",
                hoverinfo='text',
                hovertext=f"y{ion_index}<br>m/z: {mz_text}"
            ))
    
    # Draw connecting lines for peptide backbone
    for i in range(len(peptide_sequence) - 1):
        fig.add_shape(
            type="line",
            x0=i+0.2, y0=0, x1=i+0.8, y1=0,
            line=dict(color="black", width=2)
        )
    
    fig.update_layout(
        title="Peptide Fragment Map",
        template='plotly_white',
        height=300,
        showlegend=True,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        xaxis=dict(showticklabels=False, range=[-0.5, len(peptide_sequence)-0.5]),
        yaxis=dict(showticklabels=False, range=[-1, 1]),
        margin=dict(l=40, r=40, t=60, b=40)
    )
    
    return fig


if __name__ == "__main__":
    # Example usage
    mz_values = [400.1978, 400.7, 401.2, 401.7, 402.2]
    intensities = [100, 80, 50, 20, 5]
    
    fig = plot_spectrum(mz_values, intensities, title="Example Spectrum")
    fig.show()

def create_interactive_fragment_map(
    peptide_sequence: str,
    fragment_ions: Dict[str, List[Dict[str, Any]]],
    title: str = "Interactive Fragment Map"
) -> go.Figure:
    """
    Create an interactive fragment ion map where users can click on peptide bonds 
    to highlight the resulting fragments.
    
    Args:
        peptide_sequence: Amino acid sequence
        fragment_ions: Dictionary of fragment ions
        title: Plot title
        
    Returns:
        Plotly figure object
    """
    fig = go.Figure()
    
    # Extract b-ions and y-ions information
    b_ions = {ion["position"]: ion for ion in fragment_ions.get("b_ions", [])}
    y_ions = {ion["position"]: ion for ion in fragment_ions.get("y_ions", [])}
    
    # Create a blank canvas for peptide visualization
    fig.add_trace(go.Scatter(
        x=[i for i in range(len(peptide_sequence))],
        y=[0] * len(peptide_sequence),
        mode='markers+text',
        marker=dict(size=30, color='rgba(0,0,0,0)'),  # Transparent markers
        text=list(peptide_sequence),
        textfont=dict(size=18),
        hoverinfo='none',
        showlegend=False
    ))
    
    # Add clickable bond areas
    for i in range(len(peptide_sequence) - 1):
        # Add invisible marker for each bond that will be clickable
        fig.add_trace(go.Scatter(
            x=[i + 0.5],
            y=[0],
            mode='markers',
            marker=dict(size=15, color='rgba(0,0,0,0.1)', symbol='square'),
            hoverinfo='text',
            hovertext=f"Click to see fragments for {peptide_sequence[:i+1]}-{peptide_sequence[i+1:]}",
            showlegend=False,
            customdata=[i + 1]  # Store the bond position for callback
        ))
    
    # Add bond lines
    for i in range(len(peptide_sequence) - 1):
        fig.add_shape(
            type="line",
            x0=i + 0.1, y0=0,
            x1=i + 0.9, y1=0,
            line=dict(color="black", width=2)
        )
    
    # Annotations for positions
    for i in range(len(peptide_sequence)):
        fig.add_annotation(
            x=i,
            y=0.2,
            text=f"{i+1}",
            showarrow=False,
            font=dict(size=10, color="gray")
        )
    
    # Create hidden traces for b-ions and y-ions that will be revealed on click
    for position in range(1, len(peptide_sequence)):
        # B ion
        if position in b_ions:
            ion = b_ions[position]
            fig.add_trace(go.Scatter(
                x=[position/2],
                y=[0.5],
                mode='markers+text',
                marker=dict(size=0),
                text=f"b{position} (m/z: {ion['mz']:.4f})",
                textfont=dict(color="#1E88E5"),
                visible=False,
                showlegend=False,
                customdata=[f"b{position}"]
            ))
        
        # Y ion
        y_pos = len(peptide_sequence) - position
        if position in y_ions:
            ion = y_ions[position]
            fig.add_trace(go.Scatter(
                x=[(len(peptide_sequence) + position)/2],
                y=[-0.5],
                mode='markers+text',
                marker=dict(size=0),
                text=f"y{position} (m/z: {ion['mz']:.4f})",
                textfont=dict(color="#FFA000"),
                visible=False,
                showlegend=False,
                customdata=[f"y{position}"]
            ))
    
    # Update layout
    fig.update_layout(
        title=title,
        height=300,
        showlegend=False,
        plot_bgcolor='white',
        xaxis=dict(
            showticklabels=False,
            range=[-0.5, len(peptide_sequence) - 0.5],
            zeroline=False
        ),
        yaxis=dict(
            showticklabels=False,
            range=[-1, 1],
            zeroline=False
        ),
        margin=dict(l=20, r=20, t=50, b=20),
        hovermode='closest'
    )
    
    return fig

def visualize_peptide_with_modifications(
    sequence: str,
    modifications: List[Dict[str, Any]] = None,
    width: int = 700,
    height: int = 500
):
    """
    Render a 3D visualization of the peptide with highlighted modifications.
    
    Args:
        sequence: Amino acid sequence
        modifications: List of modifications with position and name
        width: Width of the visualization
        height: Height of the visualization
    """
    try:
        import py3Dmol
        from stmol import showmol
        import requests
    except ImportError:
        st.error("Required packages not installed. Please install py3Dmol and stmol: pip install py3dmol stmol")
        return
    
    # Get PDB structure (simplified approach)
    def get_pdb_from_sequence(sequence):
        """Generate a simple PDB structure from the sequence"""
        pdb_content = "HEADER    PEPTIDE\n"
        pdb_content += f"TITLE     GENERATED MODEL OF {sequence}\n"
        
        # Add atoms for each amino acid in a linear chain
        for i, aa in enumerate(sequence):
            pdb_content += f"ATOM  {(i+1):5d}  CA  {aa}   A{i+1:4d}    {i*3.8:.3f}   0.000   0.000  1.00  0.00\n"
        
        # Connect atoms with bonds
        for i in range(len(sequence)-1):
            pdb_content += f"CONECT{i+1:5d}{i+2:5d}\n"
        
        pdb_content += "END\n"
        return pdb_content
    
    # Create a py3Dmol view
    view = py3Dmol.view(width=width, height=height)
    pdb_str = get_pdb_from_sequence(sequence)
    view.addModel(pdb_str, "pdb")
    
    # Basic style for the whole peptide
    view.setStyle({'cartoon': {'color': 'lightgray'}})
    
    # Add spheres for each alpha carbon
    view.addStyle({'atom': 'CA'}, {'sphere': {'radius': 0.6, 'color': 'lightgray'}})
    
    # Highlight modified residues if any
    if modifications:
        for mod in modifications:
            if mod.get('position') and mod.get('position') != 'N-term' and mod.get('position') != 'C-term':
                pos = int(mod['position'])
                # Select the residue
                selection = {'resi': pos}
                
                # Different colors for different modification types
                color = {
                    'Phospho': 'orange',
                    'Oxidation': 'red',
                    'Carbamidomethyl': 'green',
                    'Acetyl': 'blue',
                    'Methyl': 'purple',
                    'Dimethyl': 'magenta'
                }.get(mod.get('name'), 'yellow')
                
                # Style the modified residue
                view.addStyle(selection, {'cartoon': {'color': color}})
                view.addStyle({'resi': pos, 'atom': 'CA'}, 
                             {'sphere': {'radius': 1.0, 'color': color}})
                
                # Add a label
                view.addLabel(f"{mod.get('name')}", 
                             {'position': {'resi': pos, 'atom': 'CA'}, 
                              'backgroundColor': color,
                              'fontColor': 'white',
                              'fontSize': 12})
    
    # N-terminal and C-terminal modifications
    if modifications:
        for mod in modifications:
            if mod.get('position') == 'N-term':
                view.addStyle({'resi': 1}, {'cartoon': {'color': 'cyan'}})
                view.addStyle({'resi': 1, 'atom': 'CA'}, 
                             {'sphere': {'radius': 1.0, 'color': 'cyan'}})
                view.addLabel("N-term:" + mod.get('name', ''),
                             {'position': {'resi': 1, 'atom': 'CA'},
                              'backgroundColor': 'cyan',
                              'fontColor': 'black'})
            
            if mod.get('position') == 'C-term':
                view.addStyle({'resi': len(sequence)}, {'cartoon': {'color': 'pink'}})
                view.addStyle({'resi': len(sequence), 'atom': 'CA'}, 
                             {'sphere': {'radius': 1.0, 'color': 'pink'}})
                view.addLabel("C-term:" + mod.get('name', ''), 
                             {'position': {'resi': len(sequence), 'atom': 'CA'},
                              'backgroundColor': 'pink',
                              'fontColor': 'black'})
    
    # Final view settings
    view.zoomTo()
    view.setViewStyle({'style': 'outline', 'color': 'black', 'width': 0.1})
    view.setBackgroundColor('white')
    
    # Show in Streamlit
    showmol(view, height=height, width=width)