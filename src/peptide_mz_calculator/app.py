import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pyopenms
from pyopenms import AASequence, EmpiricalFormula
import plotly.graph_objects as go
import io
import base64
from PIL import Image
import tempfile
import os
import threading
import socket
import json

# Import from other modules
from calculator import calculate_mz, generate_isotope_pattern, generate_fragment_ions
from parser import parse_fasta, parse_file, detect_file_format
from export import export_to_csv, export_to_excel, export_summary_statistics
from utils import get_aa_composition, validate_sequence, format_formula
from visualization import (
    plot_spectrum, plot_fragment_ions, plot_aa_composition, create_interactive_fragment_map,
    visualize_peptide_with_modifications
)
from modifications import (
    load_modifications, render_modification_selector, visualize_modified_peptide,
    render_modification_set_selector, register_custom_modification_with_pyopenms
)
from openms_integration import OpenMSIntegration
import api as api

def get_base64_encoded_image(image_path):
    """
    Returns the base64 encoded string of an image
    """
    with open(image_path, "rb") as img_file:
        return base64.b64encode(img_file.read()).decode()
    
# Set page configuration
st.set_page_config(
    page_title="Peptide m/z Calculator",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for professional look
st.markdown("""
<style>
    .main-head{
        text-align: center;
    }
    .main-header {
        font-size: 2.5rem;
        font-weight: 700;
        text-align: center;
        margin-bottom: 1rem;
        padding: 10px;
        background-clip: text;
        -webkit-background-clip: text;
        color: transparent;
        background-image: linear-gradient(45deg, 
            #FFAE01, #FF9601, #FF8600, #FF7901, 
            #FF6124, #FF4B4C, #FF3B6A, #FF16A9, 
            #FF05C8, #EC02D9, #C903F1, #B703FD, 
            #772AF6, #414DEC, #3157E9);
        background-size: 300% 300%;
        animation: gradientAnimation 10s ease infinite;
        /* Add text shadow to improve readability */
        text-shadow: 0px 0px 1px rgba(255,255,255,0.1);
    }
    @keyframes gradientAnimation {
        0% {
            background-position: 0% 50%;
        }
        50% {
            background-position: 100% 50%;
        }
        100% {
            background-position: 0% 50%;
        }
    }
    .sub-header {
        font-size: 1.5rem;
        font-weight: 600;
        color: #0D47A1;
        margin-top: 1.5rem;
        margin-bottom: 1rem;
    }
    .info-box {
        background-color: #E3F2FD;
        padding: 1rem;
        border-radius: 0.5rem;
        border-left: 0.5rem solid #1E88E5;
        margin-bottom: 1rem;
    }
    .result-box {
        background-color: #E8F5E9;
        padding: 1.5rem;
        border-radius: 0.5rem;
        border-left: 0.5rem solid #43A047;
        margin-top: 1rem;
        margin-bottom: 1rem;
    }
    .warning-box {
        background-color: #FFF8E1;
        padding: 1rem;
        border-radius: 0.5rem;
        border-left: 0.5rem solid #FFA000;
        margin-bottom: 1rem;
    }
    .footer {
        text-align: center;
        color: #616161;
        font-size: 0.8rem;
        margin-top: 3rem;
        padding-top: 1rem;
        border-top: 1px solid #EEEEEE;
    }
    /* Custom select widget styles */
    div[data-baseweb="select"] {
        margin-top: 0.5rem;
        margin-bottom: 1.5rem;
    }
</style>
""", unsafe_allow_html=True)

st.markdown("""
<style>
    /* Reset Streamlit's defaults */
    .stApp {
        margin: 0;
        padding: 0;
    }
    
    /* Remove padding from container */
    .block-container {
        padding-top: 0 !important;
    }
    
    /* Custom header positioning */
    .custom-header {
        position: absolute;
        top: 1cm;  /* Adjust this value - smaller = higher */
        left: 0;
        right: 0;
        text-align: center;
        z-index: 1000;
    }
    
    /* Add padding to the main content to prevent overlap */
    .main-content-padding {
        padding-top: 5cm;
    }
</style>

<div class="custom-header">
""", unsafe_allow_html=True)

# Logo and header
st.markdown("""
<div style='display: flex; flex-direction: column; align-items: center; justify-content: center; text-align: center; margin: 0; padding: 0;'>
""", unsafe_allow_html=True)

# Load and display the logo centered and smaller
try:
    logo_path = os.path.join(os.path.dirname(__file__), "assets", "logo.png")
    if os.path.exists(logo_path):
        st.markdown(f"""
        <div style='display: flex; justify-content: center; margin-bottom: -10px;'>
            <img src="data:image/png;base64,{get_base64_encoded_image(logo_path)}" width="100px">
        </div>
        <div class="main-header" style="margin-top: -5px;">Peptide m/z Calculator</div>
        <p style="font-style: italic; text-align: center; margin-top: 0;">A professional tool for mass spectrometry analysis</p>
        """, unsafe_allow_html=True)
    else:
        st.warning("Logo file not found. Please place your logo at: " + logo_path)
except Exception as e:
    st.error(f"Error loading logo: {str(e)}")

st.markdown("</div>", unsafe_allow_html=True)
st.markdown('<div class="main-content-padding"></div>', unsafe_allow_html=True)
# Sidebar for application navigation
with st.sidebar:
    # First add your custom logo
    try:
        logo_path = os.path.join(os.path.dirname(__file__), "assets", "logo.png")
        if os.path.exists(logo_path):
            st.image(logo_path, width=150)
        else:
            st.warning("Logo file not found. Please place your logo at: " + logo_path)
    except Exception as e:
        st.error(f"Error loading logo: {str(e)}")
    st.markdown("## Navigation")
    
    app_mode = st.radio(
        "Select Mode",
        options=["Calculator", "Batch Processing", "About"],
        index=0
    )
    
    st.markdown("---")
    st.markdown("## Settings")
    
    mass_type = st.selectbox(
        "Mass Type",
        options=["Monoisotopic", "Average"],
        index=0
    )
    
    st.markdown("### Advanced Settings")
    
    ppm_tolerance = st.slider(
        "Mass Tolerance (ppm)",
        min_value=1,
        max_value=100,
        value=10
    )
    
    min_intensity = st.slider(
        "Minimum Intensity (%)",
        min_value=1,
        max_value=50,
        value=5
    )
    
    st.markdown("### Advanced Features")
    
    enable_3d = st.checkbox(
        "Enable 3D Visualization",
        value=False,
        help="Show 3D visualization of peptides"
    )
    
    enable_api = st.checkbox(
        "Enable API Server",
        value=False,
        help="Start the API server alongside the UI"
    )
    
    st.markdown("### Integration")
    
    integration_options = st.multiselect(
        "Integration Options",
        options=["Export to OpenMS TraML", "Predict Retention Time", "Match Against Spectral Library"],
        default=[]
    )


# Main application logic
if app_mode == "Calculator":
    st.markdown('<div class="sub-header">Peptide m/z Calculator</div>', unsafe_allow_html=True)
    
    # Input method selection
    input_method = st.radio(
        "Input Method",
        options=["Sequence Entry", "FASTA Upload", "Preset Peptides"],
        horizontal=True
    )
    
    # Input based on selected method
    sequence = ""
    if input_method == "Sequence Entry":
        sequence = st.text_input(
            "Enter Peptide Sequence (single-letter amino acid code)",
            value="PEPTIDE",
            help="Example: ACDEFGHIKLMNPQRSTVWY"
        )
        
    elif input_method == "FASTA Upload":
        fasta_file = st.file_uploader("Upload FASTA File", type=['fasta', 'fa', 'txt'])
        
        if fasta_file is not None:
            fasta_content = fasta_file.read().decode('utf-8')
            sequences = parse_fasta(fasta_content)
            
            if sequences:
                sequence_options = [f"{seq['header']} | {seq['sequence'][:20]}..." 
                                   if len(seq['sequence']) > 20 else f"{seq['header']} | {seq['sequence']}" 
                                   for seq in sequences]
                
                selected_sequence = st.selectbox(
                    "Select Sequence",
                    options=sequence_options
                )
                
                selected_idx = sequence_options.index(selected_sequence)
                sequence = sequences[selected_idx]['sequence']
                
                st.markdown(f"<div class='info-box'>Selected sequence: {sequence}</div>", unsafe_allow_html=True)
            else:
                st.warning("No sequences found in the uploaded file.")
                
    elif input_method == "Preset Peptides":
        preset_peptides = {
            "Angiotensin I": "DRVYIHPFHL",
            "Bradykinin": "RPPGFSPFR",
            "Substance P": "RPKPQQFFGLM",
            "ACTH (1-10)": "SYSMEHFRWG",
            "Insulin Chain B Oxidized": "FVNQHLCGSHLVEALYLVCGERGFFYTPKA",
            "Bombesin": "EQRLGNQWAVGHLM"
        }
        
        selected_preset = st.selectbox(
            "Select Preset Peptide",
            options=list(preset_peptides.keys())
        )
        
        sequence = preset_peptides[selected_preset]
        st.markdown(f"<div class='info-box'>Sequence: {sequence}</div>", unsafe_allow_html=True)
    
    # Charge state selection
    charge_state = st.selectbox(
        "Charge State",
        options=[1, 2, 3, 4, 5, 6],
        index=1,  # Default to +2
        help="Select the charge state (z) for the peptide"
    )
    
    # Modifications
    st.markdown('<div class="sub-header">Modifications</div>', unsafe_allow_html=True)
    
    add_modifications = st.checkbox("Add Modifications", value=False)
    
    # Initialize modifications list
    modifications = []
    
    # Load modifications
    modifications_data = load_modifications("data/modifications.json")
    
    # Expander for predefined modification sets
    with st.expander("Predefined Modification Sets", expanded=False):
        preset_mods = render_modification_set_selector(sequence, modifications_data)
        
        if preset_mods and st.button("Apply This Set"):
            modifications = preset_mods
            st.success(f"Applied modification set with {len(modifications)} modifications.")
    
    # Expander for custom modifications
    with st.expander("Modification Editor", expanded=add_modifications):
        modifications = render_modification_selector(sequence, modifications_data, modifications)
    
    # Visualize the modified peptide
    if sequence and modifications:
        visualize_modified_peptide(sequence, modifications)
    
    # 3D Visualization if enabled
    if enable_3d and sequence:
        st.markdown('<div class="sub-header">3D Visualization</div>', unsafe_allow_html=True)
        visualize_peptide_with_modifications(sequence, modifications)
    
    # Calculate button
    st.markdown("---")
    calculate_button = st.button("Calculate m/z", type="primary")
    
    if calculate_button and sequence:
        # Preprocess the sequence (remove spaces and convert to uppercase)
        sequence = sequence.replace(" ", "").upper()
        
        # Check for invalid amino acids
        valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
        invalid_aa = [aa for aa in sequence if aa not in valid_aa]
        
        if invalid_aa:
            st.warning(f"Warning: Sequence contains invalid amino acids: {', '.join(invalid_aa)}")
        else:
            # Perform the calculation
            isotope_type = "monoisotopic" if mass_type == "Monoisotopic" else "average"
            result = calculate_mz(sequence, charge_state, modifications, isotope_type)
            
            if result:
                # Display results
                st.markdown('<div class="result-box">', unsafe_allow_html=True)
                
                col1, col2 = st.columns([1, 1])
                
                with col1:
                    st.markdown(f"**Peptide Sequence:** {result['sequence']}")
                    st.markdown(f"**Charge State:** +{result['charge']}")
                    st.markdown(f"**Mass Type:** {result['isotope_type'].capitalize()}")
                    
                with col2:
                    st.markdown(f"**m/z Ratio:** {result['mz']:.4f}")
                    st.markdown(f"**Mass:** {result['mass']:.4f} Da")
                    st.markdown(f"**Chemical Formula:** {result['formula']}")
                
                st.markdown('</div>', unsafe_allow_html=True)
                
                # Generate isotope pattern
                mz_values, intensities = generate_isotope_pattern(sequence, charge_state)
                
                if mz_values and intensities:
                    st.markdown('<div class="sub-header">Theoretical Isotope Pattern</div>', unsafe_allow_html=True)
                    
                    # Normalize intensities
                    max_intensity = max(intensities)
                    normalized_intensities = [i * 100 / max_intensity for i in intensities]
                    
                    # Create plot
                    spectrum_fig = plot_spectrum(
                        mz_values, 
                        normalized_intensities,
                        f"Theoretical Isotope Pattern for {sequence} (Charge: +{charge_state})"
                    )
                    st.plotly_chart(spectrum_fig, use_container_width=True)
                
                # Add interactive fragment map
                st.markdown('<div class="sub-header">Interactive Fragment Map</div>', unsafe_allow_html=True)
                
                fragment_map_fig = create_interactive_fragment_map(
                    sequence, 
                    result['fragment_ions'],
                    f"Fragment Map for {sequence} (Charge: +{charge_state})"
                )
                st.plotly_chart(fragment_map_fig, use_container_width=True)
                
                # Show fragment ions
                with st.expander("Fragment Ions", expanded=False):
                    st.markdown("### b-ions")
                    b_ions_df = pd.DataFrame(result['fragment_ions']['b_ions'])
                    st.dataframe(b_ions_df)
                    
                    st.markdown("### y-ions")
                    y_ions_df = pd.DataFrame(result['fragment_ions']['y_ions'])
                    st.dataframe(y_ions_df)
                
                # Amino acid composition
                with st.expander("Amino Acid Composition", expanded=False):
                    composition = get_aa_composition(sequence)
                    comp_df = pd.DataFrame({
                        'Amino Acid': list(composition.keys()),
                        'Count': list(composition.values())
                    })
                    st.dataframe(comp_df)
                
                # OpenMS integration features
                if integration_options:
                    st.markdown('<div class="sub-header">OpenMS Integration</div>', unsafe_allow_html=True)
                    
                    if "Export to OpenMS TraML" in integration_options:
                        # Create TraML export
                        try:
                            traml_file = OpenMSIntegration.export_to_traml(
                                [{
                                    "sequence": sequence,
                                    "charge": charge_state,
                                    "mz": result["mz"],
                                    "fragments": [
                                        {"mz": ion["mz"], "type": f"b{ion['position']}", "charge": 1}
                                        for ion in result["fragment_ions"]["b_ions"]
                                    ] + [
                                        {"mz": ion["mz"], "type": f"y{ion['position']}", "charge": 1}
                                        for ion in result["fragment_ions"]["y_ions"]
                                    ]
                                }]
                            )
                            
                            with open(traml_file, 'r') as f:
                                traml_content = f.read()
                                
                            st.download_button(
                                "Download TraML File",
                                data=traml_content,
                                file_name=f"{sequence}_transitions.traML",
                                mime="text/xml"
                            )
                        except Exception as e:
                            st.error(f"Error creating TraML file: {str(e)}")
                    
                    if "Predict Retention Time" in integration_options:
                        try:
                            rt_results = OpenMSIntegration.predict_retention_time([sequence])
                            
                            if sequence in rt_results:
                                st.info(f"Predicted Retention Time: {rt_results[sequence]:.2f} minutes")
                            else:
                                st.warning("Unable to predict retention time for this sequence.")
                        except Exception as e:
                            st.error(f"Error predicting retention time: {str(e)}")
                    
                    if "Match Against Spectral Library" in integration_options:
                        st.markdown("#### Spectral Library Matching")
                        
                        # Upload spectral library
                        spectral_lib = st.file_uploader(
                            "Upload Spectral Library (MSP/MGF format)",
                            type=["msp", "mgf"]
                        )
                        
                        if spectral_lib is not None:
                            # Save to temp file
                            with tempfile.NamedTemporaryFile(delete=False, suffix=f".{spectral_lib.name.split('.')[-1]}") as tmp:
                                tmp.write(spectral_lib.getvalue())
                                lib_path = tmp.name
                            
                            # Generate theoretical spectrum
                            mz_values, intensities, _ = OpenMSIntegration.generate_theoretical_spectrum(
                                sequence, charge_state
                            )
                            
                            # Match against library
                            if st.button("Search Library"):
                                try:
                                    matches = OpenMSIntegration.match_against_spectral_library(
                                        mz_values, intensities, lib_path, result["mz"], charge_state
                                    )
                                    
                                    if matches:
                                        st.success(f"Found {len(matches)} matches in the spectral library.")
                                        
                                        # Display top matches
                                        match_df = pd.DataFrame(matches[:5])  # Top 5 matches
                                        st.dataframe(match_df)
                                    else:
                                        st.info("No matches found in the spectral library.")
                                except Exception as e:
                                    st.error(f"Error matching against spectral library: {str(e)}")
                
                # Export options
                st.markdown('<div class="sub-header">Export Results</div>', unsafe_allow_html=True)
                
                col1, col2 = st.columns([1, 1])
                
                with col1:
                    # Export as CSV
                    csv_data = export_to_csv(result)
                    st.download_button(
                        label="Download Results (CSV)",
                        data=csv_data,
                        file_name=f"{sequence}_mz_{result['charge']}+.csv",
                        mime="text/csv"
                    )
                
                with col2:
                    # Export plot as image
                    if mz_values and intensities:
                        # Convert Plotly figure to image
                        img_bytes = spectrum_fig.to_image(format="png", scale=3)
                        st.download_button(
                            label="Download Spectrum (PNG)",
                            data=img_bytes,
                            file_name=f"{sequence}_spectrum_{result['charge']}+.png",
                            mime="image/png"
                        )
    
    elif calculate_button and not sequence:
        st.warning("Please enter a peptide sequence or select from presets.")

elif app_mode == "Batch Processing":
    st.markdown('<div class="sub-header">Batch Processing</div>', unsafe_allow_html=True)
    
    st.markdown("""
    <div class="info-box">
        Upload a file containing multiple peptide sequences for batch processing.
        Supported formats: CSV, TSV, or TXT with one sequence per line.
    </div>
    """, unsafe_allow_html=True)
    
    # File upload
    upload_method = st.radio(
        "Input Method",
        options=["File Upload", "Multiple Sequence Entry"],
        horizontal=True
    )
    
    sequences = []
    
    if upload_method == "File Upload":
        batch_file = st.file_uploader(
            "Upload Sequence File",
            type=['csv', 'tsv', 'txt', 'fasta', 'fa']
        )
        
        if batch_file is not None:
            file_extension = batch_file.name.split('.')[-1].lower()
            
            if file_extension in ['fasta', 'fa']:
                # Parse FASTA file
                fasta_content = batch_file.read().decode('utf-8')
                parsed_fasta = parse_fasta(fasta_content)
                
                if parsed_fasta:
                    sequences = [seq['sequence'] for seq in parsed_fasta]
                    headers = [seq['header'] for seq in parsed_fasta]
                    
                    st.success(f"Successfully loaded {len(sequences)} sequences from FASTA file.")
                else:
                    st.warning("No sequences found in the uploaded FASTA file.")
            
            elif file_extension in ['csv', 'tsv']:
                # Parse CSV/TSV file
                separator = ',' if file_extension == 'csv' else '\t'
                
                try:
                    df = pd.read_csv(batch_file, sep=separator)
                    
                    # Try to find sequence column
                    possible_column_names = ['Sequence', 'sequence', 'SEQUENCE', 'peptide', 'Peptide', 'PEPTIDE']
                    found_column = None
                    
                    for col_name in possible_column_names:
                        if col_name in df.columns:
                            found_column = col_name
                            break
                    
                    if found_column:
                        sequences = df[found_column].tolist()
                    else:
                        st.warning("Could not find a sequence column. Using the first column.")
                        sequences = df.iloc[:, 0].tolist()
                    
                    st.success(f"Successfully loaded {len(sequences)} sequences from file.")
                
                except Exception as e:
                    st.error(f"Error reading file: {str(e)}")
            
            else:  # Plain text file
                content = batch_file.read().decode('utf-8')
                sequences = [line.strip() for line in content.split('\n') if line.strip()]
                
                st.success(f"Successfully loaded {len(sequences)} sequences from text file.")
    
    else:  # Multiple Sequence Entry
        sequence_text = st.text_area(
            "Enter multiple sequences (one per line)",
            height=150,
            help="Enter one peptide sequence per line"
        )
        
        if sequence_text:
            sequences = [line.strip() for line in sequence_text.split('\n') if line.strip()]
    
    # Batch processing settings
    if sequences:
        st.markdown('<div class="sub-header">Batch Processing Settings</div>', unsafe_allow_html=True)
        
        col1, col2 = st.columns([1, 1])
        
        with col1:
            batch_charge = st.selectbox(
                "Charge State",
                options=[1, 2, 3, 4, 5, 6],
                index=1
            )
        
        with col2:
            batch_mass_type = st.selectbox(
                "Mass Type",
                options=["Monoisotopic", "Average"],
                index=0
            )
        
        # Process button
        process_button = st.button("Process Sequences", type="primary")
        
        if process_button:
            # Show progress bar
            progress_bar = st.progress(0)
            status_text = st.empty()
            
            # List to store results
            results = []
            
            # Process each sequence
            for i, seq in enumerate(sequences):
                status_text.text(f"Processing sequence {i+1} of {len(sequences)}: {seq}")
                progress_bar.progress((i + 1) / len(sequences))
                
                # Clean sequence
                cleaned_seq = seq.replace(" ", "").upper()
                
                # Check for valid sequence
                valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
                if all(aa in valid_aa for aa in cleaned_seq):
                    # Calculate m/z
                    isotope_type = "monoisotopic" if batch_mass_type == "Monoisotopic" else "average"
                    try:
                        result = calculate_mz(cleaned_seq, batch_charge, [], isotope_type)
                        
                        if result:
                            results.append({
                                "Sequence": cleaned_seq,
                                "Mass (Da)": result["mass"],
                                "m/z": result["mz"],
                                "Charge": f"+{batch_charge}",
                                "Formula": result["formula"]
                            })
                    except Exception as e:
                        st.error(f"Error processing sequence {cleaned_seq}: {str(e)}")
            
            # Display results
            if results:
                st.markdown('<div class="sub-header">Results</div>', unsafe_allow_html=True)
                
                # Convert to DataFrame
                results_df = pd.DataFrame(results)
                st.dataframe(results_df)
                
                # Export options
                csv_data = results_df.to_csv(index=False)
                st.download_button(
                    label="Download Results (CSV)",
                    data=csv_data,
                    file_name="batch_results.csv",
                    mime="text/csv"
                )
                
                # Excel export
                excel_buffer = io.BytesIO()
                results_df.to_excel(excel_buffer, index=False, engine='openpyxl')
                excel_data = excel_buffer.getvalue()
                
                st.download_button(
                    label="Download Results (Excel)",
                    data=excel_data,
                    file_name="batch_results.xlsx",
                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                )
                
                # Summary statistics
                st.markdown('<div class="sub-header">Summary Statistics</div>', unsafe_allow_html=True)
                
                stats_df = pd.DataFrame({
                    "Statistic": ["Count", "Mean Mass", "Median Mass", "Min Mass", "Max Mass", 
                                 "Mean m/z", "Median m/z", "Min m/z", "Max m/z"],
                    "Value": [
                        len(results),
                        f"{results_df['Mass (Da)'].mean():.4f} Da",
                        f"{results_df['Mass (Da)'].median():.4f} Da",
                        f"{results_df['Mass (Da)'].min():.4f} Da",
                        f"{results_df['Mass (Da)'].max():.4f} Da",
                        f"{results_df['m/z'].mean():.4f}",
                        f"{results_df['m/z'].median():.4f}",
                        f"{results_df['m/z'].min():.4f}",
                        f"{results_df['m/z'].max():.4f}"
                    ]
                })
                
                st.table(stats_df)
                
                # Mass distribution plot
                st.markdown('<div class="sub-header">Mass Distribution</div>', unsafe_allow_html=True)
                
                fig = go.Figure()
                
                fig.add_trace(go.Histogram(
                    x=results_df['Mass (Da)'],
                    nbinsx=20,
                    marker_color='#1E88E5',
                    opacity=0.7
                ))
                
                fig.update_layout(
                    title="Distribution of Peptide Masses",
                    xaxis_title="Mass (Da)",
                    yaxis_title="Count",
                    template="plotly_white"
                )
                
                st.plotly_chart(fig, use_container_width=True)

elif app_mode == "About":
    st.markdown('<div class="sub-header">About Peptide m/z Calculator</div>', unsafe_allow_html=True)
    
    st.markdown("""
    <div class="info-box">
        <h3>Peptide m/z Calculator</h3>
        <p>A professional tool for mass spectrometry researchers to calculate mass-to-charge ratios of peptides and other biomolecules without requiring programming knowledge.</p>
        <p>Built on the OpenMS framework and deployed using Streamlit, this tool bridges the gap between wet lab scientists and computational proteomics.</p>
    </div>
    """, unsafe_allow_html=True)
    
    st.markdown('<div class="sub-header">Features</div>', unsafe_allow_html=True)
    
    col1, col2 = st.columns([1, 1])
    
    with col1:
        st.markdown("""
        - **Peptide m/z Calculation**
          - Multiple input methods
          - Charge state options
          - Modification support
          - Isotope options
        
        - **Visualization**
          - Theoretical spectrum
          - Interactive fragment maps
          - 3D peptide visualization
          - Interactive plots
          - Export options
        """)
    
    with col2:
        st.markdown("""
        - **Batch Processing**
          - Multiple sequence handling
          - Bulk calculations
          - Results export
        
        - **Integration Features**
          - REST API for programmatic access
          - OpenMS integration (TraML, RT prediction)
          - Spectral library matching
          - Advanced modification handling
        """)
    
    st.markdown('<div class="sub-header">How It Works</div>', unsafe_allow_html=True)
    
    st.markdown("""
    <div class="info-box">
        <h4>Calculation Methodology</h4>
        <p>The calculator uses PyOpenMS to compute m/z values based on the amino acid sequence:</p>
        <ol>
            <li>The peptide sequence is parsed and validated</li>
            <li>Each amino acid's mass is calculated using standard values</li>
            <li>Modifications are applied if specified</li>
            <li>The mass is divided by the charge state to obtain m/z</li>
            <li>Theoretical isotope patterns are generated</li>
            <li>Fragment ions are calculated for structural analysis</li>
        </ol>
    </div>
    """, unsafe_allow_html=True)
    
    st.markdown('<div class="sub-header">Technical Information</div>', unsafe_allow_html=True)
    
    st.markdown("""
    <div class="info-box">
        <h4>Technical Details</h4>
        <p>This application is built using:</p>
        <ul>
            <li><strong>PyOpenMS</strong>: Core calculation engine</li>
            <li><strong>Streamlit</strong>: Web interface</li>
            <li><strong>Pandas</strong>: Data handling</li>
            <li><strong>Plotly</strong>: Interactive visualizations</li>
            <li><strong>py3Dmol/stmol</strong>: 3D peptide visualization</li>
            <li><strong>FastAPI</strong>: REST API functionality</li>
        </ul>
        <p>For more information, see the <a href="https://github.com/OpenMS/OpenMS" target="_blank">OpenMS GitHub repository</a>.</p>
    </div>
    """, unsafe_allow_html=True)
    
    st.markdown('<div class="sub-header">Project Information</div>', unsafe_allow_html=True)
    
    st.markdown("""
    <div class="info-box">
        <h4>GSoC Project</h4>
        <p>This application was developed as part of the Google Summer of Code (GSoC) program for the OpenMS organization.</p>
        <p>The goal was to create a user-friendly interface for MS proteomics calculations to help wet lab scientists without programming experience.</p>
    </div>
    """, unsafe_allow_html=True)
    
    # Add information about the API
    if enable_api:
        st.markdown('<div class="sub-header">API Documentation</div>', unsafe_allow_html=True)
        
        st.markdown("""
        <div class="info-box">
            <h4>REST API</h4>
            <p>The Peptide m/z Calculator provides a REST API for programmatic access to its functionality. 
            When enabled, the API server runs alongside the web interface and provides the following endpoints:</p>
            <ul>
                <li><strong>/calculate</strong>: Calculate m/z for a single peptide</li>
                <li><strong>/batch-calculate</strong>: Process multiple peptides</li>
                <li><strong>/parse-fasta</strong>: Parse FASTA format</li>
            </ul>
            <p>API documentation is available at the /docs endpoint when the API server is running.</p>
        </div>
        """, unsafe_allow_html=True)

# Footer
st.markdown("""
<div class="footer">
    <p>Peptide m/z Calculator | Developed for OpenMS | Google Summer of Code</p>
    <p>Powered by PyOpenMS and Streamlit</p>
</div>
""", unsafe_allow_html=True)

# Start API server if enabled
if enable_api:
    def is_port_in_use(port):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            return s.connect_ex(('localhost', port)) == 0
    
    # Find an available port
    api_port = 8000
    while is_port_in_use(api_port):
        api_port += 1
    
    # Start API server in a separate thread
    def run_api():
        api.run_api_server(host="0.0.0.0", port=api_port)
    
    api_thread = threading.Thread(target=run_api, daemon=True)
    api_thread.start()
    
    st.sidebar.success(f"API server running at http://localhost:{api_port}")
    st.sidebar.markdown(f"Documentation available at http://localhost:{api_port}/docs")