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
    .main-header {
        font-size: 2.5rem;
        font-weight: 700;
        color: #1E88E5;
        text-align: center;
        margin-bottom: 1rem;
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

# Logo and header
col1, col2, col3 = st.columns([1, 2, 1])
with col2:
    st.markdown('<div class="main-header">Peptide m/z Calculator</div>', unsafe_allow_html=True)
    st.markdown('*A professional tool for mass spectrometry analysis*', unsafe_allow_html=True)


# Function to calculate m/z from peptide sequence
def calculate_mz(sequence, charge_state=1, modifications=None, isotope_type="monoisotopic"):
    """
    Calculate the m/z ratio of a peptide sequence using PyOpenMS.
    
    Args:
        sequence (str): Amino acid sequence in single-letter code
        charge_state (int): Charge state of the peptide
        modifications (list): List of modifications to apply
        isotope_type (str): Type of isotope mass to use
        
    Returns:
        dict: Dictionary containing m/z value, mass, and other calculated properties
    """
    try:
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
        if isotope_type == "monoisotopic":
            mass = peptide.getMonoWeight()
        else:
            mass = peptide.getAverageWeight()
        
        mz = (mass + charge_state * 1.007276466) / charge_state  # Add proton mass for each charge
        
        # Get formula
        formula = peptide.getFormula()
        formula_string = str(formula)
        
        # Get fragment ions
        fragment_ions = {}
        # b-ions
        b_ions = []
        for i in range(1, len(sequence)):
            b_fragment = peptide.getPrefix(i)
            b_mass = b_fragment.getMonoWeight(pyopenms.Residue.ResidueType.BIon, charge_state)
            b_ions.append({"position": i, "mass": b_mass, "mz": b_mass/charge_state})
        fragment_ions["b_ions"] = b_ions
        
        # y-ions
        y_ions = []
        for i in range(1, len(sequence)):
            y_fragment = peptide.getSuffix(i)
            y_mass = y_fragment.getMonoWeight(pyopenms.Residue.ResidueType.YIon, charge_state)
            y_ions.append({"position": i, "mass": y_mass, "mz": y_mass/charge_state})
        fragment_ions["y_ions"] = y_ions
        
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
        st.error(f"Error calculating m/z: {str(e)}")
        return None


# Function to generate theoretical isotope pattern
def generate_isotope_pattern(peptide_sequence, charge_state=1):
    """
    Generate theoretical isotope pattern for a peptide sequence.
    
    Args:
        peptide_sequence (str): Amino acid sequence
        charge_state (int): Charge state of the peptide
        
    Returns:
        tuple: (mz_values, intensities) arrays for plotting
    """
    try:
        peptide = AASequence.fromString(peptide_sequence)
        formula = peptide.getFormula()
        
        isotope_dist = pyopenms.CoarseIsotopePatternGenerator(3)
        iso_pattern = pyopenms.IsotopeDistribution()
        isotope_dist.estimateFromPeptideWeight(peptide.getMonoWeight())
        
        masses = []
        intensities = []
        
        for i in range(iso_pattern.size()):
            masses.append(iso_pattern.getMassAbsolute(i) / charge_state)
            intensities.append(iso_pattern.getIntensity(i))
            
        # If empty, use theoretical approximation
        if not masses:
            mono_mass = peptide.getMonoWeight(pyopenms.Residue.ResidueType.Full, charge_state)
            mz = mono_mass / charge_state
            
            # Approximate isotope pattern
            masses = [mz + (i * 1.003) / charge_state for i in range(5)]
            # Typical intensity distribution for peptides
            intensities = [100, 80 - (charge_state * 10), 60 - (charge_state * 15), 
                           30 - (charge_state * 5), 10]
            intensities = [max(i, 0) for i in intensities]
        
        return masses, intensities
    
    except Exception as e:
        st.error(f"Error generating isotope pattern: {str(e)}")
        return [], []


# Function to parse FASTA file
def parse_fasta(fasta_content):
    """
    Parse a FASTA file and extract sequences.
    
    Args:
        fasta_content (str): Content of the FASTA file
        
    Returns:
        list: List of dictionaries with header and sequence
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


# Function to plot mass spectrum
def plot_spectrum(mz_values, intensities, title="Theoretical Mass Spectrum"):
    """
    Create an interactive plot of a mass spectrum.
    
    Args:
        mz_values (list): m/z values
        intensities (list): Corresponding intensity values
        title (str): Plot title
        
    Returns:
        plotly.graph_objects.Figure: Plotly figure object
    """
    fig = go.Figure()
    
    fig.add_trace(go.Scatter(
        x=mz_values,
        y=intensities,
        mode='lines',
        name='Isotope Pattern',
        line=dict(color='#1E88E5', width=2)
    ))
    
    fig.update_layout(
        title=title,
        xaxis_title='m/z',
        yaxis_title='Relative Intensity (%)',
        template='plotly_white',
        height=500,
        margin=dict(l=0, r=0, t=40, b=0),
        showlegend=False
    )
    
    return fig


# Function to calculate amino acid composition
def get_aa_composition(sequence):
    """
    Calculate the amino acid composition of a peptide.
    
    Args:
        sequence (str): Amino acid sequence
        
    Returns:
        dict: Dictionary of amino acids and their counts
    """
    composition = {}
    for aa in sequence:
        if aa in composition:
            composition[aa] += 1
        else:
            composition[aa] = 1
    
    return composition


# Function to export results as CSV
def export_to_csv(results):
    """
    Convert results to CSV format.
    
    Args:
        results (dict): Results from calculation
        
    Returns:
        str: CSV content as string
    """
    output = io.StringIO()
    writer = pd.DataFrame([{
        'Sequence': results['sequence'],
        'Charge State': results['charge'],
        'Mass (Da)': results['mass'],
        'm/z': results['mz'],
        'Formula': results['formula'],
        'Isotope Type': results['isotope_type']
    }])
    
    writer.to_csv(output, index=False)
    return output.getvalue()


# Sidebar for application navigation
with st.sidebar:
    st.image("https://www.openms.de/wp-content/uploads/2016/06/OpenMS.png", width=200)
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
    
    modifications = []
    if add_modifications:
        mod_expander = st.expander("Modification Panel", expanded=True)
        
        with mod_expander:
            # Common modifications
            common_mods = {
                "None": None,
                "Phosphorylation (S,T,Y)": "Phospho",
                "Oxidation (M)": "Oxidation",
                "Carbamidomethyl (C)": "Carbamidomethyl",
                "Acetylation (K, N-term)": "Acetyl",
                "Methylation (K, R)": "Methyl",
                "Dimethylation (K, R)": "Dimethyl"
            }
            
            selected_mod = st.selectbox(
                "Common Modifications",
                options=list(common_mods.keys()),
                index=0
            )
            
            if selected_mod != "None":
                mod_name = common_mods[selected_mod]
                
                # Find applicable positions based on the modification
                applicable_positions = []
                
                if "N-term" in selected_mod:
                    applicable_positions.append("N-term")
                
                if "C-term" in selected_mod:
                    applicable_positions.append("C-term")
                
                for i, aa in enumerate(sequence):
                    if ("(S,T,Y)" in selected_mod and aa in "STY") or \
                       ("(M)" in selected_mod and aa == "M") or \
                       ("(C)" in selected_mod and aa == "C") or \
                       ("(K" in selected_mod and aa == "K") or \
                       ("(K, R)" in selected_mod and aa in "KR"):
                        applicable_positions.append(str(i + 1))
                
                if applicable_positions:
                    mod_position = st.selectbox(
                        "Position",
                        options=applicable_positions
                    )
                    
                    add_mod_button = st.button("Add Modification")
                    
                    if add_mod_button:
                        modifications.append({
                            "name": mod_name,
                            "position": mod_position
                        })
            
            # Display current modifications
            if modifications:
                st.markdown("**Current Modifications:**")
                for i, mod in enumerate(modifications):
                    st.markdown(f"{i+1}. {mod['name']} at position {mod['position']}")
                
                if st.button("Clear All Modifications"):
                    modifications = []
    
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
                    result = calculate_mz(cleaned_seq, batch_charge, [], isotope_type)
                    
                    if result:
                        results.append({
                            "Sequence": cleaned_seq,
                            "Mass (Da)": result["mass"],
                            "m/z": result["mz"],
                            "Charge": f"+{batch_charge}",
                            "Formula": result["formula"]
                        })
            
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
          - Interactive plots
          - Export options
        """)
    
    with col2:
        st.markdown("""
        - **Batch Processing**
          - Multiple sequence handling
          - Bulk calculations
          - Results export
        
        - **Scientific Accuracy**
          - Powered by PyOpenMS
          - Industry-standard calculations
          - Detailed results
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

# Footer
st.markdown("""
<div class="footer">
    <p>Peptide m/z Calculator | Developed for OpenMS | Google Summer of Code</p>
    <p>Powered by PyOpenMS and Streamlit</p>
</div>
""", unsafe_allow_html=True)