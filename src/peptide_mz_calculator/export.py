"""
Export module for generating output files.

This module provides functions for exporting calculation results in
various formats like CSV, Excel, and image files.
"""

import pandas as pd
import io
import base64
from datetime import datetime
import json
from typing import Dict, List, Any, Union, Optional


def export_to_csv(results: Dict[str, Any]) -> str:
    """
    Convert calculation results to CSV format.
    
    Args:
        results: Results from m/z calculation
        
    Returns:
        CSV content as string
    """
    output = io.StringIO()
    
    # Create a DataFrame with the basic results
    result_df = pd.DataFrame([{
        'Sequence': results['sequence'],
        'Charge State': results['charge'],
        'Mass (Da)': results['mass'],
        'm/z': results['mz'],
        'Formula': results['formula'],
        'Isotope Type': results['isotope_type']
    }])
    
    # Export to CSV
    result_df.to_csv(output, index=False)
    
    return output.getvalue()


def export_to_excel(results_df: pd.DataFrame) -> bytes:
    """
    Export DataFrame to Excel format.
    
    Args:
        results_df: DataFrame containing results
        
    Returns:
        Excel file content as bytes
    """
    output = io.BytesIO()
    
    # Create Excel writer
    with pd.ExcelWriter(output, engine='openpyxl') as writer:
        results_df.to_excel(writer, sheet_name='Results', index=False)
        
        # Add metadata sheet
        metadata = pd.DataFrame([{
            'Generated': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'Tool': 'Peptide m/z Calculator',
            'Version': '1.0.0'
        }])
        metadata.to_excel(writer, sheet_name='Metadata', index=False)
    
    # Get the content
    output.seek(0)
    
    return output.getvalue()


def export_to_json(results: Dict[str, Any]) -> str:
    """
    Export results to JSON format.
    
    Args:
        results: Calculation results
        
    Returns:
        JSON string
    """
    # Create a copy of the results to avoid modifying the original
    export_data = dict(results)
    
    # Format timestamp
    export_data['timestamp'] = datetime.now().isoformat()
    
    # Convert to JSON
    return json.dumps(export_data, indent=2)


def export_batch_results(batch_results: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Export batch processing results in multiple formats.
    
    Args:
        batch_results: List of calculation results
        
    Returns:
        Dictionary with different export formats
    """
    # Create a DataFrame from the results
    df = pd.DataFrame(batch_results)
    
    # Generate exports
    csv_data = df.to_csv(index=False)
    
    excel_buffer = io.BytesIO()
    with pd.ExcelWriter(excel_buffer, engine='openpyxl') as writer:
        df.to_excel(writer, sheet_name='Results', index=False)
        
        # Add summary sheet
        summary = pd.DataFrame({
            'Total Sequences': [len(batch_results)],
            'Average Mass': [df['Mass (Da)'].mean() if 'Mass (Da)' in df else None],
            'Average m/z': [df['m/z'].mean() if 'm/z' in df else None],
            'Generated': [datetime.now().strftime('%Y-%m-%d %H:%M:%S')]
        })
        summary.to_excel(writer, sheet_name='Summary', index=False)
    
    excel_buffer.seek(0)
    excel_data = excel_buffer.getvalue()
    
    # Generate JSON
    json_data = json.dumps({
        'results': batch_results,
        'summary': {
            'count': len(batch_results),
            'generated': datetime.now().isoformat()
        }
    }, indent=2)
    
    return {
        'csv': csv_data,
        'excel': excel_data,
        'json': json_data
    }


def export_fragment_ions(fragment_ions: Dict[str, List[Dict[str, Any]]]) -> str:
    """
    Export fragment ions to CSV format.
    
    Args:
        fragment_ions: Dictionary of fragment ions
        
    Returns:
        CSV content as string
    """
    output = io.StringIO()
    
    # Extract b-ions and y-ions
    b_ions = fragment_ions.get('b_ions', [])
    y_ions = fragment_ions.get('y_ions', [])
    
    # Create DataFrames
    b_df = pd.DataFrame(b_ions) if b_ions else pd.DataFrame()
    y_df = pd.DataFrame(y_ions) if y_ions else pd.DataFrame()
    
    # Add ion type column
    if not b_df.empty:
        b_df['ion_type'] = 'b'
    if not y_df.empty:
        y_df['ion_type'] = 'y'
    
    # Combine and sort by position
    combined_df = pd.concat([b_df, y_df]).sort_values(['ion_type', 'position'])
    
    # Export to CSV
    combined_df.to_csv(output, index=False)
    
    return output.getvalue()


def generate_report(results: Dict[str, Any], include_plots: bool = False) -> str:
    """
    Generate a full HTML report of the calculation results.
    
    Args:
        results: Calculation results
        include_plots: Whether to include base64-encoded plots
        
    Returns:
        HTML content as string
    """
    # Start building the HTML
    html = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Peptide m/z Calculator Report</title>
        <style>
            body {{
                font-family: Arial, sans-serif;
                line-height: 1.6;
                color: #333;
                max-width: 1000px;
                margin: 0 auto;
                padding: 20px;
            }}
            h1, h2, h3 {{
                color: #1E88E5;
            }}
            .result-box {{
                background-color: #f5f5f5;
                border-left: 4px solid #1E88E5;
                padding: 15px;
                margin: 20px 0;
            }}
            .formula {{
                font-family: monospace;
                background-color: #eef;
                padding: 5px;
            }}
            table {{
                border-collapse: collapse;
                width: 100%;
                margin: 20px 0;
            }}
            th, td {{
                border: 1px solid #ddd;
                padding: 8px;
                text-align: left;
            }}
            th {{
                background-color: #f2f2f2;
            }}
            .footer {{
                margin-top: 30px;
                padding-top: 10px;
                border-top: 1px solid #ddd;
                font-size: 0.8em;
                color: #777;
            }}
            img {{
                max-width: 100%;
                height: auto;
            }}
        </style>
    </head>
    <body>
        <h1>Peptide m/z Calculator Report</h1>
        <p>Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        
        <div class="result-box">
            <h2>Basic Information</h2>
            <table>
                <tr>
                    <th>Peptide Sequence</th>
                    <td>{results['sequence']}</td>
                </tr>
                <tr>
                    <th>Charge State</th>
                    <td>+{results['charge']}</td>
                </tr>
                <tr>
                    <th>Mass Type</th>
                    <td>{results['isotope_type'].capitalize()}</td>
                </tr>
                <tr>
                    <th>Mass (Da)</th>
                    <td>{results['mass']:.4f}</td>
                </tr>
                <tr>
                    <th>m/z Ratio</th>
                    <td>{results['mz']:.4f}</td>
                </tr>
                <tr>
                    <th>Chemical Formula</th>
                    <td class="formula">{results['formula']}</td>
                </tr>
            </table>
        </div>
    """
    
    # Add fragment ions section if available
    if 'fragment_ions' in results and results['fragment_ions']:
        b_ions = results['fragment_ions'].get('b_ions', [])
        y_ions = results['fragment_ions'].get('y_ions', [])
        
        html += f"""
        <div class="result-box">
            <h2>Fragment Ions</h2>
            
            <h3>b-ions</h3>
            <table>
                <tr>
                    <th>Position</th>
                    <th>m/z</th>
                    <th>Mass (Da)</th>
                    <th>Sequence</th>
                </tr>
        """
        
        for ion in b_ions:
            html += f"""
                <tr>
                    <td>b{ion['position']}</td>
                    <td>{ion['mz']:.4f}</td>
                    <td>{ion['mass']:.4f}</td>
                    <td>{ion.get('sequence', '')}</td>
                </tr>
            """
        
        html += f"""
            </table>
            
            <h3>y-ions</h3>
            <table>
                <tr>
                    <th>Position</th>
                    <th>m/z</th>
                    <th>Mass (Da)</th>
                    <th>Sequence</th>
                </tr>
        """
        
        for ion in y_ions:
            html += f"""
                <tr>
                    <td>y{ion['position']}</td>
                    <td>{ion['mz']:.4f}</td>
                    <td>{ion['mass']:.4f}</td>
                    <td>{ion.get('sequence', '')}</td>
                </tr>
            """
        
        html += """
            </table>
        </div>
        """
    
    # Add plots if requested and available
    if include_plots and 'plots' in results:
        html += """
        <div class="result-box">
            <h2>Visualizations</h2>
        """
        
        for plot_name, plot_data in results['plots'].items():
            html += f"""
            <h3>{plot_name}</h3>
            <img src="data:image/png;base64,{plot_data}" alt="{plot_name}">
            """
        
        html += """
        </div>
        """
    
    # Footer
    html += """
        <div class="footer">
            <p>Generated by Peptide m/z Calculator v1.0.0</p>
            <p>Powered by PyOpenMS</p>
        </div>
    </body>
    </html>
    """
    
    return html


def encode_image_base64(image_bytes: bytes) -> str:
    """
    Encode image bytes as base64 string.
    
    Args:
        image_bytes: Image data as bytes
        
    Returns:
        Base64 encoded string
    """
    return base64.b64encode(image_bytes).decode('utf-8')


def get_download_link(content: Union[str, bytes], filename: str, mime_type: str) -> str:
    """
    Generate an HTML download link for content.
    
    Args:
        content: File content as string or bytes
        filename: Name of the download file
        mime_type: MIME type of the file
        
    Returns:
        HTML anchor tag with download link
    """
    if isinstance(content, str):
        content_bytes = content.encode('utf-8')
    else:
        content_bytes = content
    
    b64_content = base64.b64encode(content_bytes).decode('utf-8')
    
    return f'<a href="data:{mime_type};base64,{b64_content}" download="{filename}">Download {filename}</a>'


def export_to_msp(
    results: Dict[str, Any],
    include_fragments: bool = True
) -> str:
    """
    Export results to MSP format for spectral libraries.
    
    Args:
        results: Calculation results
        include_fragments: Whether to include fragment ions
        
    Returns:
        MSP formatted string
    """
    msp_content = f"""NAME: {results['sequence']}
PEPTIDE: {results['sequence']}
CHARGE: {results['charge']}
PRECURSORMZ: {results['mz']:.6f}
PRECURSORTYPE: [M+{results['charge']}H]+
FORMULA: {results['formula']}
EXACTMASS: {results['mass']:.6f}
"""
    
    # Add fragment peaks if requested and available
    if include_fragments and 'fragment_ions' in results:
        # Start the peaks section
        peaks = []
        
        # Add b-ions
        for ion in results['fragment_ions'].get('b_ions', []):
            peaks.append((ion['mz'], 100))  # Using arbitrary intensity
        
        # Add y-ions
        for ion in results['fragment_ions'].get('y_ions', []):
            peaks.append((ion['mz'], 100))  # Using arbitrary intensity
        
        # Sort by m/z
        peaks.sort(key=lambda x: x[0])
        
        # Add peak count and peak data
        msp_content += f"NUMPEAKS: {len(peaks)}\n"
        
        for mz, intensity in peaks:
            msp_content += f"{mz:.6f} {intensity:.1f}\n"
    
    return msp_content


def export_to_mgf(results: Dict[str, Any]) -> str:
    """
    Export results to MGF format.
    
    Args:
        results: Calculation results
        
    Returns:
        MGF formatted string
    """
    mgf_content = f"""BEGIN IONS
TITLE={results['sequence']}
PEPMASS={results['mz']:.6f}
CHARGE={results['charge']}+
"""
    
    # Add fragment peaks if available
    if 'fragment_ions' in results:
        # Collect peaks
        peaks = []
        
        # Add b-ions
        for ion in results['fragment_ions'].get('b_ions', []):
            peaks.append((ion['mz'], 100, f"b{ion['position']}"))
        
        # Add y-ions
        for ion in results['fragment_ions'].get('y_ions', []):
            peaks.append((ion['mz'], 100, f"y{ion['position']}"))
        
        # Sort by m/z
        peaks.sort(key=lambda x: x[0])
        
        # Add peak data
        for mz, intensity, label in peaks:
            mgf_content += f"{mz:.6f} {intensity:.1f} # {label}\n"
    
    mgf_content += "END IONS\n"
    
    return mgf_content


def export_summary_statistics(results_df: pd.DataFrame) -> pd.DataFrame:
    """
    Generate summary statistics from batch results.
    
    Args:
        results_df: DataFrame with batch results
        
    Returns:
        DataFrame with summary statistics
    """
    # Ensure numeric columns
    numeric_cols = results_df.select_dtypes(include=['number']).columns
    
    # Calculate statistics
    stats = {}
    
    for col in numeric_cols:
        stats[col] = {
            'Mean': results_df[col].mean(),
            'Median': results_df[col].median(),
            'Min': results_df[col].min(),
            'Max': results_df[col].max(),
            'Std Dev': results_df[col].std()
        }
    
    # Convert to DataFrame
    stats_df = pd.DataFrame()
    
    for col, values in stats.items():
        col_df = pd.DataFrame(values.items(), columns=['Statistic', col])
        
        if stats_df.empty:
            stats_df = col_df
        else:
            stats_df = pd.merge(stats_df, col_df, on='Statistic')
    
    return stats_df


def export_to_spectrast(results: Dict[str, Any]) -> str:
    """
    Export results to SpectraST .sptxt format.
    
    Args:
        results: Calculation results
        
    Returns:
        SpectraST formatted string
    """
    sptxt_content = f"""Name: {results['sequence']}/{results['charge']}
MW: {results['mass']:.6f}
Comment: CHARGE={results['charge']}+ PEPTIDE={results['sequence']} FORMULA={results['formula']}
NumPeaks: """
    
    # Add fragment peaks if available
    if 'fragment_ions' in results:
        # Collect peaks
        peaks = []
        
        # Add b-ions
        for ion in results['fragment_ions'].get('b_ions', []):
            peaks.append((ion['mz'], 100, f"b{ion['position']}"))
        
        # Add y-ions
        for ion in results['fragment_ions'].get('y_ions', []):
            peaks.append((ion['mz'], 100, f"y{ion['position']}"))
        
        # Sort by m/z
        peaks.sort(key=lambda x: x[0])
        
        # Add peak count
        sptxt_content += f"{len(peaks)}\n"
        
        # Add peak data
        for mz, intensity, label in peaks:
            sptxt_content += f"{mz:.6f}\t{intensity:.1f}\t{label}\n"
    else:
        # No peaks
        sptxt_content += "0\n"
    
    return sptxt_content


def export_to_mztab(results: Dict[str, Any]) -> str:
    """
    Export results to mzTab format.
    
    Args:
        results: Calculation results
        
    Returns:
        mzTab formatted string
    """
    # Current date and time
    current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    mztab_content = f"""MTD	mzTab-version	1.0.0
MTD	mzTab-mode	Summary
MTD	mzTab-type	Quantification
MTD	description	Peptide m/z Calculator results
MTD	software[1]	[MS, MS:1000752, TOPP software, ]
MTD	software[1]-setting[1]	Peptide m/z Calculator v1.0.0
MTD	ms_run[1]-location	virtual_run
MTD	ms_run[1]-id_format	[MS, MS:1000774, multiple peak list nativeID format, ]
MTD	fixed_mod[1]	[MS, MS:1002453, No fixed modifications searched, ]
MTD	variable_mod[1]	[MS, MS:1002454, No variable modifications searched, ]

SMH	sequence	calculated_mass	charge	theoretical_mz	formula
SML	{results['sequence']}	{results['mass']:.6f}	{results['charge']}	{results['mz']:.6f}	{results['formula']}
"""
    
    return mztab_content


def export_to_traml(results: Dict[str, Any]) -> str:
    """
    Export results to TraML format for transition lists.
    
    Args:
        results: Calculation results
        
    Returns:
        TraML formatted string
    """
    # Current date and time
    current_time = datetime.now().strftime("%Y-%m-%dT%H:%M:%SZ")
    
    traml_content = f"""<?xml version="1.0" encoding="UTF-8"?>
<TraML version="1.0.0" xmlns="http://psi.hupo.org/ms/traml" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <cvList>
    <cv id="MS" fullName="Proteomics Standards Initiative Mass Spectrometry Ontology" URI="https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo" />
    <cv id="UO" fullName="Unit Ontology" URI="https://raw.githubusercontent.com/bio-ontology-research-group/unit-ontology/master/unit.obo" />
  </cvList>
  <SourceFileList>
    <SourceFile id="SF1" name="Peptide m/z Calculator" location="virtual">
      <cvParam cvRef="MS" accession="MS:1000796" name="SourceFile type" value="Peptide m/z Calculator output" />
    </SourceFile>
  </SourceFileList>
  <ProteinList>
    <Protein id="PROT1">
      <cvParam cvRef="MS" accession="MS:1000885" name="protein accession" value="VIRTUAL_PROTEIN" />
      <Sequence>{results['sequence']}</Sequence>
    </Protein>
  </ProteinList>
  <CompoundList>
    <Peptide id="PEPTIDE1" sequence="{results['sequence']}">
      <cvParam cvRef="MS" accession="MS:1000041" name="charge state" value="{results['charge']}" />
      <cvParam cvRef="MS" accession="MS:1000905" name="chemical formula" value="{results['formula']}" />
      <cvParam cvRef="MS" accession="MS:1000895" name="SMILES string" value="N/A" />
      <ProteinRef ref="PROT1" />
    </Peptide>
  </CompoundList>
  <TransitionList>
"""
    
    transition_id = 1
    
    # Add fragment ions if available
    if 'fragment_ions' in results:
        # Add b-ions
        for ion in results['fragment_ions'].get('b_ions', []):
            traml_content += f"""    <Transition id="TR{transition_id}" peptideRef="PEPTIDE1">
      <Precursor>
        <cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="{results['mz']:.6f}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z" />
      </Precursor>
      <Product>
        <cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="{ion['mz']:.6f}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z" />
        <cvParam cvRef="MS" accession="MS:1001220" name="frag: b ion" value="{ion['position']}" />
      </Product>
      <cvParam cvRef="MS" accession="MS:1000045" name="collision energy" value="35" />
    </Transition>
"""
            transition_id += 1
        
        # Add y-ions
        for ion in results['fragment_ions'].get('y_ions', []):
            traml_content += f"""    <Transition id="TR{transition_id}" peptideRef="PEPTIDE1">
      <Precursor>
        <cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="{results['mz']:.6f}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z" />
      </Precursor>
      <Product>
        <cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="{ion['mz']:.6f}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z" />
        <cvParam cvRef="MS" accession="MS:1001220" name="frag: y ion" value="{ion['position']}" />
      </Product>
      <cvParam cvRef="MS" accession="MS:1000045" name="collision energy" value="35" />
    </Transition>
"""
            transition_id += 1
    
    # Close the XML
    traml_content += """  </TransitionList>
</TraML>
"""
    
    return traml_content


def export_to_skyline(results: Dict[str, Any]) -> str:
    """
    Export results to a Skyline transition list.
    
    Args:
        results: Calculation results
        
    Returns:
        Skyline transition list formatted string
    """
    skyline_content = "Protein,Peptide Sequence,Precursor Charge,Fragment Ion,Fragment Charge,m/z\n"
    
    # Add precursor
    skyline_content += f"VIRTUAL_PROTEIN,{results['sequence']},{results['charge']},precursor,{results['charge']},{results['mz']:.6f}\n"
    
    # Add fragment ions if available
    if 'fragment_ions' in results:
        # Add b-ions
        for ion in results['fragment_ions'].get('b_ions', []):
            fragment_charge = 1  # Assuming singly charged fragments
            skyline_content += f"VIRTUAL_PROTEIN,{results['sequence']},{results['charge']},b{ion['position']},{fragment_charge},{ion['mz']:.6f}\n"
        
        # Add y-ions
        for ion in results['fragment_ions'].get('y_ions', []):
            fragment_charge = 1  # Assuming singly charged fragments
            skyline_content += f"VIRTUAL_PROTEIN,{results['sequence']},{results['charge']},y{ion['position']},{fragment_charge},{ion['mz']:.6f}\n"
    
    return skyline_content


if __name__ == "__main__":
    # Example usage
    example_results = {
        'sequence': 'PEPTIDE',
        'charge': 2,
        'mass': 799.3593,
        'mz': 400.1869,
        'formula': 'C37H56N8O15',
        'isotope_type': 'monoisotopic',
        'fragment_ions': {
            'b_ions': [
                {'position': 1, 'mass': 98.06, 'mz': 98.06, 'sequence': 'P'},
                {'position': 2, 'mass': 227.1, 'mz': 227.1, 'sequence': 'PE'}
            ],
            'y_ions': [
                {'position': 1, 'mass': 147.11, 'mz': 147.11, 'sequence': 'E'},
                {'position': 2, 'mass': 262.14, 'mz': 262.14, 'sequence': 'DE'}
            ]
        }
    }
    
    # Export to various formats
    csv_output = export_to_csv(example_results)
    json_output = export_to_json(example_results)
    msp_output = export_to_msp(example_results)
    mgf_output = export_to_mgf(example_results)
    skyline_output = export_to_skyline(example_results)
    
    print("CSV Output:")
    print(csv_output)
    
    print("\nJSON Output:")
    print(json_output)
    
    print("\nMSP Output:")
    print(msp_output)
    
    print("\nMGF Output:")
    print(mgf_output)
    
    print("\nSkyline Output:")
    print(skyline_output)