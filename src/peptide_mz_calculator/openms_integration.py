"""
Create a new file called openms_integration.py to integrate with other OpenMS tools.
"""

import os
import subprocess
import tempfile
import pandas as pd
from typing import Dict, List, Any, Optional, Tuple
import pyopenms
import logging

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class OpenMSIntegration:
    """Integration with OpenMS tools and workflows."""
    
    @staticmethod
    def export_to_idxml(sequence: str, charge_state: int) -> str:
        """
        Export a peptide as an idXML file for use with OpenMS tools.
        
        Args:
            sequence: Peptide sequence
            charge_state: Charge state
            
        Returns:
            Path to the generated idXML file
        """
        # Create a temporary file
        fd, idxml_path = tempfile.mkstemp(suffix=".idXML")
        os.close(fd)
        
        try:
            # Create a protein identification
            protein_id = pyopenms.ProteinIdentification()
            protein_id.setIdentifier("Peptide_Calculator")
            
            # Set search parameters
            search_params = protein_id.getSearchParameters()
            search_params.db = "Peptide Calculator DB"
            protein_id.setSearchParameters(search_params)
            
            # Create a hit for this peptide
            peptide_id = pyopenms.PeptideIdentification()
            peptide_id.setIdentifier("Peptide_Calculator")
            
            # Create peptide hit
            hit = pyopenms.PeptideHit()
            hit.setSequence(pyopenms.AASequence.fromString(sequence))
            hit.setCharge(charge_state)
            hit.setScore(100.0)  # Perfect match
            
            # Add hit to identification
            peptide_id.insertHit(hit)
            
            # Create id data vector
            protein_ids = [protein_id]
            peptide_ids = [peptide_id]
            
            # Write to idXML
            pyopenms.IdXMLFile().store(idxml_path, protein_ids, peptide_ids)
            
            return idxml_path
        
        except Exception as e:
            logger.error(f"Error exporting to idXML: {str(e)}")
            if os.path.exists(idxml_path):
                os.remove(idxml_path)
            raise RuntimeError(f"Failed to export to idXML: {str(e)}")
    
    @staticmethod
    def run_peptide_indexer(idxml_path: str, fasta_path: str) -> str:
        """
        Run PeptideIndexer to map peptides to proteins.
        
        Args:
            idxml_path: Path to idXML file
            fasta_path: Path to FASTA database
            
        Returns:
            Path to output idXML file
        """
        # Create output path
        output_path = idxml_path.replace(".idXML", "_indexed.idXML")
        
        try:
            # Run PeptideIndexer
            cmd = [
                "PeptideIndexer",
                "-in", idxml_path,
                "-fasta", fasta_path,
                "-out", output_path,
                "-allow_unmatched", "true",
                "-missing_decoy_action", "warn"
            ]
            
            process = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            if process.returncode != 0:
                logger.error(f"PeptideIndexer failed: {process.stderr}")
                raise RuntimeError(f"PeptideIndexer failed: {process.stderr}")
            
            return output_path
        
        except subprocess.CalledProcessError as e:
            logger.error(f"PeptideIndexer process error: {e.stderr}")
            raise RuntimeError(f"PeptideIndexer process error: {e.stderr}")
        
        except Exception as e:
            logger.error(f"Error running PeptideIndexer: {str(e)}")
            raise RuntimeError(f"Error running PeptideIndexer: {str(e)}")
    
    @staticmethod
    def generate_theoretical_spectrum(
        sequence: str,
        charge_state: int,
        activation_method: str = "HCD"
    ) -> Tuple[List[float], List[float]]:
        """
        Generate a theoretical spectrum using OpenMS TheoreticalSpectrumGenerator.
        
        Args:
            sequence: Peptide sequence
            charge_state: Charge state
            activation_method: Fragmentation method
            
        Returns:
            Tuple of m/z values and intensities
        """
        try:
            peptide = pyopenms.AASequence.fromString(sequence)
            
            # Create generator
            tsg = pyopenms.TheoreticalSpectrumGenerator()
            
            # Set parameters based on activation method
            params = tsg.getParameters()
            if activation_method == "HCD":
                params.setValue("add_b_ions", "true")
                params.setValue("add_y_ions", "true")
                params.setValue("add_metainfo", "true")
            elif activation_method == "ETD":
                params.setValue("add_c_ions", "true")
                params.setValue("add_z_ions", "true")
                params.setValue("add_metainfo", "true")
            elif activation_method == "CID":
                params.setValue("add_b_ions", "true")
                params.setValue("add_y_ions", "true")
                params.setValue("add_a_ions", "true")
                params.setValue("add_metainfo", "true")
            
            tsg.setParameters(params)
            
            # Generate spectrum
            spectrum = pyopenms.MSSpectrum()
            tsg.getSpectrum(spectrum, peptide, min_charge=1, max_charge=charge_state)
            
            # Extract m/z and intensity values
            mz_values = []
            intensities = []
            annotations = []
            
            for peak in spectrum:
                mz_values.append(peak.getMZ())
                intensities.append(peak.getIntensity())
                
                # Get annotation if available
                if peak.metaValueExists("IonName"):
                    annotations.append(peak.getMetaValue("IonName"))
                else:
                    annotations.append("")
            
            # Normalize intensities
            if intensities:
                max_intensity = max(intensities)
                intensities = [i / max_intensity * 100 for i in intensities]
            
            return mz_values, intensities, annotations
        
        except Exception as e:
            logger.error(f"Error generating theoretical spectrum: {str(e)}")
            raise RuntimeError(f"Error generating theoretical spectrum: {str(e)}")
    
    @staticmethod
    def export_to_traml(
        peptide_data: List[Dict[str, Any]],
        output_path: Optional[str] = None
    ) -> str:
        """
        Export peptide data to TraML format for targeted proteomics.
        
        Args:
            peptide_data: List of peptide dictionaries with sequence, precursor m/z, fragments
            output_path: Optional path for output file
            
        Returns:
            Path to the generated TraML file
        """
        if output_path is None:
            fd, output_path = tempfile.mkstemp(suffix=".traML")
            os.close(fd)
        
        try:
            # Create TraML document
            traml = pyopenms.TraMLFile()
            targeted_exp = pyopenms.TargetedExperiment()
            
            # Set CV terms
            cv_terms = targeted_exp.getCVTerms()
            
            # Create protein and peptide lists
            proteins = []
            compounds = []
            transitions = []
            
            for i, peptide in enumerate(peptide_data):
                sequence = peptide["sequence"]
                charge = peptide.get("charge", 2)
                mz = peptide.get("mz", 0.0)
                rt = peptide.get("rt", -1.0)
                fragments = peptide.get("fragments", [])
                
                # Create protein
                protein = pyopenms.TargetedExperiment.Protein()
                protein.id = f"PROT_{i+1}"
                protein.sequence = sequence
                proteins.append(protein)
                
                # Create peptide/compound
                compound = pyopenms.TargetedExperiment.Peptide()
                compound.id = f"PEP_{i+1}"
                compound.sequence = sequence
                compound.charge = charge
                if rt > 0:
                    rt_term = pyopenms.CVTerm()
                    rt_term.setCVIdentifierRef("MS")
                    rt_term.setAccession("MS:1000895")
                    rt_term.setName("local retention time")
                    rt_term.setValue(pyopenms.DataValue(rt))
                    compound.addCVTerm(rt_term)
                
                # Add protein reference
                ref = pyopenms.TargetedExperiment.PeptideCompoundRef()
                ref.compound_ref = protein.id
                compound.protein_refs.append(ref)
                
                compounds.append(compound)
                
                # Create transitions
                for j, fragment in enumerate(fragments):
                    transition = pyopenms.TargetedExperiment.Transition()
                    transition.id = f"TRANSITION_{i+1}_{j+1}"
                    transition.peptide_ref = compound.id
                    
                    # Set precursor
                    precursor = pyopenms.TargetedExperiment.Precursor()
                    precursor.mz = mz
                    precursor.charge = charge
                    transition.setPrecursor(precursor)
                    
                    # Set product
                    product = pyopenms.TargetedExperiment.Product()
                    product.mz = fragment.get("mz", 0.0)
                    product.charge = fragment.get("charge", 1)
                    
                    # Set interpretation (ion type)
                    interpretation = pyopenms.TargetedExperiment.Interpretation()
                    interpretation.rank = 1
                    
                    ion_type = fragment.get("type", "")
                    if ion_type.startswith("b"):
                        interpretation.ordinal = int(ion_type[1:])
                        interpretation.iontype = pyopenms.IonType.BIon
                    elif ion_type.startswith("y"):
                        interpretation.ordinal = int(ion_type[1:])
                        interpretation.iontype = pyopenms.IonType.YIon
                    elif ion_type.startswith("a"):
                        interpretation.ordinal = int(ion_type[1:])
                        interpretation.iontype = pyopenms.IonType.AIon
                    
                    product.interpretations.append(interpretation)
                    transition.setProduct(product)
                    
                    transitions.append(transition)
            
            # Set proteins, compounds, and transitions
            targeted_exp.setProteins(proteins)
            targeted_exp.setPeptides(compounds)
            targeted_exp.setTransitions(transitions)
            
            # Store TraML file
            traml.store(output_path, targeted_exp)
            
            return output_path
        
        except Exception as e:
            logger.error(f"Error exporting to TraML: {str(e)}")
            if os.path.exists(output_path):
                os.remove(output_path)
            raise RuntimeError(f"Failed to export to TraML: {str(e)}")
    
    @staticmethod
    def predict_retention_time(
        sequences: List[str],
        model_path: Optional[str] = None
    ) -> Dict[str, float]:
        """
        Predict retention times for peptides using OpenMS RTModel.
        
        Args:
            sequences: List of peptide sequences
            model_path: Path to trained model file (SSRCalibration)
            
        Returns:
            Dictionary mapping sequences to predicted retention times
        """
        try:
            # Create RTModel
            rt_model = pyopenms.RTModel()
            
            # Load model if provided, otherwise use default model
            if model_path and os.path.exists(model_path):
                param = pyopenms.Param()
                param.setValue("filename", model_path)
                rt_model.setParameters(param)
            
            # Initialize model
            rt_model.init(False)
            
            # Predict retention times
            results = {}
            for sequence in sequences:
                peptide = pyopenms.AASequence.fromString(sequence)
                predicted_rt = rt_model.computeRT(peptide)
                results[sequence] = predicted_rt
            
            return results
        
        except Exception as e:
            logger.error(f"Error predicting retention times: {str(e)}")
            raise RuntimeError(f"Error predicting retention times: {str(e)}")
    
    @staticmethod
    def match_against_spectral_library(
        mz_values: List[float],
        intensities: List[float],
        library_path: str,
        precursor_mz: float,
        precursor_charge: int,
        tolerance: float = 0.01,
        tolerance_unit: str = "Da"
    ) -> List[Dict[str, Any]]:
        """
        Match a spectrum against a spectral library.
        
        Args:
            mz_values: Experimental m/z values
            intensities: Experimental intensities
            library_path: Path to spectral library (MSP or MGF)
            precursor_mz: Precursor m/z
            precursor_charge: Precursor charge
            tolerance: Mass tolerance for matching
            tolerance_unit: Unit for tolerance ("Da" or "ppm")
            
        Returns:
            List of matched spectra with scores
        """
        try:
            # Create MSPFile parser
            spec_lib = pyopenms.MSPFile()
            
            # Create the experimental spectrum
            exp_spectrum = pyopenms.MSSpectrum()
            for mz, intensity in zip(mz_values, intensities):
                peak = pyopenms.Peak1D()
                peak.setMZ(mz)
                peak.setIntensity(intensity)
                exp_spectrum.push_back(peak)
            
            # Set precursor information
            precursor = pyopenms.Precursor()
            precursor.setMZ(precursor_mz)
            precursor.setCharge(precursor_charge)
            exp_spectrum.setPrecursors([precursor])
            
            # Sort by m/z
            exp_spectrum.sortByPosition()
            
            # Load library spectra
            library_spectra = []
            spec_lib.load(library_path, library_spectra)
            
            # Match against library
            results = []
            for lib_spectrum in library_spectra:
                # Check if precursor m/z is within tolerance
                lib_precursors = lib_spectrum.getPrecursors()
                if not lib_precursors:
                    continue
                
                lib_mz = lib_precursors[0].getMZ()
                
                # Calculate mass difference
                if tolerance_unit.lower() == "ppm":
                    mass_diff = abs(lib_mz - precursor_mz) / precursor_mz * 1e6
                else:  # Da
                    mass_diff = abs(lib_mz - precursor_mz)
                
                # Skip if precursor doesn't match
                if mass_diff > tolerance:
                    continue
                
                # Calculate similarity score
                score = pyopenms.SpectrumAlignment.getSpectrumAlignment(
                    exp_spectrum, lib_spectrum)
                
                # Get spectrum title/name if available
                spectrum_name = "Unknown"
                if lib_spectrum.metaValueExists("Name"):
                    spectrum_name = lib_spectrum.getMetaValue("Name")
                
                results.append({
                    "name": spectrum_name,
                    "score": score,
                    "precursor_mz": lib_mz,
                    "mass_diff": mass_diff
                })
            
            # Sort by score (descending)
            results.sort(key=lambda x: x["score"], reverse=True)
            
            return results
        
        except Exception as e:
            logger.error(f"Error matching against spectral library: {str(e)}")
            raise RuntimeError(f"Error matching against spectral library: {str(e)}")