"""
Create a new file called api.py to provide API endpoints for your application.
"""

import json
from typing import Dict, List, Any, Optional, Union
from fastapi import FastAPI, HTTPException, Query, Depends, Request
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field
import uvicorn
import logging
from contextlib import asynccontextmanager

# Import your existing calculator functionality
from calculator import calculate_mz, generate_isotope_pattern
from parser import parse_fasta, parse_file


# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler()]
)
logger = logging.getLogger(__name__)


# Define data models
class Modification(BaseModel):
    name: str = Field(..., description="Modification name (e.g., 'Phospho')")
    position: str = Field(..., description="Position in the peptide (e.g., '3' or 'N-term')")


class PeptideRequest(BaseModel):
    sequence: str = Field(..., description="Peptide sequence in single-letter amino acid code")
    charge_state: int = Field(2, description="Charge state of the peptide")
    modifications: Optional[List[Modification]] = Field(None, description="List of modifications to apply")
    isotope_type: str = Field("monoisotopic", description="Type of isotope mass to use")


class BatchPeptideRequest(BaseModel):
    sequences: List[str] = Field(..., description="List of peptide sequences")
    charge_state: int = Field(2, description="Charge state to apply to all peptides")
    isotope_type: str = Field("monoisotopic", description="Type of isotope mass to use")


# Setup the FastAPI app
@asynccontextmanager
async def lifespan(app: FastAPI):
    logger.info("Starting peptide m/z calculator API...")
    yield
    logger.info("Shutting down peptide m/z calculator API...")


app = FastAPI(
    title="Peptide m/z Calculator API",
    description="API for calculating m/z ratios and other properties of peptides",
    version="1.0.0",
    lifespan=lifespan
)

# Enable CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # In production, you would specify the allowed origins
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.get("/")
async def root():
    return {"message": "Welcome to the Peptide m/z Calculator API", 
            "docs_url": "/docs"}


@app.post("/calculate")
async def calculate_peptide(request: PeptideRequest):
    """
    Calculate m/z and other properties for a peptide sequence.
    """
    try:
        result = calculate_mz(
            sequence=request.sequence,
            charge_state=request.charge_state,
            modifications=request.modifications,
            isotope_type=request.isotope_type
        )
        
        # Get isotope pattern
        if result:
            masses, intensities = generate_isotope_pattern(
                request.sequence,
                request.charge_state
            )
            
            result["isotope_pattern"] = {
                "masses": masses,
                "intensities": intensities
            }
        
        return result
    
    except Exception as e:
        logger.error(f"Error calculating m/z: {str(e)}")
        raise HTTPException(status_code=400, detail=str(e))


@app.post("/batch-calculate")
async def batch_calculate(request: BatchPeptideRequest):
    """
    Calculate m/z and other properties for multiple peptide sequences.
    """
    results = []
    errors = []
    
    for sequence in request.sequences:
        try:
            result = calculate_mz(
                sequence=sequence,
                charge_state=request.charge_state,
                modifications=None,  # No modifications for batch processing
                isotope_type=request.isotope_type
            )
            results.append(result)
        
        except Exception as e:
            errors.append({
                "sequence": sequence,
                "error": str(e)
            })
    
    return {
        "results": results,
        "errors": errors,
        "total": len(request.sequences),
        "successful": len(results),
        "failed": len(errors)
    }


@app.post("/parse-fasta")
async def parse_fasta_endpoint(fasta_content: str = Query(..., description="FASTA format content")):
    """
    Parse FASTA content and return sequences.
    """
    try:
        sequences = parse_fasta(fasta_content)
        return {"sequences": sequences}
    
    except Exception as e:
        logger.error(f"Error parsing FASTA: {str(e)}")
        raise HTTPException(status_code=400, detail=str(e))


# Add additional endpoints as needed...


def run_api_server(host: str = "0.0.0.0", port: int = 8000):
    """
    Run the FastAPI server.
    """
    uvicorn.run(app, host=host, port=port)


if __name__ == "__main__":
    run_api_server()