"""Schémas Pydantic pour validation"""
from pydantic import BaseModel, Field
from typing import List, Literal, Optional
from api.config import config

class AnalysisRequest(BaseModel):
    """Requête d'analyse - MODIFIÉE pour accepter proteinId"""
    proteinId: str = Field(..., description="UniProt accession (ex: P01189)")
    mode: Literal["strict", "permissive"] = Field(default="strict")
    signalPeptideLength: int = Field(default=config.DEFAULT_SIGNAL_PEPTIDE_LENGTH, ge=10, le=50)
    minCleavageSites: int = Field(default=config.DEFAULT_MIN_CLEAVAGE_SITES, ge=2, le=10)
    minCleavageSpacing: int = Field(default=config.DEFAULT_MIN_CLEAVAGE_SPACING, ge=1, le=20)
    maxPeptideLength: int = Field(default=100, ge=10, le=500, description="Longueur maximale des peptides (aa)")

class CleavageSite(BaseModel):
    """Site de clivage"""
    position: int
    motif: str
    index: int

class PeptideResult(BaseModel):
    """Peptide prédit"""
    sequence: str
    start: int
    end: int
    length: int
    inRange: bool
    cleavageMotif: str
    bioactivityScore: float = Field(ge=0, le=100)
    bioactivitySource: Literal["api", "heuristic", "none"]
    uniprotStatus: Literal["exact", "partial", "unknown"] = "unknown"
    uniprotName: Optional[str] = None
    uniprotNote: Optional[str] = None
    uniprotAccession: Optional[str] = None

class AnalysisResponse(BaseModel):
    """Résultat d'analyse"""
    sequenceLength: int
    cleavageSitesCount: int
    peptides: List[PeptideResult]
    peptidesInRange: int
    topPeptides: List[PeptideResult]
    cleavageSites: List[CleavageSite] = []
    mode: str
    proteinId: str = ""
    geneName: str = ""
    proteinName: str = ""

class HealthResponse(BaseModel):
    """État de santé"""
    status: Literal["ok", "degraded", "error"]
    timestamp: str
    version: str
    peptideranker_available: bool = True