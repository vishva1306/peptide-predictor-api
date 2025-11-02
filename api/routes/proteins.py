"""Routes API pour la recherche de protéines"""
from fastapi import APIRouter, HTTPException, Query
from typing import Literal
import aiohttp

from api.services.protein_db import protein_db

router = APIRouter()


@router.get("/proteins/search")
async def search_proteins(
    q: str = Query(..., min_length=2, description="Gene name or UniProt ID"),
    type: Literal["gene_name", "accession"] = Query("gene_name", description="Search type"),
    limit: int = Query(10, ge=1, le=20, description="Max results")
):
    """
    Recherche de protéines sécrétées humaines
    
    Args:
        q: Gene name (ex: "POMC", "INS") ou UniProt ID (ex: "P01189")
        type: "gene_name" ou "accession"
        limit: Nombre maximum de résultats (1-20)
    
    Returns:
        Liste de protéines matchant la requête
        
    Examples:
        - /proteins/search?q=POMC&type=gene_name
        - /proteins/search?q=P01189&type=accession
    """
    
    async with aiohttp.ClientSession() as session:
        proteins = await protein_db.search_proteins(q, type, session, limit)
    
    if not proteins:
        return []
    
    # Retourner format complet pour sélection
    return [
        {
            "accession": p["accession"],
            "geneName": p["geneName"],
            "proteinName": p["proteinName"],
            "length": p["length"],
            "signalPeptideEnd": p["signalPeptideEnd"],
            "fastaHeader": p["fastaHeader"]
        }
        for p in proteins
    ]


@router.get("/proteins/{accession}")
async def get_protein(accession: str):
    """Récupère les détails complets d'une protéine"""
    
    async with aiohttp.ClientSession() as session:
        protein = await protein_db.get_protein(accession, session)
    
    if not protein:
        raise HTTPException(
            status_code=404,
            detail=f"Protein {accession} not found or not secreted"
        )
    
    return protein