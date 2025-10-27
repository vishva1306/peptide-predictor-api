"""API FastAPI refactorisée"""
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from datetime import datetime
import aiohttp

from api.config import config
from api.models.schemas import (
    AnalysisRequest,
    AnalysisResponse,
    PeptideResult,
    HealthResponse
)
from api.services import (
    SequenceValidator,
    CleavageDetector,
    PeptideExtractor,
    BioactivityPredictor
)

# ==================== APP ====================

app = FastAPI(
    title=config.API_TITLE,
    version=config.API_VERSION,
    description=config.API_DESCRIPTION
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=config.CORS_ORIGINS,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# ==================== ROUTES ====================

@app.get("/")
async def root():
    return {
        "message": f"{config.API_TITLE} v{config.API_VERSION}",
        "version": config.API_VERSION,
        "paper": "Coassolo et al. Nature 2025",
        "endpoints": {
            "analyze": "POST /analyze",
            "health": "GET /health",
            "stats": "GET /stats",
            "docs": "/docs"
        }
    }


@app.get("/health", response_model=HealthResponse)
async def health():
    """Vérifie l'état de santé de l'API"""
    return HealthResponse(
        status="ok",
        timestamp=datetime.now().isoformat(),
        version=config.API_VERSION,
        peptideranker_available=True
    )


@app.get("/stats")
async def stats():
    """Statistiques de l'API"""
    return {
        "api": config.API_TITLE,
        "version": config.API_VERSION,
        "paper": "Coassolo et al. Nature 2025 - doi:10.1038/s41586-025-08683-y",
        "modes": {
            "strict": "Regex complète du papier (recommandé)",
            "permissive": "Détection simplifiée (plus sensible)"
        },
        "bioactivity": {
            "primary": "PeptideRanker API (parallel)",
            "fallback": "Heuristic scoring",
            "note": "Le papier original utilise RT-qPCR expérimental"
        },
        "default_params": {
            "signalPeptideLength": config.DEFAULT_SIGNAL_PEPTIDE_LENGTH,
            "minCleavageSites": config.DEFAULT_MIN_CLEAVAGE_SITES,
            "minCleavageSpacing": config.DEFAULT_MIN_CLEAVAGE_SPACING
        },
        "optimal_range": f"{config.OPTIMAL_PEPTIDE_MIN_LENGTH}-{config.OPTIMAL_PEPTIDE_MAX_LENGTH} aa"
    }


@app.post("/analyze", response_model=AnalysisResponse)
async def analyze(request: AnalysisRequest):
    """
    Analyse une séquence protéique pour prédire les peptides bioactifs.
    
    Pipeline:
    1. Validation et nettoyage
    2. Détection sites de clivage PCSK1/3
    3. Extraction peptides
    4. Prédiction bioactivité (parallèle)
    5. Tri et sélection top candidats
    """
    try:
        # 1. Nettoyer et valider
        clean_seq = SequenceValidator.clean_sequence(request.sequence)
        SequenceValidator.validate_characters(clean_seq)
        
        min_length = request.signalPeptideLength + 10
        SequenceValidator.validate_length(clean_seq, min_length)
        
        # 2. Détecter sites de clivage
        cleavage_sites = CleavageDetector.find_sites(
            sequence=clean_seq,
            mode=request.mode,
            signal_length=request.signalPeptideLength,
            min_spacing=request.minCleavageSpacing
        )
        
        # 3. Extraire peptides
        peptides = PeptideExtractor.extract(
            sequence=clean_seq,
            cleavage_sites=cleavage_sites,
            signal_length=request.signalPeptideLength,
            min_spacing=request.minCleavageSpacing,
            min_sites=request.minCleavageSites,
            mode=request.mode  # ⭐ AJOUT DU PARAMÈTRE MODE
        )
        
        # 4. Calculer bioactivité (parallèle)
        async with aiohttp.ClientSession() as session:
            bioactivity_results = await BioactivityPredictor.predict_batch(
                [p['sequence'] for p in peptides],
                session
            )
        
        # Assigner scores
        for peptide, (score, source) in zip(peptides, bioactivity_results):
            peptide['bioactivityScore'] = score
            peptide['bioactivitySource'] = source
        
        # 5. Trier par bioactivité
        peptides.sort(key=lambda x: x['bioactivityScore'], reverse=True)
        
        # 6. Top 5 peptides
        top_peptides = peptides[:5]
        
        # 7. Stats
        peptides_in_range = sum(1 for p in peptides if p['inRange'])
        
        # 8. Construire réponse
        return AnalysisResponse(
            sequenceLength=len(clean_seq),
            cleavageSitesCount=len(cleavage_sites),
            peptides=[PeptideResult(**p) for p in peptides],
            peptidesInRange=peptides_in_range,
            topPeptides=[PeptideResult(**p) for p in top_peptides],
            cleavageSites=cleavage_sites,
            mode=request.mode
        )
    
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Erreur interne: {str(e)}"
        )


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)