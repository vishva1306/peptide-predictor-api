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
    PTMResult,
    HealthResponse
)
from api.services import (
    SequenceValidator,
    CleavageDetector,
    PeptideExtractor,
    BioactivityPredictor,
    UniProtChecker,
    protein_db,
    ptm_detector
)
from api.routes import proteins

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

app.include_router(proteins.router, prefix="/api", tags=["proteins"])

# ==================== ROUTES ====================

@app.get("/")
async def root():
    return {
        "message": f"{config.API_TITLE} v{config.API_VERSION}",
        "version": config.API_VERSION,
        "paper": "Coassolo et al. Nature 2025",
        "endpoints": {
            "search_proteins": "GET /api/proteins/search?q=POMC",
            "get_protein": "GET /api/proteins/P01189",
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
            "fallback": "Lab ML Bioactivity Model",
            "note": "Le papier original utilise RT-qPCR expérimental"
        },
        "ptms": {
            "enabled": True,
            "count": 6,
            "types": [
                "C-terminal amidation",
                "N-terminal pyroglutamate",
                "Disulfide bonds",
                "Ghrelin acylation",
                "Tyrosine O-sulfation",
                "N-glycosylation"
            ]
        },
        "uniprot": {
            "enabled": True,
            "description": "Vérification automatique des peptides connus",
            "database": "UniProtKB/Swiss-Prot (reviewed entries only)",
            "statuses": ["exact", "partial", "unknown"]
        },
        "protein_search": {
            "enabled": True,
            "description": "Recherche par gene name ou UniProt ID",
            "scope": "Human secreted proteins only",
            "cache": "24 hours"
        },
        "default_params": {
            "signalPeptideLength": config.DEFAULT_SIGNAL_PEPTIDE_LENGTH,
            "minCleavageSites": config.DEFAULT_MIN_CLEAVAGE_SITES,
            "minCleavageSpacing": config.DEFAULT_MIN_CLEAVAGE_SPACING,
            "maxPeptideLength": 100
        },
        "optimal_range": f"{config.OPTIMAL_PEPTIDE_MIN_LENGTH}-{config.OPTIMAL_PEPTIDE_MAX_LENGTH} aa"
    }


@app.post("/analyze", response_model=AnalysisResponse)
async def analyze(request: AnalysisRequest):
    """
    Analyse une protéine pour prédire les peptides bioactifs
    
    Pipeline:
    1. Récupération de la protéine depuis UniProt (avec cache)
    2. Validation de la séquence
    3. Détection sites de clivage PCSK1/3
    4. Extraction peptides (avec filtre maxPeptideLength)
    5. Prédiction bioactivité (parallèle)
    6. Vérification UniProt (parallèle)
    7. Détection PTMs
    8. Tri et sélection top candidats
    """
    try:
        # 1. Récupérer la protéine depuis UniProt
        async with aiohttp.ClientSession() as session:
            protein = await protein_db.get_protein(request.proteinId, session)
        
        if not protein:
            raise HTTPException(
                status_code=404,
                detail=f"Protein {request.proteinId} not found or not secreted"
            )
        
        # 2. Extraire la séquence
        clean_seq = protein["sequence"]
        gene_name = protein["geneName"]
        protein_name = protein["proteinName"]
        
        protein_id_header = f"SP|{protein['accession']}|{gene_name}_HUMAN {protein_name}"
        
        print(f"\n🧬 Analyzing protein: {gene_name} ({protein['accession']})")
        print(f"📊 Length: {protein['length']} aa")
        print(f"📋 Parameters: signal={request.signalPeptideLength}, min_sites={request.minCleavageSites}, spacing={request.minCleavageSpacing}, max_len={request.maxPeptideLength}")
        
        # 3. Valider la séquence
        SequenceValidator.validate_characters(clean_seq)
        min_length = request.signalPeptideLength + 10
        SequenceValidator.validate_length(clean_seq, min_length)
        
        # 4. Détecter sites de clivage
        cleavage_sites = CleavageDetector.find_sites(
            sequence=clean_seq,
            mode=request.mode,
            signal_length=request.signalPeptideLength,
            min_spacing=request.minCleavageSpacing
        )
        
        # 5. Extraire peptides
        peptides = PeptideExtractor.extract(
            sequence=clean_seq,
            cleavage_sites=cleavage_sites,
            signal_length=request.signalPeptideLength,
            min_spacing=request.minCleavageSpacing,
            min_sites=request.minCleavageSites,
            mode=request.mode
        )
        
        # 5.5. Filtrer par maxPeptideLength
        peptides_filtered = [
            p for p in peptides
            if p['length'] <= request.maxPeptideLength
        ]
        
        print(f"📊 Peptides before max filter: {len(peptides)}")
        print(f"📊 Peptides after max filter: {len(peptides_filtered)}")
        
        peptides = peptides_filtered
        
        # 6. Session aiohttp pour bioactivité + UniProt
        async with aiohttp.ClientSession() as session:
            # Calculer bioactivité (parallèle)
            bioactivity_results = await BioactivityPredictor.predict_batch(
                [p['sequence'] for p in peptides],
                session
            )
            
            # Vérifier UniProt (parallèle)
            uniprot_results = await UniProtChecker.check_batch(
                [p['sequence'] for p in peptides],
                session,
                protein_id=protein_id_header
            )
        
        # 7. Assigner scores bioactivité
        for peptide, (score, source) in zip(peptides, bioactivity_results):
            peptide['bioactivityScore'] = score
            peptide['bioactivitySource'] = source
        
        # 8. Assigner données UniProt
        for peptide, uniprot_data in zip(peptides, uniprot_results):
            peptide['uniprotStatus'] = uniprot_data['uniprotStatus']
            peptide['uniprotName'] = uniprot_data['uniprotName']
            peptide['uniprotNote'] = uniprot_data['uniprotNote']
            peptide['uniprotAccession'] = uniprot_data['uniprotAccession']
        
        # 9. Détecter PTMs pour chaque peptide
        print(f"\n🔬 Detecting PTMs for {len(peptides)} peptides...")
        for idx, peptide in enumerate(peptides, 1):
            try:
                # Vérifier que start et end sont des int
                if not isinstance(peptide['start'], int) or not isinstance(peptide['end'], int):
                    print(f"  ⚠️ Peptide #{idx}: Invalid start/end types")
                    peptide['ptms'] = []
                    peptide['modifiedSequence'] = None
                    continue
                
                detected_ptms = ptm_detector.detect_all_ptms(
                    peptide_sequence=peptide['sequence'],
                    full_protein_sequence=clean_seq,
                    peptide_start=peptide['start'],
                    peptide_end=peptide['end']
                )
                
                peptide['ptms'] = detected_ptms
                
                # Générer séquence modifiée si PTMs détectées
                if detected_ptms:
                    peptide['modifiedSequence'] = ptm_detector.generate_modified_sequence(
                        peptide['sequence'],
                        detected_ptms
                    )
                    print(f"  ✅ Peptide #{idx}: {len(detected_ptms)} PTMs detected")
                else:
                    peptide['modifiedSequence'] = None
                    
            except Exception as e:
                print(f"  ❌ Peptide #{idx}: PTM detection error: {e}")
                import traceback
                traceback.print_exc()
                peptide['ptms'] = []
                peptide['modifiedSequence'] = None
        
        # 10. Trier par bioactivité
        peptides.sort(key=lambda x: x['bioactivityScore'], reverse=True)
        
        # 11. Top 5 peptides
        top_peptides = peptides[:5]
        
        # 12. Stats
        peptides_in_range = sum(1 for p in peptides if p['inRange'])
        
        # 13. Construire réponse
        return AnalysisResponse(
            sequenceLength=len(clean_seq),
            cleavageSitesCount=len(cleavage_sites),
            peptides=[PeptideResult(**p) for p in peptides],
            peptidesInRange=peptides_in_range,
            topPeptides=[PeptideResult(**p) for p in top_peptides],
            cleavageSites=cleavage_sites,
            mode=request.mode,
            proteinId=protein['accession'],
            geneName=gene_name,
            proteinName=protein_name
        )
    
    except HTTPException:
        raise
    except Exception as e:
        print(f"❌ Analysis error: {e}")
        import traceback
        traceback.print_exc()
        raise HTTPException(
            status_code=500,
            detail=f"Erreur interne: {str(e)}"
        )


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)