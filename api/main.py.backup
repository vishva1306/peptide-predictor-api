from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import List, Optional
import re
import json
from datetime import datetime
import asyncio
import aiohttp
from concurrent.futures import ThreadPoolExecutor

app = FastAPI(
    title="Peptide Predictor API",
    version="2.0.0",
    description="API pour prédire les peptides bioactifs avec bioactivité ML"
)

# CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# ==================== MODELS ====================

class AnalysisRequest(BaseModel):
    sequence: str
    mode: str = "strict"
    signalPeptideLength: int = 20
    minCleavageSites: int = 4
    minCleavageSpacing: int = 5

class PeptideResult(BaseModel):
    sequence: str
    start: int
    end: int
    length: int
    inRange: bool
    cleavageMotif: str
    bioactivityScore: float
    bioactivitySource: str  # "api" ou "heuristic"

class AnalysisResponse(BaseModel):
    sequenceLength: int
    cleavageSitesCount: int
    peptides: List[PeptideResult]
    peptidesInRange: int
    topPeptides: List[PeptideResult]

# ==================== BIOACTIVITÉ ====================

async def get_peptideranker_score(peptide: str, session: aiohttp.ClientSession) -> Optional[dict]:
    """
    Appelle l'API PeptideRanker en ligne pour prédire la bioactivité
    Retourne un dict avec le score et la source
    """
    if len(peptide) < 2:
        return None
    
    try:
        url = "http://peptideranker.ilincs.org/api/predict"
        
        payload = {
            "sequence": peptide
        }
        
        async with session.post(url, json=payload, timeout=aiohttp.ClientTimeout(total=10)) as response:
            if response.status == 200:
                data = await response.json()
                # PeptideRanker retourne un score 0-1
                score = data.get('score', 0)
                return {
                    "score": float(score * 100),
                    "source": "api"
                }
    except asyncio.TimeoutError:
        print(f"Timeout calling PeptideRanker API for peptide: {peptide}")
    except Exception as e:
        print(f"Error calling PeptideRanker API: {e}")
    
    return None

def calculate_bioactivity_heuristic(peptide: str) -> float:
    """
    Fallback heuristique si API échoue
    Basé sur propriétés physicochimiques
    """
    score = 0.0
    
    if len(peptide) == 0:
        return 0.0
    
    # Hydrophobicité (30%)
    hydrophobic = ['A', 'L', 'I', 'V', 'M', 'F', 'W', 'P']
    hydro_count = sum(1 for aa in peptide if aa in hydrophobic)
    hydro_ratio = hydro_count / len(peptide)
    score += hydro_ratio * 30
    
    # Charge (20%)
    positive = sum(1 for aa in peptide if aa in ['K', 'R', 'H'])
    negative = sum(1 for aa in peptide if aa in ['D', 'E'])
    
    if positive > 0:
        score += 10
    if negative > 0:
        score += 10
    
    # Longueur gamme optimal (35%)
    if 5 <= len(peptide) <= 25:
        score += 35
    elif len(peptide) < 5:
        score -= 10
    elif len(peptide) > 100:
        score -= 15
    
    # Stabilité (15%)
    if 'C' in peptide:
        score += 8
    
    if peptide.count('P') <= 2:
        score += 7
    else:
        score -= 5
    
    unique_aa = len(set(peptide))
    if unique_aa >= 6:
        score += 5
    
    return min(max(score, 0), 100)

async def calculate_bioactivity_async(peptide: str, session: aiohttp.ClientSession) -> tuple:
    """
    Calcule bioactivité de manière asynchrone
    Retourne (score, source)
    """
    # Essayer l'API PeptideRanker d'abord
    api_result = await get_peptideranker_score(peptide, session)
    
    if api_result:
        return api_result["score"], api_result["source"]
    
    # Fallback sur heuristique
    heuristic_score = calculate_bioactivity_heuristic(peptide)
    return heuristic_score, "heuristic"

# ==================== DÉTECTION SITES CLIVAGE ====================

def find_cleavage_sites(seq: str, mode: str, params: dict) -> List[dict]:
    """
    Détecte les sites de clivage PCSK1/3
    Mode STRICT: regex complète du papier Nature
    Mode PERMISSIF: regex simplifiée (plus sensible)
    """
    sites = []
    
    try:
        if mode == "strict":
            pattern = (
                f"(?<=.{{{params['signalPeptideLength']}}})(?<!K|R)"
                f"(?:KK|KR|RR|RK)(?=[^RKILPVH]|(?<=KR)H|$)"
                f"(?=(?:(?!(?R)).){{{params['minCleavageSpacing']},}}|$))"
            )
        else:  # permissive
            pattern = (
                f"(?<=.{{{params['signalPeptideLength']}}})(?:KK|KR|RR|RK)"
                f"(?=[^RKILPVH]|(?<=KR)H|$)"
            )
        
        for match in re.finditer(pattern, seq):
            sites.append({
                'position': match.start() + 2,
                'motif': match.group(),
                'index': match.start()
            })
    
    except Exception as e:
        print(f"Erreur regex: {e}")
        return []
    
    return sites

# ==================== EXTRACTION PEPTIDES ====================

def extract_peptides(seq: str, cleavage_sites: List[dict], params: dict) -> List[dict]:
    """
    Extrait les peptides entre les sites de clivage
    Retourne liste de dicts (avant calcul bioactivité)
    """
    if len(cleavage_sites) < params['minCleavageSites']:
        return []
    
    peptides = []
    prev_position = params['signalPeptideLength']
    
    for site in cleavage_sites:
        current_pos = site['position']
        distance = current_pos - prev_position
        
        if distance >= params['minCleavageSpacing']:
            pep_seq = seq[prev_position:current_pos]
            
            if len(pep_seq) > 3:
                peptides.append({
                    'sequence': pep_seq,
                    'start': prev_position,
                    'end': current_pos,
                    'length': len(pep_seq),
                    'inRange': (5 <= len(pep_seq) <= 25),
                    'cleavageMotif': site['motif'],
                    'bioactivityScore': 0,  # À calculer
                    'bioactivitySource': ''  # À calculer
                })
            
            prev_position = current_pos
    
    # Dernier peptide
    if len(seq) - prev_position > 3:
        last_seq = seq[prev_position:]
        peptides.append({
            'sequence': last_seq,
            'start': prev_position,
            'end': len(seq),
            'length': len(last_seq),
            'inRange': (5 <= len(last_seq) <= 25),
            'cleavageMotif': 'END',
            'bioactivityScore': 0,  # À calculer
            'bioactivitySource': ''  # À calculer
        })
    
    return peptides

# ==================== ROUTES ====================

@app.get("/")
async def root():
    return {
        "message": "Peptide Predictor API v2.0",
        "version": "2.0.0",
        "features": ["REGEX detection", "Parallel bioactivity prediction", "PeptideRanker API integration"],
        "endpoints": {
            "analyze": "POST /analyze",
            "health": "GET /health",
            "docs": "/docs"
        }
    }

@app.get("/health")
async def health():
    return {
        "status": "ok",
        "timestamp": datetime.now().isoformat()
    }

@app.post("/analyze", response_model=AnalysisResponse)
async def analyze(request: AnalysisRequest):
    """
    Endpoint principal pour analyser une séquence protéique
    Utilise requêtes parallèles pour la bioactivité
    """
    try:
        # 1. Valider et nettoyer la séquence
        seq = request.sequence.strip()
        
        if seq.startswith('>'):
            seq = '\n'.join(seq.split('\n')[1:])
        
        seq = seq.replace('\n', '').replace(' ', '').upper()
        
        if not re.match(r'^[ACDEFGHIKLMNPQRSTVWY*]+$', seq):
            raise HTTPException(
                status_code=400,
                detail="Séquence invalide. Utilisez uniquement les codes standard des acides aminés."
            )
        
        if len(seq) < request.signalPeptideLength + 10:
            raise HTTPException(
                status_code=400,
                detail=f"Séquence trop courte. Min: {request.signalPeptideLength + 10} aa"
            )
        
        # 2. Préparer paramètres
        params = {
            'signalPeptideLength': request.signalPeptideLength,
            'minCleavageSites': request.minCleavageSites,
            'minCleavageSpacing': request.minCleavageSpacing
        }
        
        # 3. Détecter sites de clivage
        cleavage_sites = find_cleavage_sites(seq, request.mode, params)
        
        # 4. Extraire peptides
        peptides = extract_peptides(seq, cleavage_sites, params)
        
        # 5. Calculer bioactivité en PARALLÈLE
        async with aiohttp.ClientSession() as session:
            bioactivity_tasks = [
                calculate_bioactivity_async(p['sequence'], session)
                for p in peptides
            ]
            
            # Exécuter toutes les requêtes en parallèle
            bioactivity_results = await asyncio.gather(*bioactivity_tasks)
        
        # 6. Assigner les scores
        for peptide, (score, source) in zip(peptides, bioactivity_results):
            peptide['bioactivityScore'] = score
            peptide['bioactivitySource'] = source
        
        # 7. Trier par bioactivité
        peptides.sort(key=lambda x: x['bioactivityScore'], reverse=True)
        
        # 8. Top peptides
        top_peptides = peptides[:5]
        
        # 9. Créer résultats
        results = AnalysisResponse(
            sequenceLength=len(seq),
            cleavageSitesCount=len(cleavage_sites),
            peptides=[PeptideResult(**p) for p in peptides],
            peptidesInRange=sum(1 for p in peptides if p['inRange']),
            topPeptides=[PeptideResult(**p) for p in top_peptides]
        )
        
        return results
    
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Erreur serveur: {str(e)}"
        )

@app.get("/stats")
async def stats():
    return {
        "api": "Peptide Predictor",
        "version": "2.0.0",
        "modes": ["strict", "permissive"],
        "bioactivity": {
            "primary": "PeptideRanker API (parallel requests)",
            "fallback": "Heuristic scoring"
        },
        "default_params": {
            "signalPeptideLength": 20,
            "minCleavageSites": 4,
            "minCleavageSpacing": 5
        },
        "peptide_range_optimal": "5-25 aa"
    }

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)