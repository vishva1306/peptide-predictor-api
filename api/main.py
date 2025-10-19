from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import List, Optional
import re
import json
from datetime import datetime

app = FastAPI(
    title="Peptide Predictor API",
    version="1.0.0",
    description="API pour prédire les peptides bioactifs"
)

# CORS - permet au frontend React de communiquer
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

class AnalysisResponse(BaseModel):
    sequenceLength: int
    cleavageSitesCount: int
    peptides: List[PeptideResult]
    peptidesInRange: int
    topPeptides: List[PeptideResult]

# ==================== BIOACTIVITÉ ====================

def calculate_bioactivity(peptide: str) -> float:
    """
    Calcule le score de bioactivité (0-100) basé sur propriétés physicochimiques
    Basé sur les critères du papier Nature Coassolo et al. 2025
    """
    score = 0.0
    
    if len(peptide) == 0:
        return 0.0
    
    # 1. Hydrophobicité (30%)
    hydrophobic = ['A', 'L', 'I', 'V', 'M', 'F', 'W', 'P']
    hydro_count = sum(1 for aa in peptide if aa in hydrophobic)
    hydro_ratio = hydro_count / len(peptide)
    score += hydro_ratio * 30
    
    # 2. Charge positive/négative (20%)
    positive = sum(1 for aa in peptide if aa in ['K', 'R', 'H'])
    negative = sum(1 for aa in peptide if aa in ['D', 'E'])
    
    if positive > 0:
        score += 10
    if negative > 0:
        score += 10
    
    # 3. Longueur gamme optimal 5-25 aa (35%)
    if 5 <= len(peptide) <= 25:
        score += 35
    elif len(peptide) < 5:
        score -= 10
    elif len(peptide) > 100:
        score -= 15
    
    # 4. Stabilité et structure (15%)
    # Cystéine pour ponts disulfure
    if 'C' in peptide:
        score += 8
    
    # Trop de proline = déstabilise la structure
    if peptide.count('P') <= 2:
        score += 7
    else:
        score -= 5
    
    # Diversité en acides aminés (bonne pour activité)
    unique_aa = len(set(peptide))
    if unique_aa >= 6:
        score += 5
    
    return min(max(score, 0), 100)

# ==================== DÉTECTION SITES CLIVAGE ====================

def find_cleavage_sites(seq: str, mode: str, params: dict) -> List[dict]:
    """
    Détecte les sites de clivage PCSK1/3 en utilisant regex
    Mode STRICT: regex complète du papier Nature
    Mode PERMISSIF: regex simplifiée (plus sensible)
    """
    sites = []
    
    try:
        if mode == "strict":
            # Regex stricte: lookbehind (?<!K|R) + vérification espacement
            pattern = (
                f"(?<=.{{{params['signalPeptideLength']}}})(?<!K|R)"
                f"(?:KK|KR|RR|RK)(?=[^RKILPVH]|(?<=KR)H|$)"
                f"(?=(?:(?!(?R)).){{{params['minCleavageSpacing']},}}|$)"
            )
        else:  # permissive
            # Regex permissive: pas de lookbehind strict, pas de vérification espacement
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

def extract_peptides(seq: str, cleavage_sites: List[dict], params: dict) -> List[PeptideResult]:
    """
    Extrait les peptides entre les sites de clivage
    """
    if len(cleavage_sites) < params['minCleavageSites']:
        return []
    
    peptides = []
    prev_pos = params['signalPeptideLength']
    
    for site in cleavage_sites:
        current_pos = site['position']
        distance = current_pos - prev_pos
        
        # Vérifier espacement minimum
        if distance >= params['minCleavageSpacing']:
            pep_seq = seq[prev_pos:current_pos]
            
            # Inclure peptides > 3 aa
            if len(pep_seq) > 3:
                peptides.append(PeptideResult(
                    sequence=pep_seq,
                    start=prev_pos,
                    end=current_pos,
                    length=len(pep_seq),
                    inRange=(5 <= len(pep_seq) <= 25),
                    cleavageMotif=site['motif'],
                    bioactivityScore=calculate_bioactivity(pep_seq)
                ))
            
            prev_pos = current_pos
    
    # Dernier peptide jusqu'à fin de séquence
    if len(seq) - prev_pos > 3:
        last_seq = seq[prev_pos:]
        peptides.append(PeptideResult(
            sequence=last_seq,
            start=prev_pos,
            end=len(seq),
            length=len(last_seq),
            inRange=(5 <= len(last_seq) <= 25),
            cleavageMotif="END",
            bioactivityScore=calculate_bioactivity(last_seq)
        ))
    
    return peptides

# ==================== ROUTES ====================

@app.get("/")
async def root():
    """Route de bienvenue"""
    return {
        "message": "Peptide Predictor API",
        "version": "1.0.0",
        "endpoints": {
            "analyze": "POST /analyze",
            "health": "GET /health",
            "docs": "/docs"
        }
    }

@app.get("/health")
async def health():
    """Vérifier que le serveur fonctionne"""
    return {
        "status": "ok",
        "timestamp": datetime.now().isoformat()
    }

@app.post("/analyze", response_model=AnalysisResponse)
async def analyze(request: AnalysisRequest):
    """
    Endpoint principal pour analyser une séquence protéique
    
    Paramètres:
    - sequence: séquence FASTA ou protéique
    - mode: "strict" ou "permissive"
    - signalPeptideLength: longueur du peptide signal (défaut: 20)
    - minCleavageSites: min sites de clivage (défaut: 4)
    - minCleavageSpacing: min espacement (défaut: 5)
    """
    try:
        # 1. Valider et nettoyer la séquence
        seq = request.sequence.strip()
        
        # Supprimer en-têtes FASTA
        if seq.startswith('>'):
            seq = '\n'.join(seq.split('\n')[1:])
        
        # Supprimer espaces et convertir en majuscules
        seq = seq.replace('\n', '').replace(' ', '').upper()
        
        # Valider caractères
        if not re.match(r'^[ACDEFGHIKLMNPQRSTVWY*]+$', seq):
            raise HTTPException(
                status_code=400,
                detail="Séquence invalide. Utilisez uniquement les codes standard des acides aminés."
            )
        
        # Vérifier longueur
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
        
        # 5. Trier par bioactivité
        peptides.sort(key=lambda x: x.bioactivityScore, reverse=True)
        
        # 6. Top peptides (top 5)
        top_peptides = peptides[:5]
        
        return AnalysisResponse(
            sequenceLength=len(seq),
            cleavageSitesCount=len(cleavage_sites),
            peptides=peptides,
            peptidesInRange=sum(1 for p in peptides if p.inRange),
            topPeptides=top_peptides
        )
    
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Erreur serveur: {str(e)}"
        )

@app.get("/stats")
async def stats():
    """Retourner des statistiques sur l'API"""
    return {
        "api": "Peptide Predictor",
        "version": "1.0.0",
        "modes": ["strict", "permissive"],
        "default_params": {
            "signalPeptideLength": 20,
            "minCleavageSites": 4,
            "minCleavageSpacing": 5
        },
        "peptide_range_optimal": "5-25 aa"
    }

# Pour développement local
if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)