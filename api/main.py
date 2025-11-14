"""FastAPI application principale"""
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field, field_validator
from typing import Optional, Union, List
import asyncio
import aiohttp

from api.services import (
    SequenceValidator,
    CleavageDetector,
    PeptideExtractor,
    BioactivityPredictor,
    UniProtChecker,
    protein_db,
    ptm_detector,
    batch_analyzer,
    fasta_parser
)

app = FastAPI(title="Peptide Predictor API")

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
    """Modèle pour la requête d'analyse"""
    proteinId: Optional[Union[str, List[str]]] = None
    fastaSequence: Optional[str] = None
    fastaHeader: Optional[str] = None
    mode: str = Field(default="permissive")
    signalPeptideLength: int = Field(default=20, ge=0, le=100)
    minCleavageSites: int = Field(default=4, ge=1, le=10)
    minCleavageSpacing: int = Field(default=5, ge=1, le=20)
    maxPeptideLength: int = Field(default=100, ge=10, le=500)

    @field_validator('mode')
    @classmethod
    def validate_mode(cls, v):
        if v not in ['strict', 'permissive', 'ultra-permissive']:
            raise ValueError('Mode must be "strict" ,"permissive", or "ultra-permissive"')
        return v
    
    @field_validator('fastaSequence')
    @classmethod
    def validate_input(cls, v, info):
        # Au moins un des deux doit être fourni
        if not v and not info.data.get('proteinId'):
            raise ValueError('Either proteinId or fastaSequence must be provided')
        return v


# ==================== ENDPOINTS ====================

@app.get("/")
async def root():
    return {
        "message": "Peptide Predictor API",
        "version": "2.3.0",
        "endpoints": {
            "/analyze": "POST - Analyze protein(s) or FASTA sequence",
            "/api/proteins/search": "GET - Search proteins",
            "/api/proteins/{accession}": "GET - Get protein details"
        }
    }


@app.post("/analyze")
async def analyze_protein(request: AnalysisRequest):
    """
    Analyse une ou plusieurs protéines depuis UniProt OU une séquence FASTA
    
    Modes:
    1. Single UniProt: proteinId = "P01189"
    2. Batch UniProt: proteinId = ["P01189", "P01308", "Q9UBU3"]
    3. FASTA: fastaSequence = "MALWMR..." (avec fastaHeader optionnel)
    """
    
    # ==================== MODE FASTA ====================
    if request.fastaSequence:
        print("\n🧬 FASTA MODE DETECTED")
        
        try:
            # Parser la séquence FASTA
            fasta_data = fasta_parser.parse(request.fastaSequence)
            sequence = fasta_data['sequence']
            
            print(f"📊 Sequence length: {len(sequence)} aa")
            print(f"📋 Header: {fasta_data['header']}")
            print(f"🆔 ID: {fasta_data['id']}")
            print(f"📝 Name: {fasta_data['name']}")
            
            # Valider la séquence
            is_valid, error_msg = fasta_parser.validate_sequence(sequence)
            if not is_valid:
                raise HTTPException(status_code=400, detail=error_msg)
            
            # Détection des sites de clivage
            cleavage_sites = CleavageDetector.find_sites(
                sequence=sequence,
                mode=request.mode,
                signal_length=request.signalPeptideLength,
                min_spacing=request.minCleavageSpacing
            )
            
            print(f"✂️ Found {len(cleavage_sites)} cleavage sites")
            
            # Extraction des peptides
            peptides = PeptideExtractor.extract(
                sequence=sequence,
                cleavage_sites=cleavage_sites,
                signal_length=request.signalPeptideLength,
                min_spacing=request.minCleavageSpacing,
                min_sites=request.minCleavageSites,
                mode=request.mode
            )
            
            # Filtrer par longueur max
            peptides = [p for p in peptides if p['length'] <= request.maxPeptideLength]
            
            print(f"🧬 Extracted {len(peptides)} peptides")
            
            if len(peptides) == 0:
                return {
                    "sequenceLength": len(sequence),
                    "cleavageSitesCount": len(cleavage_sites),
                    "peptides": [],
                    "peptidesInRange": 0,
                    "proteinId": fasta_data['id'] or "Custom",
                    "geneName": fasta_data['name'] or "FASTA",
                    "proteinName": fasta_data['name'] or "Custom FASTA Sequence",
                    "mode": request.mode,
                    "isFasta": True,
                    "fastaHeader": fasta_data['header']
                }
            
            # ⭐ MODIFIÉ : Prédiction bioactivité AVEC CONTEXTE (parallèle)
            async with aiohttp.ClientSession() as session:
                bioactivity_results = await BioactivityPredictor.predict_batch(
                    peptides=[p['sequence'] for p in peptides],
                    session=session,
                    cleavage_motifs=[p['cleavageMotif'] for p in peptides],
                    full_protein_sequence=sequence,
                    peptide_end_positions=[p['end'] for p in peptides]
                )
            
            # Assigner scores
            for peptide, (score, source) in zip(peptides, bioactivity_results):
                peptide['bioactivityScore'] = score
                peptide['bioactivitySource'] = source
                # Pas de vérification UniProt pour FASTA
                peptide['uniprotStatus'] = 'n/a'
                peptide['uniprotName'] = None
                peptide['uniprotNote'] = None
                peptide['uniprotAccession'] = None
            
            # Détection PTMs
            print(f"🔬 Detecting PTMs for {len(peptides)} peptides...")
            for idx, peptide in enumerate(peptides, 1):
                try:
                    detected_ptms = ptm_detector.detect_all_ptms(
                        peptide_sequence=peptide['sequence'],
                        full_protein_sequence=sequence,
                        peptide_start=peptide['start'],
                        peptide_end=peptide['end']
                    )
                    
                    peptide['ptms'] = detected_ptms
                    
                    if detected_ptms:
                        peptide['modifiedSequence'] = ptm_detector.generate_modified_sequence(
                            peptide['sequence'],
                            detected_ptms
                        )
                    else:
                        peptide['modifiedSequence'] = None
                        
                except Exception as e:
                    print(f"❌ PTM detection error for peptide {idx}: {e}")
                    peptide['ptms'] = []
                    peptide['modifiedSequence'] = None
            
            # Trier par bioactivité
            peptides.sort(key=lambda x: x['bioactivityScore'], reverse=True)
            
            # Stats
            peptides_in_range = sum(1 for p in peptides if p['inRange'])
            
            return {
                "sequenceLength": len(sequence),
                "cleavageSitesCount": len(cleavage_sites),
                "peptides": peptides,
                "peptidesInRange": peptides_in_range,
                "proteinId": fasta_data['id'] or "Custom",
                "geneName": fasta_data['name'] or "FASTA",
                "proteinName": fasta_data['name'] or "Custom FASTA Sequence",
                "cleavageSites": [
                    {"position": site.position, "motif": site.motif, "index": site.index}
                    for site in cleavage_sites
                ],
                "mode": request.mode,
                "isFasta": True,
                "fastaHeader": fasta_data['header']
            }
            
        except Exception as e:
            print(f"❌ FASTA analysis error: {e}")
            import traceback
            traceback.print_exc()
            raise HTTPException(status_code=500, detail=str(e))
    
    # ==================== MODE UNIPROT (BATCH) ====================
    elif isinstance(request.proteinId, list):
        print(f"\n📦 BATCH MODE: {len(request.proteinId)} proteins")
        
        result = await batch_analyzer.analyze_batch(
            protein_ids=request.proteinId,
            mode=request.mode
        )
        
        return result
    
    # ==================== MODE UNIPROT (SINGLE) ====================
    elif request.proteinId:
        print(f"\n🔬 SINGLE MODE: {request.proteinId}")
        
        async with aiohttp.ClientSession() as session:
            # Récupérer la protéine
            protein = await protein_db.get_protein(request.proteinId, session)
            
            if not protein:
                raise HTTPException(
                    status_code=404,
                    detail=f"Protein {request.proteinId} not found or not secreted"
                )
            
            clean_seq = protein["sequence"]
            gene_name = protein["geneName"]
            protein_name = protein["proteinName"]
            accession = protein["accession"]
            
            protein_id = f"SP|{accession}|{gene_name}_HUMAN {protein_name}"
            
            # Validation
            SequenceValidator.validate_characters(clean_seq)
            min_seq_length = request.signalPeptideLength + 10
            SequenceValidator.validate_length(clean_seq, min_seq_length)
            
            # Détection sites de clivage
            cleavage_sites = CleavageDetector.find_sites(
                sequence=clean_seq,
                mode=request.mode,
                signal_length=request.signalPeptideLength,
                min_spacing=request.minCleavageSpacing
            )
            
            # Extraction peptides
            peptides = PeptideExtractor.extract(
                sequence=clean_seq,
                cleavage_sites=cleavage_sites,
                signal_length=request.signalPeptideLength,
                min_spacing=request.minCleavageSpacing,
                min_sites=request.minCleavageSites,
                mode=request.mode
            )
            
            # Filtrer par longueur max
            peptides = [p for p in peptides if p['length'] <= request.maxPeptideLength]
            
            if len(peptides) == 0:
                return {
                    "sequenceLength": len(clean_seq),
                    "cleavageSitesCount": len(cleavage_sites),
                    "peptides": [],
                    "peptidesInRange": 0,
                    "proteinId": protein_id,
                    "geneName": gene_name,
                    "proteinName": protein_name,
                    "mode": request.mode
                }
            
            # ⭐ MODIFIÉ : Bioactivité AVEC CONTEXTE
            bioactivity_results = await BioactivityPredictor.predict_batch(
                peptides=[p['sequence'] for p in peptides],
                session=session,
                cleavage_motifs=[p['cleavageMotif'] for p in peptides],
                full_protein_sequence=clean_seq,
                peptide_end_positions=[p['end'] for p in peptides]
            )
            
            # UniProt check
            uniprot_results = await UniProtChecker.check_batch(
                [p['sequence'] for p in peptides],
                session,
                protein_id=protein_id
            )
            
            # Assigner
            for peptide, (score, source), uniprot_data in zip(peptides, bioactivity_results, uniprot_results):
                peptide['bioactivityScore'] = score
                peptide['bioactivitySource'] = source
                peptide['uniprotStatus'] = uniprot_data['uniprotStatus']
                peptide['uniprotName'] = uniprot_data['uniprotName']
                peptide['uniprotNote'] = uniprot_data['uniprotNote']
                peptide['uniprotAccession'] = uniprot_data['uniprotAccession']
            
            # PTMs
            print(f"🔬 Detecting PTMs for {len(peptides)} peptides...")
            for idx, peptide in enumerate(peptides, 1):
                try:
                    detected_ptms = ptm_detector.detect_all_ptms(
                        peptide_sequence=peptide['sequence'],
                        full_protein_sequence=clean_seq,
                        peptide_start=peptide['start'],
                        peptide_end=peptide['end']
                    )
                    
                    peptide['ptms'] = detected_ptms
                    
                    if detected_ptms:
                        peptide['modifiedSequence'] = ptm_detector.generate_modified_sequence(
                            peptide['sequence'],
                            detected_ptms
                        )
                    else:
                        peptide['modifiedSequence'] = None
                        
                except Exception as e:
                    print(f"❌ PTM detection error for peptide {idx}: {e}")
                    peptide['ptms'] = []
                    peptide['modifiedSequence'] = None
            
            # Trier
            peptides.sort(key=lambda x: x['bioactivityScore'], reverse=True)
            
            peptides_in_range = sum(1 for p in peptides if p['inRange'])
            
            return {
                "sequenceLength": len(clean_seq),
                "cleavageSitesCount": len(cleavage_sites),
                "peptides": peptides,
                "peptidesInRange": peptides_in_range,
                "proteinId": protein_id,
                "geneName": gene_name,
                "proteinName": protein_name,
                "cleavageSites": [
                    {"position": site.position, "motif": site.motif, "index": site.index}
                    for site in cleavage_sites
                ],
                "mode": request.mode,
                "isFasta": False
            }
    
    else:
        raise HTTPException(
            status_code=400,
            detail="Either proteinId or fastaSequence must be provided"
        )


@app.get("/api/proteins/search")
async def search_proteins(q: str, type: str = "gene_name", limit: int = 10):
    """Recherche de protéines dans UniProt"""
    async with aiohttp.ClientSession() as session:
        results = await protein_db.search_proteins(q, type, limit, session)
        return results


@app.get("/api/proteins/{accession}")
async def get_protein(accession: str):
    """Récupère les détails d'une protéine"""
    async with aiohttp.ClientSession() as session:
        protein = await protein_db.get_protein(accession, session)
        if not protein:
            raise HTTPException(status_code=404, detail="Protein not found")
        return protein


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000, reload=True)