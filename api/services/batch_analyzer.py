"""Service d'analyse batch pour plusieurs protéines"""
import asyncio
import aiohttp
from typing import List, Dict, Optional
from api.services import (
    SequenceValidator,
    CleavageDetector,
    PeptideExtractor,
    BioactivityPredictor,
    UniProtChecker,
    protein_db,
    ptm_detector
)


class BatchAnalyzer:
    """Analyseur batch pour plusieurs protéines"""
    
    @staticmethod
    async def analyze_single_protein(
        protein_id: str,
        mode: str,
        session: aiohttp.ClientSession,
        progress_callback: Optional[callable] = None
    ) -> Dict:
        """
        Analyse une seule protéine avec paramètres auto-recommandés
        
        Args:
            protein_id: UniProt accession
            mode: strict ou permissive
            session: Session aiohttp
            progress_callback: Fonction callback pour progression
        
        Returns:
            Dictionnaire avec résultats ou erreur
        """
        try:
            print(f"\n🔬 Starting analysis for: {protein_id}")
            
            # 1. Récupérer la protéine depuis UniProt
            protein = await protein_db.get_protein(protein_id, session)
            
            if not protein:
                return {
                    "status": "error",
                    "proteinId": protein_id,
                    "error": f"Protein {protein_id} not found or not secreted"
                }
            
            # 2. Extraire infos
            clean_seq = protein["sequence"]
            gene_name = protein["geneName"]
            protein_name = protein["proteinName"]
            accession = protein["accession"]
            
            protein_id_header = f"SP|{accession}|{gene_name}_HUMAN {protein_name}"
            
            # 3. Paramètres recommandés
            recommended_params = protein.get("recommendedParams", {
                "signalPeptideLength": 20,
                "minCleavageSites": 4,
                "minCleavageSpacing": 5,
                "maxPeptideLength": 100
            })
            
            signal_length = recommended_params["signalPeptideLength"]
            min_sites = recommended_params["minCleavageSites"]
            min_spacing = recommended_params["minCleavageSpacing"]
            max_length = recommended_params["maxPeptideLength"]
            
            print(f"📊 Using recommended params: signal={signal_length}, sites={min_sites}, spacing={min_spacing}")
            
            # 4. Valider séquence
            SequenceValidator.validate_characters(clean_seq)
            min_seq_length = signal_length + 10
            SequenceValidator.validate_length(clean_seq, min_seq_length)
            
            # 5. Détecter sites de clivage
            cleavage_sites = CleavageDetector.find_sites(
                sequence=clean_seq,
                mode=mode,
                signal_length=signal_length,
                min_spacing=min_spacing
            )
            
            # 6. Extraire peptides
            peptides = PeptideExtractor.extract(
                sequence=clean_seq,
                cleavage_sites=cleavage_sites,
                signal_length=signal_length,
                min_spacing=min_spacing,
                min_sites=min_sites,
                mode=mode
            )
            
            # 7. Filtrer par maxPeptideLength
            peptides_filtered = [
                p for p in peptides
                if p['length'] <= max_length
            ]
            
            print(f"📊 Peptides: {len(peptides)} → {len(peptides_filtered)} after filter")
            
            peptides = peptides_filtered
            
            # ⭐ Si 0 peptides, retourner succès mais avec liste vide
            if len(peptides) == 0:
                print(f"⚠️ No peptides found for {gene_name}")
                return {
                    "status": "success",
                    "proteinId": protein_id_header,
                    "geneName": gene_name,
                    "proteinName": protein_name,
                    "accession": accession,
                    "sequence": clean_seq,
                    "length": len(clean_seq),
                    "signalPeptideEnd": signal_length,
                    "fastaHeader": protein["fastaHeader"],
                    "sequenceLength": len(clean_seq),
                    "cleavageSitesCount": len(cleavage_sites),
                    "peptides": [],
                    "peptidesInRange": 0,
                    "topPeptides": [],
                    "cleavageSites": [
                        {"position": site.position, "motif": site.motif, "index": site.index}
                        for site in cleavage_sites
                    ],
                    "mode": mode
                }
            
            # ⭐ MODIFIÉ : 8. Calculer bioactivité AVEC CONTEXTE (parallèle)
            bioactivity_results = await BioactivityPredictor.predict_batch(
                peptides=[p['sequence'] for p in peptides],
                session=session,
                cleavage_motifs=[p['cleavageMotif'] for p in peptides],
                full_protein_sequence=clean_seq,
                peptide_end_positions=[p['end'] for p in peptides]
            )
            
            # 9. Vérifier UniProt (parallèle)
            uniprot_results = await UniProtChecker.check_batch(
                [p['sequence'] for p in peptides],
                session,
                protein_id=protein_id_header
            )
            
            # 10. Assigner scores bioactivité
            for peptide, (score, source) in zip(peptides, bioactivity_results):
                peptide['bioactivityScore'] = score
                peptide['bioactivitySource'] = source
            
            # 11. Assigner données UniProt
            for peptide, uniprot_data in zip(peptides, uniprot_results):
                peptide['uniprotStatus'] = uniprot_data['uniprotStatus']
                peptide['uniprotName'] = uniprot_data['uniprotName']
                peptide['uniprotNote'] = uniprot_data['uniprotNote']
                peptide['uniprotAccession'] = uniprot_data['uniprotAccession']
            
            # 12. Détecter PTMs
            print(f"🔬 Detecting PTMs for {len(peptides)} peptides...")
            for idx, peptide in enumerate(peptides, 1):
                try:
                    if not isinstance(peptide['start'], int) or not isinstance(peptide['end'], int):
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
            
            # 13. Trier par bioactivité
            peptides.sort(key=lambda x: x['bioactivityScore'], reverse=True)
            
            # 14. Top 5 peptides
            top_peptides = peptides[:5]
            
            # 15. Stats
            peptides_in_range = sum(1 for p in peptides if p['inRange'])
            
            # 16. Callback progression
            if progress_callback:
                await progress_callback(protein_id, "completed")
            
            print(f"✅ Analysis completed for {gene_name}: {len(peptides)} peptides")
            
            # 17. Retourner résultat
            return {
                "status": "success",
                "proteinId": protein_id_header,
                "geneName": gene_name,
                "proteinName": protein_name,
                "accession": accession,
                "sequence": clean_seq,
                "length": len(clean_seq),
                "signalPeptideEnd": signal_length,
                "fastaHeader": protein["fastaHeader"],
                "sequenceLength": len(clean_seq),
                "cleavageSitesCount": len(cleavage_sites),
                "peptides": peptides,
                "peptidesInRange": peptides_in_range,
                "topPeptides": top_peptides,
                "cleavageSites": [
                    {"position": site.position, "motif": site.motif, "index": site.index}
                    for site in cleavage_sites
                ],
                "mode": mode
            }
        
        except Exception as e:
            print(f"❌ Error analyzing {protein_id}: {e}")
            import traceback
            traceback.print_exc()
            
            if progress_callback:
                await progress_callback(protein_id, "error")
            
            return {
                "status": "error",
                "proteinId": protein_id,
                "error": str(e)
            }
    
    @staticmethod
    async def analyze_batch(
        protein_ids: List[str],
        mode: str,
        progress_callback: Optional[callable] = None
    ) -> Dict:
        """
        Analyse plusieurs protéines en séquentiel
        
        Args:
            protein_ids: Liste d'UniProt IDs
            mode: strict ou permissive
            progress_callback: Fonction callback pour progression
        
        Returns:
            Dictionnaire avec résultats batch
        """
        print(f"\n🚀 Starting batch analysis for {len(protein_ids)} proteins")
        print(f"📋 Proteins: {', '.join(protein_ids)}")
        
        # ⭐ ÉTAPE 0 : Dédupliquer les IDs
        unique_protein_ids = list(dict.fromkeys(protein_ids))  # Garde l'ordre
        
        if len(unique_protein_ids) < len(protein_ids):
            duplicates_count = len(protein_ids) - len(unique_protein_ids)
            print(f"⚠️ Removed {duplicates_count} duplicate(s). Analyzing {len(unique_protein_ids)} unique proteins.")
        
        results = []
        not_found = []
        
        async with aiohttp.ClientSession() as session:
            # Vérifier d'abord quelles protéines existent
            print("\n🔍 Step 1: Checking which proteins exist...")
            for protein_id in unique_protein_ids:
                protein = await protein_db.get_protein(protein_id, session)
                if not protein:
                    not_found.append(protein_id)
                    print(f"❌ Protein not found: {protein_id}")
            
            # Filtrer les protéines trouvées
            valid_protein_ids = [pid for pid in unique_protein_ids if pid not in not_found]
            
            print(f"\n✅ Found {len(valid_protein_ids)} proteins")
            if not_found:
                print(f"❌ Not found: {', '.join(not_found)}")
            
            # Analyser séquentiellement
            print("\n🔬 Step 2: Analyzing proteins sequentially...")
            for idx, protein_id in enumerate(valid_protein_ids, 1):
                print(f"\n─── Protein {idx}/{len(valid_protein_ids)} ───")
                
                if progress_callback:
                    await progress_callback(protein_id, "analyzing", idx, len(valid_protein_ids))
                
                result = await BatchAnalyzer.analyze_single_protein(
                    protein_id=protein_id,
                    mode=mode,
                    session=session,
                    progress_callback=None
                )
                
                results.append(result)
                
                # Petite pause entre analyses
                await asyncio.sleep(0.5)
        
        # Compter succès/erreurs
        successful = sum(1 for r in results if r.get("status") == "success")
        failed = sum(1 for r in results if r.get("status") == "error")
        
        print(f"\n✅ Batch analysis completed:")
        print(f"   - Total: {len(protein_ids)} proteins requested")
        print(f"   - Unique: {len(unique_protein_ids)} proteins")
        print(f"   - Successful: {successful}")
        print(f"   - Failed: {failed}")
        print(f"   - Not found: {len(not_found)}")
        
        return {
            "totalProteins": len(protein_ids),
            "successfulProteins": successful,
            "failedProteins": failed,
            "results": results,
            "notFound": not_found,
            "mode": mode
        }


# Instance globale
batch_analyzer = BatchAnalyzer()