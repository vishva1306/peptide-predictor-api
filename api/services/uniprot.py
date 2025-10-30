"""Vérification des peptides connus dans UniProt - Version 3 statuts"""
import aiohttp
import asyncio
from typing import Dict, Optional, List

class UniProtChecker:
    """Vérificateur de peptides connus dans UniProt"""
    
    BASE_URL = "https://rest.uniprot.org/uniprotkb"
    
    @staticmethod
    async def get_protein_features(
        protein_id: str,
        session: aiohttp.ClientSession
    ) -> List[Dict]:
        """Récupère tous les peptides annotés d'une protéine"""
        
        print(f"\n🔍 Récupération features pour protéine : {protein_id}")
        
        # Nettoyer l'ID (enlever le prefix sp|tr| si présent)
        clean_id = protein_id
        if '|' in clean_id:
            parts = clean_id.split('|')
            clean_id = parts[1] if len(parts) > 1 else parts[0]
        
        url = f"{UniProtChecker.BASE_URL}/{clean_id}"
        
        params = {
            "format": "json",
            "fields": "sequence,ft_peptide,ft_propep"
        }
        
        try:
            async with session.get(
                url,
                params=params,
                timeout=aiohttp.ClientTimeout(total=15)
            ) as response:
                
                print(f"✅ Status code : {response.status}")
                
                if response.status != 200:
                    print(f"❌ Erreur HTTP {response.status}")
                    return []
                
                data = await response.json()
                
                # Récupérer la séquence complète de la protéine
                full_sequence = data.get("sequence", {}).get("value", "")
                
                if not full_sequence:
                    print(f"❌ Pas de séquence trouvée")
                    return []
                
                print(f"📊 Séquence protéine : {len(full_sequence)} aa")
                
                # Extraire les features
                features = data.get("features", [])
                peptide_features = []
                
                for feat in features:
                    ftype = feat.get("type", "")
                    
                    if ftype in ["Peptide", "Propeptide"]:
                        location = feat.get("location", {})
                        start = location.get("start", {}).get("value", 0)
                        end = location.get("end", {}).get("value", 0)
                        
                        if start > 0 and end > 0 and start <= len(full_sequence):
                            # Extraire la séquence du peptide (attention: indices 1-based)
                            peptide_seq = full_sequence[start-1:end]
                            
                            peptide_features.append({
                                "type": ftype,
                                "description": feat.get("description", "Peptide connu"),
                                "start": start,
                                "end": end,
                                "sequence": peptide_seq
                            })
                
                print(f"✅ {len(peptide_features)} peptides annotés trouvés")
                return peptide_features
        
        except asyncio.TimeoutError:
            print(f"⏱️ Timeout UniProt")
            return []
        except Exception as e:
            print(f"❌ Erreur UniProt: {e}")
            return []
    
    @staticmethod
    def find_matching_peptide(
        peptide_seq: str,
        annotated_peptides: List[Dict]
    ) -> Optional[Dict]:
        """
        Cherche un match avec 3 niveaux de précision
        
        Returns:
            {
                "match_type": "exact" | "partial" | None,
                "description": "Nom du peptide",
                "note": "Information additionnelle sur le type de fragment"
            }
        """
        
        for annotated in annotated_peptides:
            annotated_seq = annotated["sequence"]
            
            # 1. EXACT MATCH ✅
            if peptide_seq == annotated_seq:
                print(f"✅ EXACT MATCH : {annotated['description']}")
                return {
                    "match_type": "exact",
                    "description": annotated['description'],
                    "note": None
                }
            
            # 2. PARTIAL MATCH - Fragment (peptide détecté DANS peptide annoté) ⚠️
            if peptide_seq in annotated_seq:
                # Déterminer quelle partie du peptide annoté
                start_pos = annotated_seq.index(peptide_seq)
                end_pos = start_pos + len(peptide_seq)
                
                if start_pos == 0:
                    fragment_type = "N-terminal fragment"
                elif end_pos == len(annotated_seq):
                    fragment_type = "C-terminal fragment"
                else:
                    fragment_type = "Internal fragment"
                
                print(f"⚠️ PARTIAL MATCH (fragment) : {fragment_type} of {annotated['description']}")
                
                return {
                    "match_type": "partial",
                    "description": annotated['description'],
                    "note": fragment_type
                }
            
            # 3. PARTIAL MATCH - Extension (peptide annoté DANS peptide détecté) ⚠️
            if annotated_seq in peptide_seq:
                print(f"⚠️ PARTIAL MATCH (extension) : Extended form of {annotated['description']}")
                
                return {
                    "match_type": "partial",
                    "description": annotated['description'],
                    "note": "Extended form"
                }
        
        return None
    
    @classmethod
    async def check_batch(
        cls,
        peptides: List[str],
        session: aiohttp.ClientSession,
        protein_id: Optional[str] = None
    ) -> List[Dict]:
        """
        Vérifie un batch de peptides contre les annotations UniProt
        
        Returns:
            Liste de dictionnaires avec :
            - uniprotStatus: "exact" | "partial" | "unknown"
            - uniprotName: Nom du peptide (si trouvé)
            - uniprotNote: Note additionnelle (type de fragment)
            - uniprotAccession: Accession UniProt
        """
        
        print(f"\n🚀 Début vérification UniProt pour {len(peptides)} peptides...")
        
        results = []
        
        # Si pas d'ID protéine fourni, impossible de vérifier
        if not protein_id or protein_id == "N/A":
            print(f"⚠️ Pas d'ID protéine fourni - skip vérification UniProt")
            return [
                {
                    "uniprotStatus": "unknown",
                    "uniprotName": None,
                    "uniprotNote": None,
                    "uniprotAccession": None
                }
                for _ in peptides
            ]
        
        # Récupérer tous les peptides annotés de la protéine (1 seule requête)
        annotated_peptides = await cls.get_protein_features(protein_id, session)
        
        if not annotated_peptides:
            print(f"⚠️ Aucun peptide annoté trouvé pour {protein_id}")
            return [
                {
                    "uniprotStatus": "unknown",
                    "uniprotName": None,
                    "uniprotNote": None,
                    "uniprotAccession": None
                }
                for _ in peptides
            ]
        
        # Extraire l'accession UniProt propre
        clean_accession = protein_id
        if '|' in clean_accession:
            parts = clean_accession.split('|')
            clean_accession = parts[1] if len(parts) > 1 else parts[0]
        
        # Comparer chaque peptide détecté avec les peptides annotés
        for i, peptide_seq in enumerate(peptides, 1):
            print(f"\n--- Peptide {i}/{len(peptides)} : {peptide_seq[:30]}... ---")
            
            match = cls.find_matching_peptide(peptide_seq, annotated_peptides)
            
            if match:
                results.append({
                    "uniprotStatus": match["match_type"],
                    "uniprotName": match["description"],
                    "uniprotNote": match["note"],
                    "uniprotAccession": clean_accession
                })
            else:
                print(f"❌ Aucun match trouvé")
                results.append({
                    "uniprotStatus": "unknown",
                    "uniprotName": None,
                    "uniprotNote": None,
                    "uniprotAccession": None
                })
        
        exact_count = sum(1 for r in results if r['uniprotStatus'] == 'exact')
        partial_count = sum(1 for r in results if r['uniprotStatus'] == 'partial')
        unknown_count = sum(1 for r in results if r['uniprotStatus'] == 'unknown')
        
        print(f"\n✅ Vérification terminée :")
        print(f"   - {exact_count} exact matches")
        print(f"   - {partial_count} partial matches")
        print(f"   - {unknown_count} unknown")
        
        return results