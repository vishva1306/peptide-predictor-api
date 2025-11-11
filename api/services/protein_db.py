"""Service de recherche et récupération de protéines depuis UniProt avec cache"""
import aiohttp
import asyncio
from typing import Dict, List, Optional
from datetime import datetime, timedelta
import re

class ProteinDatabase:
    """
    Gestionnaire de protéines sécrétées humaines
    - Recherche par gene name ou UniProt ID
    - Cache 24h pour performance
    - Calcul automatique des paramètres recommandés
    """
    
    BASE_URL = "https://rest.uniprot.org/uniprotkb"
    CACHE_DURATION = timedelta(hours=24)
    
    def __init__(self):
        self.cache: Dict[str, Dict] = {}
    
    def _is_uniprot_id(self, query: str) -> bool:
        """Détecte si la query est un UniProt ID (format: P01189)"""
        return bool(re.match(r'^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$', query.upper()))
    
    def _get_cache(self, key: str) -> Optional[Dict]:
        """Récupère du cache si valide (<24h)"""
        if key in self.cache:
            cached = self.cache[key]
            if datetime.now() - cached['timestamp'] < self.CACHE_DURATION:
                print(f"✅ Cache HIT: {key}")
                return cached['data']
        return None
    
    def _set_cache(self, key: str, data: Dict):
        """Stocke en cache"""
        self.cache[key] = {
            'data': data,
            'timestamp': datetime.now()
        }
        print(f"💾 Cache SET: {key}")
    
    async def search_proteins(
        self,
        query: str,
        search_type: str,
        limit: int,
        session: aiohttp.ClientSession
    ) -> List[Dict]:
        """
        Recherche de protéines sécrétées humaines par gene name ou ID
        
        Args:
            query: Gene name (ex: "POMC") ou UniProt ID (ex: "P01189")
            search_type: "gene_name" ou "accession"
            limit: Nombre max de résultats
            session: Session aiohttp
        
        Returns:
            Liste de protéines matchant la requête
        """
        
        cache_key = f"search_{search_type}_{query.lower()}"
        
        # Check cache
        cached = self._get_cache(cache_key)
        if cached:
            return cached
        
        print(f"\n🔍 Searching proteins for: {query} (type: {search_type})")
        
        # Construire la query UniProt
        if search_type == "accession":
            uniprot_query = f"(accession:{query.upper()})"
        else:
            uniprot_query = f"(gene:{query.upper()})"
        
        uniprot_query += " AND (organism_id:9606) AND (reviewed:true) AND (cc_subcellular_location:Secreted)"
        
        print(f"📝 UniProt query: {uniprot_query}")
        
        url = f"{self.BASE_URL}/search"
        params = {
            "query": uniprot_query,
            "format": "json",
            "size": limit,
            "fields": "accession,gene_names,protein_name,sequence,length,ft_signal,ft_peptide,ft_propep"
        }
        
        try:
            async with session.get(
                url,
                params=params,
                timeout=aiohttp.ClientTimeout(total=15)
            ) as response:
                
                print(f"📡 UniProt response status: {response.status}")
                
                if response.status != 200:
                    error_text = await response.text()
                    print(f"❌ UniProt error {response.status}: {error_text[:200]}")
                    return []
                
                data = await response.json()
                results = data.get("results", [])
                
                print(f"✅ Found {len(results)} proteins")
                
                # Parser les résultats
                proteins = []
                for entry in results:
                    # ✅ Vérifier que c'est bien un dict
                    if isinstance(entry, dict):
                        protein = self._parse_protein_entry(entry)
                        if protein:
                            proteins.append(protein)
                    else:
                        print(f"⚠️ Skipping non-dict entry: {type(entry)}")
                
                # Cache les résultats
                self._set_cache(cache_key, proteins)
                
                return proteins
        
        except asyncio.TimeoutError:
            print(f"⏱️ Timeout searching proteins")
            return []
        except Exception as e:
            print(f"❌ Error searching proteins: {e}")
            import traceback
            traceback.print_exc()
            return []
    
    async def get_protein(
        self,
        accession: str,
        session: aiohttp.ClientSession
    ) -> Optional[Dict]:
        """
        Récupère les détails complets d'une protéine
        
        Args:
            accession: UniProt ID (ex: "P01189")
            session: Session aiohttp
        
        Returns:
            Protéine complète avec paramètres recommandés
        """
        
        cache_key = f"protein_{accession}"
        
        # Check cache
        cached = self._get_cache(cache_key)
        if cached:
            return cached
        
        print(f"\n🔍 Fetching protein: {accession}")
        
        url = f"{self.BASE_URL}/{accession}"
        params = {
            "format": "json",
            "fields": "accession,gene_names,protein_name,sequence,length,ft_signal,ft_peptide,ft_propep"
        }
        
        try:
            async with session.get(
                url,
                params=params,
                timeout=aiohttp.ClientTimeout(total=15)
            ) as response:
                
                if response.status != 200:
                    print(f"❌ Protein not found: {accession}")
                    return None
                
                entry = await response.json()
                protein = self._parse_protein_entry(entry, full_details=True)
                
                if protein:
                    # Cache le résultat
                    self._set_cache(cache_key, protein)
                
                return protein
        
        except asyncio.TimeoutError:
            print(f"⏱️ Timeout fetching protein")
            return None
        except Exception as e:
            print(f"❌ Error fetching protein: {e}")
            return None
    
    def _parse_protein_entry(self, entry: Dict, full_details: bool = False) -> Optional[Dict]:
        """Parse une entrée UniProt en format standardisé"""
        
        try:
            # Accession
            accession = entry.get("primaryAccession")
            
            # Gene name
            genes = entry.get("genes", [])
            gene_name = None
            if genes:
                gene_name = genes[0].get("geneName", {}).get("value")
            
            # Protein name
            protein_desc = entry.get("proteinDescription", {})
            recommended = protein_desc.get("recommendedName", {})
            protein_name = recommended.get("fullName", {}).get("value", "Unknown protein")
            
            # Sequence
            sequence_data = entry.get("sequence", {})
            sequence = sequence_data.get("value", "")
            length = sequence_data.get("length", 0)
            
            if not accession or not sequence:
                return None
            
            # Signal peptide
            signal_end = 20
            features = entry.get("features", [])
            for feat in features:
                if feat.get("type") == "Signal":
                    location = feat.get("location", {})
                    signal_end = location.get("end", {}).get("value", 20)
                    break
            
            # Peptides annotés (si full_details)
            annotated_peptides = []
            if full_details:
                for feat in features:
                    ftype = feat.get("type", "")
                    if ftype in ["Peptide", "Propeptide"]:
                        location = feat.get("location", {})
                        start = location.get("start", {}).get("value", 0)
                        end = location.get("end", {}).get("value", 0)
                        
                        if start > 0 and end > 0 and start <= len(sequence):
                            peptide_seq = sequence[start-1:end]
                            annotated_peptides.append({
                                "name": feat.get("description", "Peptide"),
                                "start": start,
                                "end": end,
                                "sequence": peptide_seq
                            })
            
            # Calculer paramètres recommandés
            recommended_params = self.calculate_recommended_params(
                length=length,
                signal_end=signal_end,
                num_peptides=len(annotated_peptides)
            )
            
            # Header FASTA
            fasta_header = f">sp|{accession}|{gene_name or 'UNKN'}_HUMAN {protein_name}"
            
            protein = {
                "accession": accession,
                "geneName": gene_name or "Unknown",
                "proteinName": protein_name,
                "length": length,
                "sequence": sequence,
                "signalPeptideEnd": signal_end,
                "recommendedParams": recommended_params,
                "fastaHeader": fasta_header
            }
            
            if full_details:
                protein["annotatedPeptides"] = annotated_peptides
            
            return protein
        
        except Exception as e:
            print(f"❌ Error parsing protein entry: {e}")
            import traceback
            traceback.print_exc()
            return None
    
    @staticmethod
    def calculate_recommended_params(
        length: int,
        signal_end: int,
        num_peptides: int
    ) -> Dict:
        """Calcule les paramètres recommandés pour une protéine"""
        
        signal_peptide_length = signal_end
        
        estimated_sites = num_peptides * 1.5
        
        if estimated_sites > 12:
            min_cleavage_sites = 5
        elif estimated_sites > 8:
            min_cleavage_sites = 4
        elif estimated_sites > 5:
            min_cleavage_sites = 3
        else:
            min_cleavage_sites = 2
        
        if length < 150:
            min_cleavage_spacing = 3
        elif length < 300:
            min_cleavage_spacing = 4
        else:
            min_cleavage_spacing = 5
        
        max_peptide_length = 100
        
        return {
            "signalPeptideLength": signal_peptide_length,
            "minCleavageSites": min_cleavage_sites,
            "minCleavageSpacing": min_cleavage_spacing,
            "maxPeptideLength": max_peptide_length
        }


# Instance globale
protein_db = ProteinDatabase()