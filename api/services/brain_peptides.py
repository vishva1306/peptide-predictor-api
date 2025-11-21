"""
Brain Peptides Detection Service
Compare les peptides prédits avec le dataset de peptides détectés dans le cerveau humain
(Nature Communications 2016, Zougman et al.)
"""

import json
from pathlib import Path
from typing import Dict, Optional

class BrainPeptidesChecker:
    """Service pour vérifier si un peptide a été détecté dans le cerveau"""
    
    def __init__(self):
        self.brain_peptides: Dict[str, dict] = {}
        self.loaded = False
        self.total_count = 0
        self.metadata = {}
        self._load_peptides()
    
    def _load_peptides(self):
        """Charge les peptides du cerveau en mémoire (une seule fois au démarrage)"""
        try:
            data_file = Path("data/brain_peptides.json")
            
            if not data_file.exists():
                print("⚠️ Brain peptides dataset not found. Feature disabled.")
                print("   Expected location: data/brain_peptides.json")
                return
            
            print(f"📖 Loading brain peptides from {data_file}...")
            
            with open(data_file, 'r') as f:
                data = json.load(f)
            
            self.brain_peptides = data['peptides']
            self.total_count = data['total_peptides']
            self.metadata = {
                'source': data['source'],
                'doi': data['doi'],
                'reference': data['reference'],
                'statistics': data.get('statistics', {})
            }
            self.loaded = True
            
            # Calculer taille mémoire approximative
            memory_mb = len(str(self.brain_peptides)) / 1024 / 1024
            
            print(f"✅ Loaded {self.total_count:,} brain peptides into memory (~{memory_mb:.2f} MB)")
            print(f"   Source: {self.metadata['reference']}")
            print(f"   Pro-hormone peptides: {self.metadata['statistics'].get('prohormone_peptides', 0):,}")
            print(f"   Amidated peptides: {self.metadata['statistics'].get('amidated_peptides', 0):,}")
            
        except Exception as e:
            print(f"❌ Error loading brain peptides: {e}")
            import traceback
            traceback.print_exc()
            self.loaded = False
    
    def check(self, sequence: str) -> Optional[dict]:
        """
        Vérifie si une séquence peptidique a été détectée dans le cerveau
        Gère aussi les peptides amidés (avec G C-terminal retiré)
        
        Args:
            sequence: Séquence peptidique à vérifier
            
        Returns:
            Dict avec les données brain si trouvé, None sinon
            {
                'found': True,
                'isProhormone': True,
                'proteinName': 'Neuropeptide FF',
                'uniprot': 'O00471',
                'msmsCount': 156,
                'mascotScore': 78.5,
                'isAmidated': True,
                'matchNote': 'Optional note about matching' (si match avec amidation)
            }
        """
        if not self.loaded:
            return None
        
        # Normaliser la séquence (uppercase, strip)
        seq_clean = sequence.strip().upper()
        
        # ✅ ÉTAPE 1 : Vérification O(1) - séquence normale (exact match)
        if seq_clean in self.brain_peptides:
            brain_data = self.brain_peptides[seq_clean]
            return {
                'found': True,
                'isProhormone': brain_data['isProhormone'],
                'proteinName': brain_data['proteinName'],
                'uniprot': brain_data['uniprot'],
                'msmsCount': brain_data['msmsCount'],
                'mascotScore': brain_data['mascotScore'],
                'isAmidated': brain_data['isAmidated']
            }
        
        # ⭐ ÉTAPE 2 : Essayer sans le dernier G (cas amidation)
        # Les peptides amidés perdent leur G C-terminal lors de la maturation :
        # YGGFMRFG (précurseur) → YGGFMRF-NH₂ (amidé)
        # Notre extracteur trouve "YGGFMRFG" mais le brain dataset contient "YGGFMRF"
        if len(seq_clean) > 3 and seq_clean[-1] == 'G':
            seq_without_g = seq_clean[:-1]
            
            if seq_without_g in self.brain_peptides:
                brain_data = self.brain_peptides[seq_without_g]
                
                # ⚠️ Ne retourner que si le peptide brain est vraiment amidé
                # (pour éviter les faux positifs)
                if brain_data['isAmidated']:
                    return {
                        'found': True,
                        'isProhormone': brain_data['isProhormone'],
                        'proteinName': brain_data['proteinName'],
                        'uniprot': brain_data['uniprot'],
                        'msmsCount': brain_data['msmsCount'],
                        'mascotScore': brain_data['mascotScore'],
                        'isAmidated': brain_data['isAmidated'],
                        'matchNote': 'Matched after C-terminal amidation (G removed)'
                    }
        
        # ❌ Aucun match trouvé
        return None
    
    def check_batch(self, sequences: list[str]) -> list[Optional[dict]]:
        """
        Vérifie un batch de séquences
        
        Args:
            sequences: Liste de séquences à vérifier
            
        Returns:
            Liste de dicts (données brain si trouvé, None sinon)
        """
        if not self.loaded:
            return [None] * len(sequences)
        
        return [self.check(seq) for seq in sequences]
    
    def get_stats(self) -> dict:
        """Retourne les statistiques du dataset"""
        return {
            "loaded": self.loaded,
            "total_peptides": self.total_count,
            "source": self.metadata.get('source', 'Unknown'),
            "doi": self.metadata.get('doi', 'Unknown'),
            "reference": self.metadata.get('reference', 'Unknown'),
            "statistics": self.metadata.get('statistics', {})
        }


# Instance globale (singleton)
brain_checker = BrainPeptidesChecker()