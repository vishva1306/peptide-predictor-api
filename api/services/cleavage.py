"""Détection des sites de clivage PCSK1/3"""
import regex
from typing import List
from api.config import config
from api.models.schemas import CleavageSite

class CleavageDetector:
    """Détecteur de sites de clivage"""
    
    @staticmethod
    def find_sites(
        sequence: str,
        mode: str,
        signal_length: int,
        min_spacing: int
    ) -> List[CleavageSite]:
        """
        Détecte les sites de clivage
        
        STRICT: Applique toutes les contraintes (espacement minimum)
        PERMISSIVE: Détecte TOUS les sites sans vérifier l'espacement
        """
        sites = []
        
        try:
            # Récupérer le pattern
            pattern = config.get_regex_pattern(mode)
            
            # Chercher tous les sites après le peptide signal
            search_region = sequence[signal_length:]
            
            for match in regex.finditer(pattern, search_region):
                # Position absolue dans la séquence originale
                absolute_position = signal_length + match.start()
                
                # ⭐ DIFFÉRENCE ENTRE LES MODES
                if mode == "strict":
                    # Mode STRICT : Vérifier l'espacement minimum entre sites
                    if len(sites) == 0 or (absolute_position - sites[-1].position >= min_spacing):
                        site = CleavageSite(
                            position=absolute_position + 2,  # Position après le motif
                            motif=match.group(),
                            index=absolute_position
                        )
                        sites.append(site)
                else:
                    # Mode PERMISSIVE : Accepter TOUS les sites détectés
                    site = CleavageSite(
                        position=absolute_position + 2,  # Position après le motif
                        motif=match.group(),
                        index=absolute_position
                    )
                    sites.append(site)
        
        except regex.error as e:
            print(f"Erreur regex: {e}")
            return []
        
        return sites
    
    @staticmethod
    def is_prohormone(sites: List[CleavageSite], min_sites: int) -> bool:
        """Vérifie si c'est une prohormone candidate"""
        return len(sites) >= min_sites