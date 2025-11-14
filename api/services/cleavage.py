"""Détection des sites de clivage PCSK1/3"""
import regex
import re
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
        ULTRA-PERMISSIVE: Single basic + RF-amide priority scan
        """
        
        if mode == "ultra-permissive":
            return CleavageDetector._find_ultra_permissive_sites(
                sequence, 
                signal_length
            )
        
        # ==================== STRICT / PERMISSIVE (INCHANGÉ) ====================
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
    def _find_ultra_permissive_sites(
        sequence: str,
        signal_length: int
    ) -> List[CleavageSite]:
        """
        Détection ultra-permissive avec RF-amide priority
        
        Stratégie:
        1. Scan RF-amide (haute priorité) - cherche RF/RY partout
        2. Tous les R/K isolés (single basic)
        """
        sites = []
        search_region = sequence[signal_length:]
        
        print(f"\n🔍 Ultra-permissive scan on {len(search_region)} aa")
        
        # ==================== PRIORITÉ 1 : RF-AMIDE SCAN ====================
        # Chercher tous les RF, RFG, RY, RYG dans TOUTE la séquence
        rfamide_patterns = [
            (r'RF(?:G)?', 'RF'),    # RF ou RFG
            (r'RY(?:G)?', 'RY')     # RY ou RYG
        ]
        
        rfamide_sites = []
        for pattern, motif_base in rfamide_patterns:
            for match in re.finditer(pattern, search_region):
                rf_start = match.start()
                rf_end = match.end()
                rf_motif = match.group()
                absolute_rf_start = signal_length + rf_start
                
                print(f"  🟣 Found {rf_motif} at position {absolute_rf_start}")
                
                # Le R du RF est le site de clivage lui-même
                # Chercher le R/K PRÉCÉDENT (pour extraire le peptide)
                found_previous = False
                for lookback in range(1, min(51, rf_start + 1)):
                    check_pos = rf_start - lookback
                    if search_region[check_pos] in 'KR':
                        absolute_pos = signal_length + check_pos
                        
                        rfamide_sites.append({
                            'position': absolute_pos + 1,  # Après le R/K précédent
                            'motif': f"{search_region[check_pos]}...{rf_motif}",
                            'index': absolute_pos,
                            'type': 'rfamide',
                            'rf_position': absolute_rf_start,
                            'rf_end': signal_length + rf_end
                        })
                        found_previous = True
                        print(f"    ✅ RFamide site: {search_region[check_pos]} at {absolute_pos} → {rf_motif} at {absolute_rf_start}")
                        break
                
                if not found_previous:
                    # Pas de R/K avant, le R du RF est le premier site
                    # On crée quand même un site RFamide
                    rfamide_sites.append({
                        'position': absolute_rf_start + len(rf_motif),  # Après le RF
                        'motif': f"START...{rf_motif}",
                        'index': absolute_rf_start,
                        'type': 'rfamide',
                        'rf_position': absolute_rf_start,
                        'rf_end': signal_length + rf_end
                    })
                    print(f"    ⚠️ No previous R/K, using RF itself at {absolute_rf_start}")
        
        print(f"🟣 RF-amide sites found: {len(rfamide_sites)}")
        
        # ==================== PRIORITÉ 2 : TOUS LES R/K ====================
        # Détecter tous les R ou K isolés
        single_basic_count = 0
        for i, aa in enumerate(search_region):
            if aa in 'KR':
                absolute_position = signal_length + i
                
                # Vérifier que ce n'est pas déjà un site RF-amide
                is_rfamide = any(
                    s['index'] == absolute_position 
                    for s in rfamide_sites
                )
                
                if not is_rfamide:
                    sites.append(CleavageSite(
                        position=absolute_position + 1,  # Après le R/K
                        motif=aa,
                        index=absolute_position
                    ))
                    single_basic_count += 1
        
        print(f"🔵 Single basic sites found: {single_basic_count}")
        
        # Convertir RF-amide sites en CleavageSite
        for rf_site in rfamide_sites:
            sites.append(CleavageSite(
                position=rf_site['position'],
                motif=rf_site['motif'],
                index=rf_site['index']
            ))
        
        # Trier par position
        sites.sort(key=lambda s: s.index)
        
        print(f"✅ Total ultra-permissive sites: {len(sites)} ({len(rfamide_sites)} RF-amide + {single_basic_count} single basic)")
        
        return sites
    
    @staticmethod
    def is_prohormone(sites: List[CleavageSite], min_sites: int) -> bool:
        """Vérifie si c'est une prohormone candidate"""
        return len(sites) >= min_sites