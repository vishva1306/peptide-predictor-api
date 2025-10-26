"""Extraction des peptides"""
from typing import List, Dict
from api.config import config
from api.models.schemas import CleavageSite
import re

class PeptideExtractor:
    """Extracteur de peptides"""
    
    @staticmethod
    def clean_peptide(peptide: str) -> str:
        """
        Retire les résidus K et R en début et fin de peptide
        Ces résidus font partie des sites de clivage
        """
        # Retirer K/R au début
        peptide = re.sub(r'^[KR]+', '', peptide)
        
        # Retirer K/R à la fin
        peptide = re.sub(r'[KR]+$', '', peptide)
        
        return peptide
    
    @staticmethod
    def extract(
        sequence: str,
        cleavage_sites: List[CleavageSite],
        signal_length: int,
        min_spacing: int,
        min_sites: int
    ) -> List[Dict]:
        """Extrait les peptides entre les sites de clivage"""
        
        if len(cleavage_sites) < min_sites:
            return []
        
        peptides = []
        prev_position = signal_length
        
        for site in cleavage_sites:
            # Utiliser site.index (début du motif de clivage)
            current_pos = site.index
            distance = current_pos - prev_position
            
            if distance >= min_spacing:
                # Extraire le peptide brut
                pep_seq_raw = sequence[prev_position:current_pos]
                
                # ✅ NETTOYER le peptide (retirer K/R terminaux)
                pep_seq_clean = PeptideExtractor.clean_peptide(pep_seq_raw)
                
                # Vérifier qu'il reste quelque chose après nettoyage
                if len(pep_seq_clean) > 3:
                    peptides.append({
                        'sequence': pep_seq_clean,  # ✅ Séquence nettoyée
                        'start': prev_position,
                        'end': current_pos,
                        'length': len(pep_seq_clean),  # ✅ Longueur mise à jour
                        'inRange': config.OPTIMAL_PEPTIDE_MIN_LENGTH <= len(pep_seq_clean) <= config.OPTIMAL_PEPTIDE_MAX_LENGTH,
                        'cleavageMotif': site.motif,
                        'bioactivityScore': 0.0,
                        'bioactivitySource': 'none'
                    })
                
                # Passer APRÈS le motif de clivage (2 aa)
                prev_position = site.position
        
        # Dernier peptide (jusqu'à la fin de la séquence)
        if len(sequence) - prev_position > 3:
            last_seq_raw = sequence[prev_position:]
            
            # ✅ NETTOYER le dernier peptide aussi
            last_seq_clean = PeptideExtractor.clean_peptide(last_seq_raw)
            
            if len(last_seq_clean) > 3:
                peptides.append({
                    'sequence': last_seq_clean,  # ✅ Séquence nettoyée
                    'start': prev_position,
                    'end': len(sequence),
                    'length': len(last_seq_clean),  # ✅ Longueur mise à jour
                    'inRange': config.OPTIMAL_PEPTIDE_MIN_LENGTH <= len(last_seq_clean) <= config.OPTIMAL_PEPTIDE_MAX_LENGTH,
                    'cleavageMotif': 'END',
                    'bioactivityScore': 0.0,
                    'bioactivitySource': 'none'
                })
        
        return peptides