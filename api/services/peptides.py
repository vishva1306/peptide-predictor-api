"""Extraction des peptides"""
from typing import List, Dict
from api.config import config
from api.models.schemas import CleavageSite

class PeptideExtractor:
    """Extracteur de peptides"""
    
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
            current_pos = site.position
            distance = current_pos - prev_position
            
            if distance >= min_spacing:
                pep_seq = sequence[prev_position:current_pos]
                
                if len(pep_seq) > 3:
                    peptides.append({
                        'sequence': pep_seq,
                        'start': prev_position,
                        'end': current_pos,
                        'length': len(pep_seq),
                        'inRange': config.OPTIMAL_PEPTIDE_MIN_LENGTH <= len(pep_seq) <= config.OPTIMAL_PEPTIDE_MAX_LENGTH,
                        'cleavageMotif': site.motif,
                        'bioactivityScore': 0.0,
                        'bioactivitySource': 'none'
                    })
                
                prev_position = current_pos
        
        # Dernier peptide
        if len(sequence) - prev_position > 3:
            last_seq = sequence[prev_position:]
            peptides.append({
                'sequence': last_seq,
                'start': prev_position,
                'end': len(sequence),
                'length': len(last_seq),
                'inRange': config.OPTIMAL_PEPTIDE_MIN_LENGTH <= len(last_seq) <= config.OPTIMAL_PEPTIDE_MAX_LENGTH,
                'cleavageMotif': 'END',
                'bioactivityScore': 0.0,
                'bioactivitySource': 'none'
            })
        
        return peptides