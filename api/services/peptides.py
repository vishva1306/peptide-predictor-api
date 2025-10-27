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
        min_sites: int,
        mode: str = "strict"
    ) -> List[Dict]:
        """
        Extrait les peptides entre les sites de clivage
        
        STRICT: Respecte l'espacement minimum entre peptides
        PERMISSIVE: Extrait TOUS les peptides même très courts
        
        ⭐ IMPORTANT: Les sites de clivage (KK/KR/RR/RK) sont EXCLUS des peptides
        """
        
        if len(cleavage_sites) < min_sites:
            return []
        
        peptides = []
        prev_position = signal_length
        
        for site in cleavage_sites:
            # ⭐ CORRECTION: site.position est APRÈS le motif
            # On veut le peptide SANS le motif de clivage
            # Donc on prend jusqu'à site.index (début du motif)
            current_pos = site.index  # Position du DÉBUT du motif (pas après)
            distance = current_pos - prev_position
            
            # ⭐ DIFFÉRENCE ENTRE LES MODES
            if mode == "strict":
                # Mode STRICT : Vérifier l'espacement minimum
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
                    
                    # ⭐ Avancer APRÈS le motif de clivage (2 aa)
                    prev_position = site.position
            else:
                # Mode PERMISSIVE : Extraire TOUS les peptides sans vérifier l'espacement
                pep_seq = sequence[prev_position:current_pos]
                
                # Accepter même les peptides de 1 aa en mode PERMISSIVE
                if len(pep_seq) > 0:
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
                
                # ⭐ Avancer APRÈS le motif de clivage (2 aa)
                prev_position = site.position
        
        # Dernier peptide (après le dernier site de clivage)
        if len(sequence) - prev_position > 0:
            last_seq = sequence[prev_position:]
            
            # En mode STRICT : minimum 3 aa
            # En mode PERMISSIVE : accepter tout
            min_length = 3 if mode == "strict" else 0
            
            if len(last_seq) > min_length:
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