"""Calculateur de score amphipathique"""
from typing import Dict, List

class AmphipathicCalculator:
    """Calcule les propriétés amphipathiques des peptides"""
    
    # Acides aminés basiques (chargés +)
    BASIC_AA = ['K', 'R', 'H']
    
    # Acides aminés lipophiliques (hydrophobes)
    LIPOPHILIC_AA = ['A', 'V', 'L', 'I', 'M', 'F', 'W', 'Y']
    
    @staticmethod
    def calculate(sequence: str) -> Dict:
        """
        Calcule le score amphipathique d'un peptide
        
        Score = Coverage % = (basic + lipophilic) / total * 100
        
        Args:
            sequence: Séquence du peptide
            
        Returns:
            Dict avec scores et détails
        """
        if not sequence or len(sequence) == 0:
            return {
                'amphipathicScore': 0,
                'basicCount': 0,
                'lipophilicCount': 0,
                'basicRatio': 0,
                'lipophilicRatio': 0,
                'otherCount': 0,
                'otherRatio': 0,
                'basicDetails': {},
                'lipophilicDetails': {}
            }
        
        total_length = len(sequence)
        
        # Compter les acides aminés basiques
        basic_count = 0
        basic_details = {}
        for aa in AmphipathicCalculator.BASIC_AA:
            count = sequence.count(aa)
            if count > 0:
                basic_details[aa] = count
                basic_count += count
        
        # Compter les acides aminés lipophiliques
        lipophilic_count = 0
        lipophilic_details = {}
        for aa in AmphipathicCalculator.LIPOPHILIC_AA:
            count = sequence.count(aa)
            if count > 0:
                lipophilic_details[aa] = count
                lipophilic_count += count
        
        # Autres acides aminés
        other_count = total_length - basic_count - lipophilic_count
        
        # Calculer ratios
        basic_ratio = (basic_count / total_length * 100) if total_length > 0 else 0
        lipophilic_ratio = (lipophilic_count / total_length * 100) if total_length > 0 else 0
        other_ratio = (other_count / total_length * 100) if total_length > 0 else 0
        
        # Score amphipathique = Coverage (basic + lipophilic)
        amphipathic_score = basic_ratio + lipophilic_ratio
        
        return {
            'amphipathicScore': round(amphipathic_score, 1),
            'basicCount': basic_count,
            'lipophilicCount': lipophilic_count,
            'basicRatio': round(basic_ratio, 1),
            'lipophilicRatio': round(lipophilic_ratio, 1),
            'otherCount': other_count,
            'otherRatio': round(other_ratio, 1),
            'basicDetails': basic_details,
            'lipophilicDetails': lipophilic_details
        }
    
    @staticmethod
    def calculate_batch(sequences: List[str]) -> List[Dict]:
        """
        Calcule le score amphipathique pour plusieurs peptides
        
        Args:
            sequences: Liste de séquences
            
        Returns:
            Liste de dictionnaires avec scores
        """
        return [AmphipathicCalculator.calculate(seq) for seq in sequences]


# Instance globale
amphipathic_calculator = AmphipathicCalculator()