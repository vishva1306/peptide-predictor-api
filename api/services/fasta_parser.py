"""Service pour parser et valider les séquences FASTA"""
import re
from typing import Dict, Optional


class FASTAParser:
    """Parser pour séquences FASTA"""
    
    @staticmethod
    def parse(fasta_text: str) -> Dict[str, Optional[str]]:
        """
        Parse une séquence FASTA avec ou sans header
        
        Args:
            fasta_text: Texte FASTA brut
            
        Returns:
            Dict avec 'sequence', 'header', 'id', 'name'
        """
        lines = fasta_text.strip().split('\n')
        
        header = None
        protein_id = None
        protein_name = None
        sequence_lines = []
        
        for line in lines:
            line = line.strip()
            if not line:
                continue
                
            # Si commence par '>', c'est un header
            if line.startswith('>'):
                header = line[1:].strip()  # Enlever le '>'
                
                # Essayer d'extraire ID et nom du header
                # Format type: >sp|P01308|INS_HUMAN Insulin
                match = re.match(r'^(?:\w+\|)?([A-Z0-9]+)\|?([A-Z0-9_]+)?\s*(.*)?$', header)
                if match:
                    protein_id = match.group(1) or match.group(2)
                    protein_name = match.group(3).strip() if match.group(3) else None
                else:
                    # Header simple sans format standard
                    protein_name = header
            else:
                # C'est une ligne de séquence
                sequence_lines.append(line)
        
        # Joindre toutes les lignes de séquence
        sequence = ''.join(sequence_lines).upper()
        
        return {
            'sequence': sequence,
            'header': header,
            'id': protein_id,
            'name': protein_name
        }
    
    @staticmethod
    def validate_sequence(sequence: str) -> tuple[bool, Optional[str]]:
        """
        Valide qu'une séquence contient uniquement des acides aminés standards
        
        Args:
            sequence: Séquence à valider
            
        Returns:
            (is_valid, error_message)
        """
        if not sequence:
            return False, "Empty sequence"
        
        # Vérifier que la séquence contient uniquement A-Z
        if not re.match(r'^[A-Z]+$', sequence):
            invalid_chars = set(re.findall(r'[^A-Z]', sequence))
            return False, f"Invalid characters found: {', '.join(sorted(invalid_chars))}"
        
        # Vérifier longueur minimale
        if len(sequence) < 30:
            return False, f"Sequence too short ({len(sequence)} aa). Minimum 30 amino acids required."
        
        return True, None


# Instance globale
fasta_parser = FASTAParser()