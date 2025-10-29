"""Validation des séquences"""
from fastapi import HTTPException
from api.config import config

class SequenceValidator:
    """Validateur de séquences protéiques"""
    
    @staticmethod
    def clean_sequence(sequence: str) -> tuple:
        """
        Nettoie une séquence FASTA/brute
        
        Returns:
            tuple: (cleaned_sequence, protein_id)
        """
        clean = sequence.strip().upper()
        protein_id = ""  # Par défaut vide
        
        # Supprimer header FASTA et extraire l'ID
        if clean.startswith('>'):
            lines = clean.split('\n')
            # Extraire l'ID (tout ce qui suit '>' sur la première ligne)
            protein_id = lines[0][1:].strip()  # Enlève le '>' et espaces
            clean = ''.join(lines[1:])
        
        # Supprimer espaces
        clean = clean.replace('\n', '').replace(' ', '').replace('\r', '')
        
        return clean, protein_id
    
    @staticmethod
    def validate_characters(sequence: str) -> None:
        """Vérifie les caractères valides"""
        invalid = set(sequence) - config.VALID_AMINO_ACIDS
        if invalid:
            raise HTTPException(
                status_code=400,
                detail=f"Caractères invalides: {', '.join(sorted(invalid))}"
            )
    
    @staticmethod
    def validate_length(sequence: str, min_length: int) -> None:
        """Vérifie la longueur minimale"""
        if len(sequence) < min_length:
            raise HTTPException(
                status_code=400,
                detail=f"Séquence trop courte. Min: {min_length} aa (actuel: {len(sequence)} aa)"
            )