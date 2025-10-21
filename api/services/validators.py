"""Validation des séquences"""
from fastapi import HTTPException
from api.config import config

class SequenceValidator:
    """Validateur de séquences protéiques"""
    
    @staticmethod
    def clean_sequence(sequence: str) -> str:
        """Nettoie une séquence FASTA/brute"""
        clean = sequence.strip().upper()
        
        # Supprimer header FASTA
        if clean.startswith('>'):
            lines = clean.split('\n')
            clean = ''.join(lines[1:])
        
        # Supprimer espaces
        clean = clean.replace('\n', '').replace(' ', '').replace('\r', '')
        
        return clean
    
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