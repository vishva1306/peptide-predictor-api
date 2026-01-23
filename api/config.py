"""Configuration globale de l'API"""
from typing import Dict

class Config:
    """Configuration centralisée"""
    
    # API Info
    API_TITLE = "Peptide Predictor API"
    API_VERSION = "2.6.0"  # ⭐ Version bump pour PCSK5/6/7
    API_DESCRIPTION = "API pour prédire les peptides bioactifs (Coassolo et al. Nature 2025)"
    
    # CORS
    CORS_ORIGINS = ["*"]
    
    # Paramètres par défaut (du papier Nature) - PCSK1/2
    DEFAULT_SIGNAL_PEPTIDE_LENGTH = 20
    DEFAULT_MIN_CLEAVAGE_SITES = 4
    DEFAULT_MIN_CLEAVAGE_SPACING = 5
    
    # Peptides optimaux PCSK1/2
    OPTIMAL_PEPTIDE_MIN_LENGTH = 5
    OPTIMAL_PEPTIDE_MAX_LENGTH = 25
    
    # ⭐ NOUVEAU : Peptides PCSK5/6/7 (plus grands)
    PCSK567_MIN_LENGTH = 50
    PCSK567_MAX_LENGTH = 500
    
    # Bioactivité
    PEPTIDERANKER_API_URL = "http://peptideranker.ilincs.org/api/predict"
    PEPTIDERANKER_TIMEOUT = 10
    
    # Acides aminés valides
    VALID_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY*")
    
    # Patterns regex
    REGEX_PATTERNS: Dict[str, str] = {
        # MODE STRICT : Regex complète du papier Nature (inchangé)
        "strict": r"(?<!K|R)(?:KK|KR|RR|RK)(?=[^RKILPVH]|$)",
        
        # MODE PERMISSIVE : Simplifié au maximum - juste les motifs de clivage
        "permissive": r"(?:KK|KR|RR|RK)",
        
        # MODE ULTRA-PERMISSIVE (single basic sites)
        "ultra-permissive": r"(?:K|R)",
        
        # ⭐ NOUVEAU : MODE PCSK5/6/7
        # Motif : R-X-(K/R)-R où X = n'importe quel acide aminé
        # Exemples : RSKR, RKRR, RVRR, RARR, RHRR, etc.
        "pcsk567": r"R[A-Z](?:K|R)R"
    }
    
    @classmethod
    def get_regex_pattern(cls, mode: str) -> str:
        """Retourne le pattern regex"""
        return cls.REGEX_PATTERNS.get(mode, cls.REGEX_PATTERNS["strict"])

config = Config()