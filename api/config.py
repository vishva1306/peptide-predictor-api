"""Configuration globale de l'API"""
from typing import Dict

class Config:
    """Configuration centralisée"""
    
    # API Info
    API_TITLE = "Peptide Predictor API"
    API_VERSION = "2.2.0"
    API_DESCRIPTION = "API pour prédire les peptides bioactifs (Coassolo et al. Nature 2025)"
    
    # CORS
    CORS_ORIGINS = ["*"]
    
    # Paramètres par défaut (du papier Nature)
    DEFAULT_SIGNAL_PEPTIDE_LENGTH = 20
    DEFAULT_MIN_CLEAVAGE_SITES = 4
    DEFAULT_MIN_CLEAVAGE_SPACING = 5
    
    # Peptides optimaux
    OPTIMAL_PEPTIDE_MIN_LENGTH = 5
    OPTIMAL_PEPTIDE_MAX_LENGTH = 25
    
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
        "permissive": r"(?:KK|KR|RR|RK)"
    }
    
    @classmethod
    def get_regex_pattern(cls, mode: str) -> str:
        """Retourne le pattern regex"""
        return cls.REGEX_PATTERNS.get(mode, cls.REGEX_PATTERNS["strict"])

config = Config()