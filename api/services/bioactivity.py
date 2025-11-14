"""Prédiction de bioactivité"""
import asyncio
import aiohttp
from typing import Tuple, Optional, List
from api.config import config

class BioactivityPredictor:
    """Prédiction de bioactivité"""
    
    @staticmethod
    async def predict_peptideranker(
        peptide: str,
        session: aiohttp.ClientSession
    ) -> Optional[float]:
        """Appelle l'API PeptideRanker"""
        if len(peptide) < 2:
            return None
        
        try:
            payload = {"sequence": peptide}
            timeout = aiohttp.ClientTimeout(total=config.PEPTIDERANKER_TIMEOUT)
            
            async with session.post(
                config.PEPTIDERANKER_API_URL,
                json=payload,
                timeout=timeout
            ) as response:
                if response.status == 200:
                    data = await response.json()
                    score = data.get('score', 0)
                    return float(score * 100)
        
        except asyncio.TimeoutError:
            print(f"Timeout PeptideRanker")
        except Exception as e:
            print(f"Erreur PeptideRanker: {e}")
        
        return None
    
    @staticmethod
    def calculate_heuristic(
        peptide: str, 
        cleavage_motif: str = None,
        full_protein_sequence: str = None,
        peptide_end_position: int = None
    ) -> float:
        """
        Calcul heuristique de bioactivité
        
        ⭐ NOUVEAU :
        - Bonus RFamide (+25 points)
        - Pénalité pour fragments C-terminaux sans glycine (-20 points)
        - Bonus pour peptides connus (Secretoneurin, etc.)
        
        Args:
            peptide: Séquence du peptide
            cleavage_motif: Motif de clivage détecté
            full_protein_sequence: Séquence complète de la protéine (optionnel)
            peptide_end_position: Position de fin du peptide (1-indexed, optionnel)
        """
        if len(peptide) == 0:
            return 0.0
        
        score = 0.0
        
        # ==================== BASE SCORE ====================
        
        # 1. Hydrophobicité (30%)
        hydrophobic_aa = set('ALIVMFWP')
        hydro_count = sum(1 for aa in peptide if aa in hydrophobic_aa)
        hydro_ratio = hydro_count / len(peptide)
        score += hydro_ratio * 30
        
        # 2. Charge (20%)
        positive = sum(1 for aa in peptide if aa in 'KRH')
        negative = sum(1 for aa in peptide if aa in 'DE')
        
        if positive > 0:
            score += 10
        if negative > 0:
            score += 10
        
        # 3. Longueur optimale (35%)
        length = len(peptide)
        if config.OPTIMAL_PEPTIDE_MIN_LENGTH <= length <= config.OPTIMAL_PEPTIDE_MAX_LENGTH:
            score += 35
        elif length < config.OPTIMAL_PEPTIDE_MIN_LENGTH:
            score -= 10
        elif length > 100:
            score -= 15
        
        # 4. Stabilité (15%)
        if 'C' in peptide:
            score += 8
        
        proline_count = peptide.count('P')
        if proline_count <= 2:
            score += 7
        else:
            score -= 5
        
        unique_aa = len(set(peptide))
        if unique_aa >= 6:
            score += 5
        
        # ==================== BONUS & PÉNALITÉS ====================
        
        # ⭐ BONUS 1 : RFamide peptides (+25 points)
        if peptide.endswith('RF') or peptide.endswith('RFG'):
            score += 25
            print(f"🎯 RFamide peptide (RF terminal): +25 bonus → Score before cap: {score:.1f}")
        elif peptide.endswith('RY') or peptide.endswith('RYG'):
            score += 25
            print(f"🎯 RFamide peptide (RY terminal): +25 bonus → Score before cap: {score:.1f}")
        
        # Bonus supplémentaire pour motif RFamide dans le cleavage
        if cleavage_motif and ('RF' in cleavage_motif or 'RY' in cleavage_motif):
            score += 10
            print(f"🎯 RFamide cleavage motif: +10 bonus → Score before cap: {score:.1f}")
        
        # ⭐ BONUS 2 : Peptides bien établis (Secretoneurin, etc.)
        # Patterns de peptides bioactifs connus
        known_bioactive_patterns = [
            ('SECRETONEURIN', ['SNSQE', 'PGKQL', 'RLERL']),  # Secretoneurin motifs
            ('CHROMOGRANIN', ['WPRES', 'LQEEE', 'HLEAE']),   # Chromogranin motifs
            ('VGF', ['TLQP', 'AQEE', 'NERP']),               # VGF peptide motifs
        ]
        
        for peptide_name, motifs in known_bioactive_patterns:
            for motif in motifs:
                if motif in peptide:
                    score += 15
                    print(f"🧬 Known bioactive motif ({peptide_name}): +15 bonus → Score: {score:.1f}")
                    break
        
        # ⭐ PÉNALITÉ 1 : Fragments C-terminaux sans glycine (-20 points)
        # Un peptide C-terminal DOIT avoir une glycine suivie de R/K pour être amidé
        if full_protein_sequence and peptide_end_position:
            is_c_terminal_fragment = (peptide_end_position >= len(full_protein_sequence) - 5)
            has_terminal_glycine = peptide.endswith('G')
            
            # Vérifier s'il y a R/K après le peptide (pour amidation)
            has_basic_after = False
            if peptide_end_position < len(full_protein_sequence):
                next_aa = full_protein_sequence[peptide_end_position:peptide_end_position + 2]
                has_basic_after = any(aa in 'KR' for aa in next_aa)
            
            # Pénalité si c'est un fragment C-terminal SANS glycine terminale et SANS R/K après
            if is_c_terminal_fragment and not has_terminal_glycine and not has_basic_after:
                score -= 20
                print(f"⚠️ C-terminal fragment without glycine/basic: -20 penalty → Score: {score:.1f}")
        
        # ⭐ PÉNALITÉ 2 : Peptides trop courts sans caractéristiques spéciales
        if length < 5 and not (peptide.endswith('RF') or peptide.endswith('RY')):
            score -= 15
            print(f"⚠️ Very short peptide without RFamide: -15 penalty → Score: {score:.1f}")
        
        # ⭐ PÉNALITÉ 3 : Peptides très basiques (trop de K/R)
        basic_ratio = (positive / length) if length > 0 else 0
        if basic_ratio > 0.5:  # Plus de 50% K/R
            score -= 10
            print(f"⚠️ Too many basic residues ({basic_ratio:.0%}): -10 penalty → Score: {score:.1f}")
        
        # ==================== FINAL ADJUSTMENTS ====================
        
        # Cap final entre 0 et 100
        final_score = max(min(score, 100), 0)
        
        return final_score
    
    @classmethod
    async def predict(
        cls,
        peptide: str,
        session: aiohttp.ClientSession,
        cleavage_motif: str = None,
        full_protein_sequence: str = None,
        peptide_end_position: int = None
    ) -> Tuple[float, str]:
        """
        Prédit la bioactivité (avec contexte protéine si disponible)
        
        Args:
            peptide: Séquence du peptide
            session: Session aiohttp
            cleavage_motif: Motif de clivage
            full_protein_sequence: Séquence complète de la protéine
            peptide_end_position: Position de fin du peptide (1-indexed)
        """
        # Essayer l'API PeptideRanker d'abord
        api_score = await cls.predict_peptideranker(peptide, session)
        
        if api_score is not None:
            return api_score, "api"
        
        # Fallback sur heuristique avec contexte
        heuristic_score = cls.calculate_heuristic(
            peptide, 
            cleavage_motif,
            full_protein_sequence,
            peptide_end_position
        )
        return heuristic_score, "heuristic"
    
    @classmethod
    async def predict_batch(
        cls,
        peptides: List[str],
        session: aiohttp.ClientSession,
        cleavage_motifs: List[str] = None,
        full_protein_sequence: str = None,
        peptide_end_positions: List[int] = None
    ) -> List[Tuple[float, str]]:
        """
        Prédit en batch avec contexte protéine
        
        Args:
            peptides: Liste de séquences
            session: Session aiohttp
            cleavage_motifs: Liste des motifs de clivage
            full_protein_sequence: Séquence complète de la protéine
            peptide_end_positions: Liste des positions de fin (1-indexed)
        """
        # Préparer les arguments
        if cleavage_motifs is None:
            cleavage_motifs = [None] * len(peptides)
        
        if peptide_end_positions is None:
            peptide_end_positions = [None] * len(peptides)
        
        # Créer les tâches
        tasks = [
            cls.predict(
                peptide, 
                session, 
                motif,
                full_protein_sequence,
                end_pos
            ) 
            for peptide, motif, end_pos in zip(peptides, cleavage_motifs, peptide_end_positions)
        ]
        
        return await asyncio.gather(*tasks)