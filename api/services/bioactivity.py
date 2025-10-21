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
    def calculate_heuristic(peptide: str) -> float:
        """Calcul heuristique de bioactivité"""
        if len(peptide) == 0:
            return 0.0
        
        score = 0.0
        
        # Hydrophobicité (30%)
        hydrophobic_aa = set('ALIVMFWP')
        hydro_count = sum(1 for aa in peptide if aa in hydrophobic_aa)
        hydro_ratio = hydro_count / len(peptide)
        score += hydro_ratio * 30
        
        # Charge (20%)
        positive = sum(1 for aa in peptide if aa in 'KRH')
        negative = sum(1 for aa in peptide if aa in 'DE')
        
        if positive > 0:
            score += 10
        if negative > 0:
            score += 10
        
        # Longueur optimale (35%)
        length = len(peptide)
        if config.OPTIMAL_PEPTIDE_MIN_LENGTH <= length <= config.OPTIMAL_PEPTIDE_MAX_LENGTH:
            score += 35
        elif length < config.OPTIMAL_PEPTIDE_MIN_LENGTH:
            score -= 10
        elif length > 100:
            score -= 15
        
        # Stabilité (15%)
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
        
        return max(min(score, 100), 0)
    
    @classmethod
    async def predict(
        cls,
        peptide: str,
        session: aiohttp.ClientSession
    ) -> Tuple[float, str]:
        """Prédit la bioactivité"""
        api_score = await cls.predict_peptideranker(peptide, session)
        
        if api_score is not None:
            return api_score, "api"
        
        heuristic_score = cls.calculate_heuristic(peptide)
        return heuristic_score, "heuristic"
    
    @classmethod
    async def predict_batch(
        cls,
        peptides: List[str],
        session: aiohttp.ClientSession
    ) -> List[Tuple[float, str]]:
        """Prédit en batch"""
        tasks = [cls.predict(peptide, session) for peptide in peptides]
        return await asyncio.gather(*tasks)