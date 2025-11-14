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
        ULTRA-PERMISSIVE: Extrait TOUS + calcul confidence score + INCLUT RF TERMINAL
        """
        
        if mode == "ultra-permissive":
            return PeptideExtractor._extract_ultra_permissive(
                sequence,
                cleavage_sites,
                signal_length
            )
        
        # ==================== STRICT / PERMISSIVE (INCHANGÉ) ====================
        if len(cleavage_sites) < min_sites:
            return []
        
        peptides = []
        prev_position = signal_length
        
        for site in cleavage_sites:
            current_pos = site.index
            distance = current_pos - prev_position
            
            if mode == "strict":
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
                    
                    prev_position = site.position
            else:
                pep_seq = sequence[prev_position:current_pos]
                
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
                
                prev_position = site.position
        
        # Dernier peptide
        if len(sequence) - prev_position > 0:
            last_seq = sequence[prev_position:]
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
    
    @staticmethod
    def _extract_ultra_permissive(
        sequence: str,
        cleavage_sites: List[CleavageSite],
        signal_length: int
    ) -> List[Dict]:
        """
        Extraction ultra-permissive OPTIMISÉE
        
        ⭐ NOUVEAU : 
        - INCLUT le RF/RY terminal dans le peptide
        - MIN_CONFIDENCE = 30 (au lieu de 15)
        - MAX_PEPTIDES = 50 (au lieu de 200)
        - Filtre anti-doublons (overlap > 70%)
        - Priorité aux peptides courts (5-30 aa)
        """
        peptides = []
        
        if len(cleavage_sites) == 0:
            return []
        
        print(f"\n🧬 Ultra-permissive extraction from {len(cleavage_sites)} sites...")
        
        sorted_sites = sorted(cleavage_sites, key=lambda s: s.index)
        
        # ⭐ FILTRE 1 : Limites
        MAX_DISTANCE = 100  # Maximum 100 aa entre deux sites
        MAX_LENGTH = 50     # Peptides > 50 aa sont trop longs
        MIN_CONFIDENCE = 30 # Score minimum (au lieu de 15)
        
        # Extraire peptides entre chaque paire de sites
        for i in range(len(sorted_sites)):
            for j in range(i + 1, len(sorted_sites)):
                site_start = sorted_sites[i]
                site_end = sorted_sites[j]
                
                start_pos = site_start.position
                
                # ⭐ CRITIQUE : Pour les motifs RFamide, INCLURE le RF terminal
                if 'RF' in site_end.motif or 'RY' in site_end.motif:
                    # Le peptide va jusqu'à la FIN du motif RF/RY
                    end_pos = site_end.position  # Après le RF
                    print(f"   🎯 RFamide peptide: including terminal {site_end.motif}")
                else:
                    # Peptide standard : s'arrête au DÉBUT du motif
                    end_pos = site_end.index
                
                peptide_length = end_pos - start_pos
                
                # ⭐ FILTRE : Distance et longueur
                if peptide_length > MAX_DISTANCE or peptide_length > MAX_LENGTH:
                    continue
                
                if peptide_length < 3:
                    continue
                
                peptide_seq = sequence[start_pos:end_pos]
                
                # Calculer score de confiance
                confidence = PeptideExtractor._calculate_confidence(
                    peptide_seq,
                    site_start,
                    site_end,
                    peptide_length
                )
                
                # ⭐ FILTRE 2 : Confidence minimum = 30
                if confidence < MIN_CONFIDENCE:
                    continue
                
                peptides.append({
                    'sequence': peptide_seq,
                    'start': start_pos,
                    'end': end_pos,
                    'length': peptide_length,
                    'inRange': config.OPTIMAL_PEPTIDE_MIN_LENGTH <= peptide_length <= config.OPTIMAL_PEPTIDE_MAX_LENGTH,
                    'cleavageMotif': PeptideExtractor._get_cleavage_label(site_start, site_end),
                    'bioactivityScore': 0.0,
                    'bioactivitySource': 'none',
                    'confidenceScore': confidence,
                    'confidenceBadge': PeptideExtractor._get_confidence_badge(confidence)
                })
        
        # ⭐ FILTRE 3 : Anti-doublons (éliminer peptides qui se chevauchent > 70%)
        peptides = PeptideExtractor._remove_overlapping_peptides(peptides)
        
        # ⭐ FILTRE 4 : Trier et limiter à top 50
        peptides.sort(key=lambda p: (
            p['confidenceScore'],  # Priorité 1 : Score
            -p['length']          # Priorité 2 : Plus court = mieux
        ), reverse=True)
        
        MAX_PEPTIDES = 50
        if len(peptides) > MAX_PEPTIDES:
            print(f"⚠️ Truncating from {len(peptides)} to top {MAX_PEPTIDES} peptides")
            peptides = peptides[:MAX_PEPTIDES]
        
        print(f"✅ Ultra-permissive extracted: {len(peptides)} peptides")
        if len(peptides) > 0:
            print(f"   Top confidence: {peptides[0]['confidenceScore']}")
            print(f"   Top peptide: {peptides[0]['sequence'][:30]}...")
        
        return peptides
    
    @staticmethod
    def _remove_overlapping_peptides(peptides: List[Dict]) -> List[Dict]:
        """
        Élimine les peptides qui se chevauchent à plus de 70%
        Garde celui avec le meilleur score de confiance
        """
        if len(peptides) <= 1:
            return peptides
        
        # Trier par score (meilleurs d'abord)
        sorted_peptides = sorted(peptides, key=lambda p: p['confidenceScore'], reverse=True)
        
        filtered = []
        
        for peptide in sorted_peptides:
            # Vérifier si ce peptide chevauche un peptide déjà gardé
            is_overlapping = False
            
            for kept_peptide in filtered:
                overlap = PeptideExtractor._calculate_overlap(peptide, kept_peptide)
                
                # Si chevauchement > 70%, on ignore ce peptide
                if overlap > 0.7:
                    is_overlapping = True
                    print(f"   🚫 Removing overlapping: {peptide['sequence'][:20]}... (overlap {overlap:.0%})")
                    break
            
            if not is_overlapping:
                filtered.append(peptide)
        
        print(f"   🔄 Removed {len(sorted_peptides) - len(filtered)} overlapping peptides")
        return filtered
    
    @staticmethod
    def _calculate_overlap(pep1: Dict, pep2: Dict) -> float:
        """
        Calcule le pourcentage de chevauchement entre deux peptides
        Retourne 0.0 (aucun) à 1.0 (identique)
        """
        start1, end1 = pep1['start'], pep1['end']
        start2, end2 = pep2['start'], pep2['end']
        
        # Pas de chevauchement
        if end1 <= start2 or end2 <= start1:
            return 0.0
        
        # Calculer la zone de chevauchement
        overlap_start = max(start1, start2)
        overlap_end = min(end1, end2)
        overlap_length = overlap_end - overlap_start
        
        # Chevauchement par rapport au plus petit peptide
        min_length = min(end1 - start1, end2 - start2)
        
        return overlap_length / min_length if min_length > 0 else 0.0
    
    @staticmethod
    def _calculate_confidence(
        peptide_seq: str,
        site_start: CleavageSite,
        site_end: CleavageSite,
        length: int
    ) -> int:
        """
        Calcule le score de confiance (0-100)
        
        ⭐ MODIFIÉ : Bonus pour peptides courts et RFamide
        """
        score = 0
        
        # 1. Type de clivage N-terminal (15-50 pts)
        start_motif = site_start.motif
        if '...' in start_motif:  # RFamide
            score += 50
        elif start_motif in ['KK', 'KR', 'RR', 'RK']:  # Dibasic classique
            score += 50
        else:  # Single basic
            score += 15
        
        # 2. Type de clivage C-terminal (0-50 pts)
        end_motif = site_end.motif
        if '...' in end_motif:  # RFamide C-term
            score += 50
        
        # 3. Motif terminal (0-30 pts)
        terminal_motif = PeptideExtractor._get_terminal_motif(peptide_seq)
        if terminal_motif in ['RF', 'RY', 'RFG', 'RYG']:
            score += 30
            print(f"      🎯 Terminal RFamide: +30 pts")
        elif terminal_motif == 'G':
            score += 15
        
        # 4. Longueur optimale (0-20 pts) ⭐ BONUS POUR COURTS
        if 5 <= length <= 15:
            score += 20  # Peptides courts = meilleurs
        elif 15 < length <= 30:
            score += 10
        elif 30 < length <= 50:
            score += 5
        elif length > 100:
            score -= 30  # Pénalité pour trop longs
        
        # Bonus RFamide (garantir score minimum 90)
        if '...' in start_motif or '...' in end_motif:
            score = max(score, 90)
            print(f"      🎯 RFamide bonus: minimum score 90")
        
        # Plafonner à 100
        final_score = min(max(score, 0), 100)
        
        return final_score
    
    @staticmethod
    def _get_terminal_motif(sequence: str) -> str:
        """Détecte motif terminal spécial"""
        if len(sequence) < 2:
            return 'none'
        
        # Vérifier les 3 derniers aa pour RFG/RYG
        if len(sequence) >= 3:
            last_three = sequence[-3:]
            if last_three in ['RFG', 'RYG']:
                return last_three
        
        # Vérifier les 2 derniers aa pour RF/RY
        last_two = sequence[-2:]
        if last_two in ['RF', 'RY']:
            return last_two
        elif sequence[-1] == 'G':
            return 'G'
        
        return 'none'
    
    @staticmethod
    def _get_cleavage_label(site_start: CleavageSite, site_end: CleavageSite) -> str:
        """Génère label pour le motif de clivage"""
        if '...' in site_start.motif:
            return site_start.motif
        if '...' in site_end.motif:
            return site_end.motif
        
        if site_start.motif in ['KK', 'KR', 'RR', 'RK']:
            return site_start.motif
        
        return f"{site_start.motif}→{site_end.motif}"
    
    @staticmethod
    def _get_confidence_badge(score: int) -> str:
        """Retourne le badge de confiance"""
        if score >= 80:
            return "high"
        elif score >= 60:
            return "medium"
        elif score >= 40:
            return "low"
        else:
            return "very_low"