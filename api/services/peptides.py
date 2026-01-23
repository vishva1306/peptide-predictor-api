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
        
        ⭐ NOUVEAU : Capture les motifs N-terminal et C-terminal pour chaque peptide
        
        STRICT: Respecte l'espacement minimum entre peptides
        PERMISSIVE: Extrait TOUS les peptides même très courts
        ULTRA-PERMISSIVE: Extrait TOUS + calcul confidence score + INCLUT RF TERMINAL
        PCSK567: Extrait forme mature après clivage R-X-K/R-R
        """
        
        if mode == "ultra-permissive":
            return PeptideExtractor._extract_ultra_permissive(
                sequence,
                cleavage_sites,
                signal_length
            )
        
        if mode == "pcsk567":
            return PeptideExtractor._extract_pcsk567(
                sequence,
                cleavage_sites,
                signal_length
            )
        
        # ==================== STRICT / PERMISSIVE ====================
        if len(cleavage_sites) < min_sites:
            return []
        
        peptides = []
        prev_position = signal_length
        prev_motif = "SIGNAL"  # ⭐ NOUVEAU : Motif N-terminal du premier peptide
        
        for site in cleavage_sites:
            current_pos = site.index
            distance = current_pos - prev_position
            
            if mode == "strict":
                if distance >= min_spacing:
                    pep_seq = sequence[prev_position:current_pos]
                    
                    if len(pep_seq) > 3:
                        peptides.append({
                            'sequence': pep_seq,
                            'start': prev_position + 1,  # 1-indexed
                            'end': current_pos,
                            'length': len(pep_seq),
                            'inRange': config.OPTIMAL_PEPTIDE_MIN_LENGTH <= len(pep_seq) <= config.OPTIMAL_PEPTIDE_MAX_LENGTH,
                            'cleavageMotifN': prev_motif,  # ⭐ NOUVEAU
                            'cleavageMotifC': site.motif,  # ⭐ NOUVEAU
                            'cleavageMotif': site.motif,   # Compatibilité
                            'bioactivityScore': 0.0,
                            'bioactivitySource': 'none'
                        })
                    
                    prev_position = site.position
                    prev_motif = site.motif  # ⭐ NOUVEAU : Mémoriser le motif pour le prochain peptide
            else:
                pep_seq = sequence[prev_position:current_pos]
                
                if len(pep_seq) > 0:
                    peptides.append({
                        'sequence': pep_seq,
                        'start': prev_position + 1,  # 1-indexed
                        'end': current_pos,
                        'length': len(pep_seq),
                        'inRange': config.OPTIMAL_PEPTIDE_MIN_LENGTH <= len(pep_seq) <= config.OPTIMAL_PEPTIDE_MAX_LENGTH,
                        'cleavageMotifN': prev_motif,  # ⭐ NOUVEAU
                        'cleavageMotifC': site.motif,  # ⭐ NOUVEAU
                        'cleavageMotif': site.motif,   # Compatibilité
                        'bioactivityScore': 0.0,
                        'bioactivitySource': 'none'
                    })
                
                prev_position = site.position
                prev_motif = site.motif  # ⭐ NOUVEAU
        
        # Dernier peptide
        if len(sequence) - prev_position > 0:
            last_seq = sequence[prev_position:]
            min_length = 3 if mode == "strict" else 0
            
            if len(last_seq) > min_length:
                peptides.append({
                    'sequence': last_seq,
                    'start': prev_position + 1,  # 1-indexed
                    'end': len(sequence),
                    'length': len(last_seq),
                    'inRange': config.OPTIMAL_PEPTIDE_MIN_LENGTH <= len(last_seq) <= config.OPTIMAL_PEPTIDE_MAX_LENGTH,
                    'cleavageMotifN': prev_motif,  # ⭐ NOUVEAU
                    'cleavageMotifC': 'END',       # ⭐ NOUVEAU
                    'cleavageMotif': 'END',        # Compatibilité
                    'bioactivityScore': 0.0,
                    'bioactivitySource': 'none'
                })
        
        return peptides
    
    # ⭐ NOUVEAU : Extraction pour PCSK5/6/7 avec motifs N et C
    @staticmethod
    def _extract_pcsk567(
        sequence: str,
        cleavage_sites: List[CleavageSite],
        signal_length: int
    ) -> List[Dict]:
        """
        Extraction spécifique PCSK5/6/7
        
        Différences avec PCSK1/2 :
        - Un seul site de clivage suffit
        - On extrait la partie APRÈS le clivage (forme mature)
        - Peptides généralement plus grands (50-200+ aa)
        - Souvent en C-terminal
        
        ⭐ NOUVEAU : Capture motifs N et C terminal
        """
        peptides = []
        
        print(f"\n🧬 PCSK5/6/7 extraction from {len(cleavage_sites)} site(s)...")
        
        if len(cleavage_sites) == 0:
            print("   ❌ No cleavage sites found")
            return []
        
        for i, site in enumerate(cleavage_sites):
            cleavage_pos = site.position  # Position après le motif R-X-K/R-R
            
            # ==================== FORME MATURE (après clivage) ====================
            mature_seq = sequence[cleavage_pos:]
            
            if len(mature_seq) >= 10:
                peptides.append({
                    'sequence': mature_seq,
                    'start': cleavage_pos + 1,  # 1-indexed
                    'end': len(sequence),
                    'length': len(mature_seq),
                    'inRange': config.PCSK567_MIN_LENGTH <= len(mature_seq) <= config.PCSK567_MAX_LENGTH,
                    'cleavageMotifN': site.motif,  # ⭐ NOUVEAU : Le motif PCSK5/6/7 est en N-term de la forme mature
                    'cleavageMotifC': 'END',       # ⭐ NOUVEAU : C'est la fin de la protéine
                    'cleavageMotif': site.motif,   # Compatibilité
                    'bioactivityScore': 0.0,
                    'bioactivitySource': 'none',
                    'peptideType': 'mature_form',
                    'cleavedBy': 'PCSK5/6/7'
                })
                
                print(f"   ✅ Mature form: {len(mature_seq)} aa")
                print(f"      N-term motif: {site.motif}")
                print(f"      C-term motif: END")
                print(f"      Sequence: {mature_seq[:50]}...")
            
            # ==================== PRODOMAIN (avant clivage) ====================
            prodomain_start = signal_length
            prodomain_seq = sequence[prodomain_start:site.index]
            
            if len(prodomain_seq) >= 20:
                peptides.append({
                    'sequence': prodomain_seq,
                    'start': prodomain_start + 1,  # 1-indexed
                    'end': site.index,
                    'length': len(prodomain_seq),
                    'inRange': config.PCSK567_MIN_LENGTH <= len(prodomain_seq) <= config.PCSK567_MAX_LENGTH,
                    'cleavageMotifN': 'SIGNAL',    # ⭐ NOUVEAU : Commence après le signal peptide
                    'cleavageMotifC': site.motif,  # ⭐ NOUVEAU : Se termine au site PCSK5/6/7
                    'cleavageMotif': site.motif,   # Compatibilité
                    'bioactivityScore': 0.0,
                    'bioactivitySource': 'none',
                    'peptideType': 'prodomain',
                    'cleavedBy': 'PCSK5/6/7'
                })
                
                print(f"   📦 Prodomain: {len(prodomain_seq)} aa")
                print(f"      N-term motif: SIGNAL")
                print(f"      C-term motif: {site.motif}")
        
        print(f"✅ PCSK5/6/7 extracted: {len(peptides)} peptide(s)")
        
        return peptides
    
    @staticmethod
    def _extract_ultra_permissive(
        sequence: str,
        cleavage_sites: List[CleavageSite],
        signal_length: int
    ) -> List[Dict]:
        """
        Extraction ultra-permissive OPTIMISÉE
        
        ⭐ NOUVEAU : Capture motifs N et C terminal
        """
        peptides = []
        
        if len(cleavage_sites) == 0:
            return []
        
        print(f"\n🧬 Ultra-permissive extraction from {len(cleavage_sites)} sites...")
        
        sorted_sites = sorted(cleavage_sites, key=lambda s: s.index)
        
        MAX_DISTANCE = 100
        MAX_LENGTH = 50
        MIN_CONFIDENCE = 30
        
        for i in range(len(sorted_sites)):
            for j in range(i + 1, len(sorted_sites)):
                site_start = sorted_sites[i]
                site_end = sorted_sites[j]
                
                start_pos = site_start.position
                
                if 'RF' in site_end.motif or 'RY' in site_end.motif:
                    end_pos = site_end.position
                else:
                    end_pos = site_end.index
                
                peptide_length = end_pos - start_pos
                
                if peptide_length > MAX_DISTANCE or peptide_length > MAX_LENGTH:
                    continue
                
                if peptide_length < 3:
                    continue
                
                peptide_seq = sequence[start_pos:end_pos]
                
                confidence = PeptideExtractor._calculate_confidence(
                    peptide_seq,
                    site_start,
                    site_end,
                    peptide_length
                )
                
                if confidence < MIN_CONFIDENCE:
                    continue
                
                peptides.append({
                    'sequence': peptide_seq,
                    'start': start_pos + 1,
                    'end': end_pos,
                    'length': peptide_length,
                    'inRange': config.OPTIMAL_PEPTIDE_MIN_LENGTH <= peptide_length <= config.OPTIMAL_PEPTIDE_MAX_LENGTH,
                    'cleavageMotifN': site_start.motif,  # ⭐ NOUVEAU
                    'cleavageMotifC': site_end.motif,    # ⭐ NOUVEAU
                    'cleavageMotif': PeptideExtractor._get_cleavage_label(site_start, site_end),  # Compatibilité
                    'bioactivityScore': 0.0,
                    'bioactivitySource': 'none',
                    'confidenceScore': confidence,
                    'confidenceBadge': PeptideExtractor._get_confidence_badge(confidence)
                })
        
        peptides = PeptideExtractor._remove_overlapping_peptides(peptides)
        
        peptides.sort(key=lambda p: (
            p.get('confidenceScore', 0),
            -p['length']
        ), reverse=True)
        
        MAX_PEPTIDES = 50
        if len(peptides) > MAX_PEPTIDES:
            print(f"⚠️ Truncating from {len(peptides)} to top {MAX_PEPTIDES} peptides")
            peptides = peptides[:MAX_PEPTIDES]
        
        print(f"✅ Ultra-permissive extracted: {len(peptides)} peptides")
        if len(peptides) > 0:
            print(f"   Top confidence: {peptides[0].get('confidenceScore', 0)}")
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
        
        sorted_peptides = sorted(peptides, key=lambda p: p.get('confidenceScore', 0), reverse=True)
        
        filtered = []
        
        for peptide in sorted_peptides:
            is_overlapping = False
            
            for kept_peptide in filtered:
                overlap = PeptideExtractor._calculate_overlap(peptide, kept_peptide)
                
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
        
        if end1 <= start2 or end2 <= start1:
            return 0.0
        
        overlap_start = max(start1, start2)
        overlap_end = min(end1, end2)
        overlap_length = overlap_end - overlap_start
        
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
        """
        score = 0
        
        start_motif = site_start.motif
        if '...' in start_motif:
            score += 50
        elif start_motif in ['KK', 'KR', 'RR', 'RK']:
            score += 50
        else:
            score += 15
        
        end_motif = site_end.motif
        if '...' in end_motif:
            score += 50
        
        terminal_motif = PeptideExtractor._get_terminal_motif(peptide_seq)
        if terminal_motif in ['RF', 'RY', 'RFG', 'RYG']:
            score += 30
            print(f"      🎯 Terminal RFamide: +30 pts")
        elif terminal_motif == 'G':
            score += 15
        
        if 5 <= length <= 15:
            score += 20
        elif 15 < length <= 30:
            score += 10
        elif 30 < length <= 50:
            score += 5
        elif length > 100:
            score -= 30
        
        if '...' in start_motif or '...' in end_motif:
            score = max(score, 90)
            print(f"      🎯 RFamide bonus: minimum score 90")
        
        final_score = min(max(score, 0), 100)
        
        return final_score
    
    @staticmethod
    def _get_terminal_motif(sequence: str) -> str:
        """Détecte motif terminal spécial"""
        if len(sequence) < 2:
            return 'none'
        
        if len(sequence) >= 3:
            last_three = sequence[-3:]
            if last_three in ['RFG', 'RYG']:
                return last_three
        
        last_two = sequence[-2:]
        if last_two in ['RF', 'RY']:
            return last_two
        elif sequence[-1] == 'G':
            return 'G'
        
        return 'none'
    
    @staticmethod
    def _get_cleavage_label(site_start: CleavageSite, site_end: CleavageSite) -> str:
        """Génère label pour le motif de clivage (compatibilité)"""
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