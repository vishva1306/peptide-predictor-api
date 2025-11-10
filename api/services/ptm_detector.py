"""Service de détection des modifications post-traductionnelles (PTMs)"""
import re
from typing import List, Dict, Optional


class PTMDetector:
    """
    Détecte les modifications post-traductionnelles dans les peptides
    
    PTMs détectées :
    1. C-terminal amidation (peptide suivi de G[RK]{1,2})
    2. N-terminal pyroglutamate (Q/E → pGlu)
    3. Disulfide bonds (2+ Cys)
    4. Ghrelin acylation (GSSF motif)
    5. Tyrosine O-sulfation
    6. N-glycosylation (N-X-S/T)
    """
    
    @staticmethod
    def detect_all_ptms(
        peptide_sequence: str,
        full_protein_sequence: str = None,
        peptide_start: int = None,
        peptide_end: int = None
    ) -> List[Dict]:
        """
        Détecte toutes les PTMs dans une séquence
        
        Args:
            peptide_sequence: Séquence du peptide (ex: "SYSMEHFRWGKPVG")
            full_protein_sequence: Séquence complète de la protéine (optionnel)
            peptide_start: Position de début du peptide dans la protéine (1-indexed)
            peptide_end: Position de fin du peptide dans la protéine (1-indexed)
        
        Returns:
            Liste de PTMs détectées avec détails
        """
        ptms = []
        
        # 1. C-terminal amidation (besoin du contexte après le peptide)
        c_amid = PTMDetector.detect_c_terminal_amidation(
            peptide_sequence,
            full_protein_sequence,
            peptide_end
        )
        if c_amid:
            ptms.append(c_amid)
        
        # 2. N-terminal pyroglutamate
        n_pglu = PTMDetector.detect_n_terminal_pyroglu(peptide_sequence)
        if n_pglu:
            ptms.append(n_pglu)
        
        # 3. Disulfide bonds
        disulfide = PTMDetector.detect_disulfide_bonds(peptide_sequence)
        if disulfide:
            ptms.append(disulfide)
        
        # 4. Ghrelin acylation
        ghrelin = PTMDetector.detect_ghrelin_acylation(peptide_sequence)
        if ghrelin:
            ptms.append(ghrelin)
        
        # 5. Tyrosine O-sulfation
        y_sulf = PTMDetector.detect_tyrosine_sulfation(peptide_sequence)
        if y_sulf:
            ptms.extend(y_sulf)
        
        # 6. N-glycosylation
        n_glyc = PTMDetector.detect_n_glycosylation(peptide_sequence)
        if n_glyc:
            ptms.extend(n_glyc)
        
        return ptms
    
    @staticmethod
    def detect_c_terminal_amidation(
        peptide_sequence: str,
        full_protein_sequence: str = None,
        peptide_end: int = None
    ) -> Optional[Dict]:
        """
        Détecte C-terminal amidation
        
        Logique :
        - Le peptide se termine par G
        - APRÈS le peptide dans la protéine, il y a [RK]{1,2}
        - Ces résidus sont clivés par PCSK1/3
        - Le G terminal est converti en -NH₂ par PAM
        
        Args:
            peptide_sequence: Séquence du peptide
            full_protein_sequence: Séquence complète de la protéine
            peptide_end: Position de fin du peptide (1-indexed)
        
        Returns:
            PTM si détectée, None sinon
        """
        # Si pas de contexte protéine, pas de détection possible
        if not full_protein_sequence or peptide_end is None:
            return None
        
        # Vérifier que peptide_end est un int
        if not isinstance(peptide_end, int):
            return None
        
        # Vérifier que le peptide se termine par G
        if not peptide_sequence or not peptide_sequence.endswith('G'):
            return None
        
        # Convertir peptide_end (1-indexed) en index Python (0-indexed)
        after_peptide_idx = peptide_end
        
        # Sécurité : vérifier qu'on ne dépasse pas la longueur
        if after_peptide_idx >= len(full_protein_sequence):
            return None
        
        # Extraire les résidus après le peptide (max 3 aa)
        after_peptide = full_protein_sequence[after_peptide_idx:after_peptide_idx + 3]
        
        # Patterns : [RK]{1,2}
        patterns = [
            (r'^RR', 'GRR'),
            (r'^RK', 'GRK'),
            (r'^KR', 'GKR'),
            (r'^KK', 'GKK'),
            (r'^R', 'GR'),
            (r'^K', 'GK'),
        ]
        
        for pattern, motif in patterns:
            if re.match(pattern, after_peptide):
                return {
                    'type': 'C-terminal amidation',
                    'shortName': 'C-amidation',
                    'emoji': '🔵',
                    'enzyme': 'PAM',
                    'motif': motif,
                    'position': 'C-terminus',
                    'description': f'{motif} → -NH₂',
                    'removes_g': True
                }
        
        return None
    
    @staticmethod
    def detect_n_terminal_pyroglu(sequence: str) -> Optional[Dict]:
        """
        Détecte N-terminal pyroglutamate
        Motif : Q ou E au N-terminus
        Enzyme : QPCT (Q) ou QPCTL (E)
        """
        if not sequence:
            return None
        
        first_aa = sequence[0]
        
        if first_aa == 'Q':
            return {
                'type': 'N-terminal pyroglutamate',
                'shortName': 'N-pGlu',
                'emoji': '🟢',
                'enzyme': 'QPCT',
                'residue': 'Q',
                'position': 1,
                'description': 'Q → pGlu'
            }
        elif first_aa == 'E':
            return {
                'type': 'N-terminal pyroglutamate',
                'shortName': 'N-pGlu',
                'emoji': '🟢',
                'enzyme': 'QPCTL',
                'residue': 'E',
                'position': 1,
                'description': 'E → pGlu'
            }
        return None
    
    @staticmethod
    def detect_disulfide_bonds(sequence: str) -> Optional[Dict]:
        """
        Détecte potentiel de ponts disulfures
        Motif : 2+ Cystéines (C)
        Enzyme : PDI, ER oxidoreductases
        """
        cys_positions = [i + 1 for i, aa in enumerate(sequence) if aa == 'C']
        
        if len(cys_positions) >= 2:
            return {
                'type': 'Disulfide bonds',
                'shortName': 'Disulfide',
                'emoji': '🔴',
                'enzyme': 'PDI / ER oxidoreductases',
                'positions': cys_positions,
                'count': len(cys_positions) // 2,
                'description': f'{len(cys_positions)} Cys (≥{len(cys_positions)//2} bonds)'
            }
        return None
    
    @staticmethod
    def detect_ghrelin_acylation(sequence: str) -> Optional[Dict]:
        """
        Détecte acylation spécifique ghrelin
        Motif : ^GSSF avec Ser3
        Enzyme : GOAT (MBOAT4)
        """
        if sequence.startswith('GSSF'):
            return {
                'type': 'Ghrelin acylation',
                'shortName': 'Ghrelin-acyl',
                'emoji': '🟣',
                'enzyme': 'GOAT (MBOAT4)',
                'residue': 'Ser3',
                'position': 3,
                'description': 'Ser3 octanoylation'
            }
        return None
    
    @staticmethod
    def detect_tyrosine_sulfation(sequence: str) -> List[Dict]:
        """
        Détecte Tyrosine O-sulfation
        Motif : Y dans contexte acide (D/E dans ±5aa)
        Enzyme : TPST1/TPST2
        """
        sulfations = []
        
        for i, aa in enumerate(sequence):
            if aa == 'Y':
                # Extraire fenêtre ±5 résidus
                start = max(0, i - 5)
                end = min(len(sequence), i + 6)
                window = sequence[start:end]
                
                # Compter résidus acides
                acidic_count = window.count('D') + window.count('E')
                
                # Au moins 2 résidus acides dans la fenêtre
                if acidic_count >= 2:
                    sulfations.append({
                        'type': 'Tyrosine O-sulfation',
                        'shortName': 'Y-sulfation',
                        'emoji': '🟡',
                        'enzyme': 'TPST1/TPST2',
                        'residue': f'Y{i + 1}',
                        'position': i + 1,
                        'description': f'Y{i + 1} → Y(SO₃)'
                    })
        
        return sulfations
    
    @staticmethod
    def detect_n_glycosylation(sequence: str) -> List[Dict]:
        """
        Détecte N-glycosylation
        Motif : N-X-[ST] où X ≠ P
        Enzyme : Oligosaccharyltransferase
        """
        glycosylations = []
        
        # Pattern : N-X-[ST] où X n'est pas P
        pattern = r'N[^P][ST]'
        
        for match in re.finditer(pattern, sequence):
            start_pos = match.start()
            motif = match.group()
            
            glycosylations.append({
                'type': 'N-glycosylation',
                'shortName': 'N-glyco',
                'emoji': '🟠',
                'enzyme': 'Oligosaccharyltransferase',
                'motif': motif,
                'position': start_pos + 1,
                'description': f'N{start_pos + 1} glycosylation'
            })
        
        return glycosylations
    
    @staticmethod
    def generate_modified_sequence(sequence: str, ptms: List[Dict]) -> str:
        """
        Génère la séquence modifiée avec notation PTM
        
        Args:
            sequence: Séquence originale
            ptms: Liste des PTMs détectées
        
        Returns:
            Séquence modifiée avec annotations
        """
        if not ptms:
            return sequence
        
        # Guard : vérifier que sequence est un string
        if not isinstance(sequence, str):
            return str(sequence)
        
        modified = list(sequence)
        
        # ⭐ ÉTAPE 1 : Vérifier si C-terminal amidation et enlever le G
        has_c_amidation = any(
            ptm.get('type') == 'C-terminal amidation' and ptm.get('removes_g', False)
            for ptm in ptms
        )
        
        if has_c_amidation and modified and modified[-1] == 'G':
            modified = modified[:-1]
            print(f"🔵 C-amidation: G terminal enlevé → séquence devient {''.join(modified)}")
        
        # ⭐ ÉTAPE 2 : Trier les PTMs par position (pour ordre d'application)
        sorted_ptms = sorted(
            [p for p in ptms if isinstance(p.get('position'), int)],
            key=lambda x: x.get('position', 0)
        )
        
        # ⭐ ÉTAPE 3 : Appliquer les modifications avec tracking d'offset
        # Tracker les cystéines déjà modifiées pour disulfide bonds
        modified_cys_indices = set()
        
        for ptm in ptms:
            ptm_type = ptm.get('type', '')
            
            if ptm_type == 'C-terminal amidation':
                # Ajouter -NH₂ au C-terminus (le G a déjà été enlevé)
                modified.append('-NH₂')
                print(f"🔵 C-amidation: -NH₂ ajouté")
            
            elif ptm_type == 'N-terminal pyroglutamate':
                # Remplacer Q ou E par pGlu
                if len(modified) > 0:
                    old_aa = modified[0]
                    modified[0] = 'pGlu'
                    print(f"🟢 N-pGlu: {old_aa} → pGlu au N-terminus")
            
            elif ptm_type == 'Ghrelin acylation':
                # Ajouter octanoyl sur G au début
                if len(modified) > 0 and modified[0] == 'G':
                    modified[0] = 'G(C8:0)'
                    print(f"🟣 Ghrelin: G → G(C8:0)")
            
            elif ptm_type == 'Disulfide bonds':
                # Numéroter toutes les cystéines
                positions = ptm.get('positions', [])
                print(f"🔴 Disulfide: Numérotation de {len(positions)} cystéines aux positions {positions}")
                
                cys_found = 0
                for i, aa in enumerate(modified):
                    if aa == 'C' and i not in modified_cys_indices:
                        cys_found += 1
                        modified[i] = f'C{cys_found}'
                        modified_cys_indices.add(i)
                        if cys_found >= len(positions):
                            break
            
            elif ptm_type == 'Tyrosine O-sulfation':
                # Trouver la position de la tyrosine
                pos = ptm.get('position', 0) - 1  # Convertir en 0-indexed
                
                if 0 <= pos < len(modified):
                    # Vérifier si c'est bien une Y
                    if modified[pos] == 'Y':
                        modified[pos] = 'Y(SO₃)'
                        print(f"🟡 Y-sulfation: Y{pos+1} → Y(SO₃)")
            
            elif ptm_type == 'N-glycosylation':
                # Trouver la position de l'asparagine
                pos = ptm.get('position', 0) - 1  # Convertir en 0-indexed
                
                if 0 <= pos < len(modified):
                    # Vérifier si c'est bien une N
                    if modified[pos] == 'N':
                        modified[pos] = 'N(GlcNAc)'
                        print(f"🟠 N-glyco: N{pos+1} → N(GlcNAc)")
        
        result = ''.join(modified)
        print(f"✅ Séquence finale modifiée : {result}")
        return result


# Instance globale
ptm_detector = PTMDetector()