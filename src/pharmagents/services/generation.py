from __future__ import annotations

from typing import Dict, List, Sequence

from pharmagents.utils.chem import basic_descriptors, brics_recombine, diversity_select, enumerate_analogs, score_druglikeness


class StrongMoleculeGenerator:
    def __init__(self, max_candidates: int = 128):
        self.max_candidates = max_candidates

    def generate(self, seeds: Sequence[str], library: Sequence[str], requirements: Dict[str, float | bool], target_name: str) -> List[str]:
        candidates = set(seeds)
        for smi in list(seeds)[:12]:
            candidates.update(enumerate_analogs(smi, limit=20))
        candidates.update(brics_recombine(list(set(list(seeds) + list(library)[:20])), limit=64))

        scored: Dict[str, float] = {}
        for smi in candidates:
            d = basic_descriptors(smi)
            if d.get('valid', 0.0) == 0.0:
                continue
            score = score_druglikeness(smi)
            if requirements.get('prefer_low_tpsa'):
                score += max(0.0, (90 - d['tpsa']) / 100)
            if requirements.get('avoid_cns'):
                score += max(0.0, (120 - d['logp'] * 15) / 120)
            if d['mw'] <= float(requirements.get('mw_max', 520.0)):
                score += 0.15
            score += 0.05 * {'JAK1': 1.0, 'JAK2': 0.9, 'JAK3': 0.95}.get(target_name, 0.7)
            scored[smi] = float(score)
        selected = diversity_select(list(scored.keys()), scored, k=min(48, self.max_candidates), max_similarity=0.72)
        return sorted(selected, key=lambda s: scored[s], reverse=True)
