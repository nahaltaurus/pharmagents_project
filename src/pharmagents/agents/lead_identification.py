from __future__ import annotations

from typing import Dict, List

from pharmagents.core.config import PipelineConfig
from pharmagents.core.schemas import MoleculeRecord, TargetCandidate
from pharmagents.services.generation import StrongMoleculeGenerator
from pharmagents.services.knowledge_base import KnowledgeBase
from pharmagents.utils.chem import basic_descriptors, diversity_select, score_druglikeness


class LeadIdentificationAgent:
    def __init__(self, kb: KnowledgeBase, config: PipelineConfig | None = None):
        self.kb = kb
        self.config = config or PipelineConfig()
        self.generator = StrongMoleculeGenerator(max_candidates=self.config.max_generated_molecules)

    def _requirements(self, disease: str, target: TargetCandidate) -> Dict[str, float | bool]:
        d = disease.lower()
        return {
            'prefer_low_tpsa': any(x in d for x in ['brain', 'parkinson', 'alzheimer']),
            'avoid_cns': any(x in d for x in ['asthma', 'dermatitis', 'eczema', 'colitis']),
            'anti_inflammatory': any(x in d for x in ['eczema', 'dermatitis', 'asthma', 'colitis', 'psoriasis']),
            'mw_max': 520.0,
        }

    def _score(self, smiles: str, source_row: Dict[str, object], target: TargetCandidate, req: Dict[str, float | bool]) -> float:
        d = basic_descriptors(smiles)
        if d.get('valid', 0.0) == 0.0:
            return 0.0
        score = score_druglikeness(smiles)
        score += 0.25 * float(source_row.get('affinity_prior', 0.5))
        score += 0.15 * float(source_row.get('target_bias', {}).get(target.target_name, 0.3))
        if req['prefer_low_tpsa']:
            score += max(0.0, (90 - d['tpsa']) / 100)
        if req['avoid_cns']:
            score += max(0.0, (120 - d['logp'] * 15) / 120)
        score += max(0.0, (float(req['mw_max']) - d['mw']) / float(req['mw_max'])) * 0.2
        return float(score)

    def run(self, disease: str, target: TargetCandidate, n_leads: int = 12) -> List[MoleculeRecord]:
        req = self._requirements(disease, target)
        rows = self.kb.get_molecule_library()
        seed_smiles = [r['smiles'] for r in rows if target.target_name in r.get('target_bias', {})][:12]
        generated = self.generator.generate(seed_smiles, [r['smiles'] for r in rows], req, target.target_name)

        raw: List[MoleculeRecord] = []
        for row in rows:
            smi = row['smiles']
            d = basic_descriptors(smi)
            if d.get('valid', 0.0) == 0.0:
                continue
            score = self._score(smi, row, target, req)
            raw.append(MoleculeRecord(smiles=smi, source=row['source'], score=score, descriptors=d, rationale='Selected from local library using target bias and medicinal chemistry filters.', metadata={'name': row.get('name', 'unknown')}))

        for smi in generated:
            d = basic_descriptors(smi)
            source_row = {'affinity_prior': 0.55, 'target_bias': {target.target_name: 0.55}}
            score = self._score(smi, source_row, target, req)
            raw.append(MoleculeRecord(smiles=smi, source='strong_generator', score=score, descriptors=d, rationale='Generated with BRICS recombination and analog enumeration conditioned on disease/target requirements.', metadata={'generator': 'strong_local'}))

        best_per_smiles = {}
        for rec in raw:
            if rec.smiles not in best_per_smiles or rec.score > best_per_smiles[rec.smiles].score:
                best_per_smiles[rec.smiles] = rec
        dedup = list(best_per_smiles.values())
        dedup.sort(key=lambda x: x.score, reverse=True)
        smiles = diversity_select([r.smiles for r in dedup], {r.smiles: r.score for r in dedup}, k=n_leads, max_similarity=0.72)
        return [r for r in dedup if r.smiles in smiles][:n_leads]
