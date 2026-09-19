from __future__ import annotations

from typing import List, Tuple

from pharmagents.core.config import PipelineConfig
from pharmagents.core.schemas import MoleculeRecord, TargetCandidate
from pharmagents.services.docking import DockingService
from pharmagents.services.experience_db import ExperienceDB
from pharmagents.utils.chem import basic_descriptors, enumerate_analogs, tanimoto


class LeadOptimizationAgent:
    def __init__(self, exp_db: ExperienceDB, config: PipelineConfig | None = None):
        self.exp_db = exp_db
        self.config = config or PipelineConfig()
        self.docking = DockingService(self.config.docking_backend)

    def _interaction_score(self, smiles: str, target: TargetCandidate) -> Tuple[float, dict]:
        dock = self.docking.score(smiles, target)
        d = basic_descriptors(smiles)
        if d.get('valid', 0.0) == 0.0:
            return 0.0, dock
        docking_bonus = max(0.0, min(1.0, (-dock['affinity'] - 4.5) / 4.0))
        score = 0.50 * d['qed'] + 0.20 * docking_bonus + 0.20 * (1.0 / (1.0 + abs(d['logp'] - 2.2))) + 0.10 * (1.0 if d['mw'] < 550 and d['hba'] <= 10 and d['hbd'] <= 5 else 0.5)
        return float(score), dock

    def _reflection(self, old_s: str, new_s: str, old_score: float, new_score: float, dock_backend: str) -> str:
        sim = tanimoto(old_s, new_s)
        trend = 'improved' if new_score > old_score else 'did not improve'
        return f'Optimization {trend}; docking_backend={dock_backend}; structural_similarity={sim:.3f}; old_score={old_score:.3f}; new_score={new_score:.3f}.'

    def optimize_one(self, disease: str, target: TargetCandidate, mol: MoleculeRecord, rounds: int) -> MoleculeRecord:
        best_smiles = mol.smiles
        best_score, best_dock = self._interaction_score(best_smiles, target)
        prior = self.exp_db.fetch_for_target(target.target_name, limit=5)
        memory_bonus = 0.03 if prior else 0.0
        history = []

        for _ in range(rounds):
            candidates = enumerate_analogs(best_smiles, limit=20)
            if not candidates:
                break
            scored: List[Tuple[str, float, dict]] = []
            for c in candidates:
                s, dock = self._interaction_score(c, target)
                scored.append((c, s + memory_bonus, dock))
            scored.sort(key=lambda x: x[1], reverse=True)
            cand_smi, cand_score, cand_dock = scored[0]
            history.append({'candidate': cand_smi, 'score': cand_score, 'dock': cand_dock})
            if cand_score > best_score:
                note = self._reflection(best_smiles, cand_smi, best_score, cand_score, cand_dock['backend'])
                self.exp_db.add(disease, target.target_name, best_smiles, cand_smi, cand_score, {'reflection': note, 'dock': cand_dock})
                best_smiles, best_score, best_dock = cand_smi, cand_score, cand_dock

        return MoleculeRecord(
            smiles=best_smiles,
            source='optimized',
            score=best_score,
            descriptors=basic_descriptors(best_smiles),
            rationale=f'Optimized over {rounds} rounds with docking-aware interaction analysis, analog generation, and reflection memory.',
            metadata={'target': target.target_name, 'dock': best_dock, 'history': history[-5:]},
        )

    def run(self, disease: str, target: TargetCandidate, leads: List[MoleculeRecord], rounds: int = 4) -> List[MoleculeRecord]:
        return [self.optimize_one(disease, target, m, rounds) for m in leads]
