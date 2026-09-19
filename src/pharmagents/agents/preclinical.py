from __future__ import annotations

from typing import List
from pharmagents.core.config import PipelineConfig
from pharmagents.core.schemas import EvaluationReport, MoleculeRecord, TargetCandidate
from pharmagents.services.docking import DockingService
from pharmagents.services.experience_db import ExperienceDB
from pharmagents.services.knowledge_base import KnowledgeBase
from pharmagents.utils.chem import basic_descriptors, nearest_neighbors, sa_like_score, tanimoto


class PreclinicalEvaluationAgent:
    def __init__(self, kb: KnowledgeBase, exp_db: ExperienceDB, config: PipelineConfig | None = None):
        self.kb = kb
        self.exp_db = exp_db
        self.config = config or PipelineConfig()
        self.docking = DockingService(self.config.docking_backend)

    def _toxicity_risk(self, smiles: str) -> float:
        refs = [r['smiles'] for r in self.kb.tox_reference()]
        nn = nearest_neighbors(smiles, refs, top_k=5)
        if not nn:
            return 0.5
        ref_map = {r['smiles']: r for r in self.kb.tox_reference()}
        weighted = 0.0
        denom = 0.0
        for s, sim in nn:
            risk = ref_map[s]['risk']
            weighted += sim * risk
            denom += sim
        return float(weighted / max(denom, 1e-6))

    def _metabolism_risk(self, smiles: str) -> float:
        d = basic_descriptors(smiles)
        aromatic_heavy = 0.15 if d.get('rings', 0) >= 3 else 0.0
        logp_risk = 0.2 if d.get('logp', 0) > 4 else 0.0
        return min(1.0, 0.2 + aromatic_heavy + logp_risk)

    def _novelty(self, smiles: str, target_name: str) -> float:
        hist = self.exp_db.fetch_for_target(target_name, limit=30)
        if not hist:
            return 1.0
        max_sim = max((tanimoto(smiles, h['output_smiles']) for h in hist), default=0.0)
        return float(1.0 - max_sim)

    def run(self, target: TargetCandidate, molecules: List[MoleculeRecord]) -> List[EvaluationReport]:
        reports = []
        for m in molecules:
            d = basic_descriptors(m.smiles)
            tox = self._toxicity_risk(m.smiles)
            met = self._metabolism_risk(m.smiles)
            sa = sa_like_score(m.smiles)
            novelty = self._novelty(m.smiles, target.target_name)
            docking = self.docking.score(m.smiles, target)
            docking_bonus = max(0.0, min(1.0, (-docking['affinity'] - 4.5) / 4.0))
            uncertainty = float(0.4 * abs(tox - 0.5) + 0.3 * abs(sa - 0.5) + 0.3 * abs(docking_bonus - 0.5))
            synth = float(100 * sa)
            portfolio = float((d.get('qed', 0) * 0.25) + (sa * 0.20) + ((1 - tox) * 0.18) + ((1 - met) * 0.10) + (novelty * 0.10) + ((1 - uncertainty) * 0.07) + (docking_bonus * 0.10))
            recommended = portfolio >= self.config.accept_threshold and tox < self.config.toxicity_threshold and sa > 0.45
            rationale = (
                f"QED-like={d.get('qed',0):.3f}, SA-like={sa:.3f}, toxicity={tox:.3f}, metabolism={met:.3f}, "
                f"novelty={novelty:.3f}, docking={docking['affinity']:.3f} ({docking['backend']}), uncertainty={uncertainty:.3f}. "
                f"{'Recommend for wet-lab prioritization.' if recommended else 'Do not prioritize yet.'}"
            )
            reports.append(EvaluationReport(
                smiles=m.smiles,
                docking_like_score=float(docking['affinity']),
                docking_backend=docking['backend'],
                qed_like_score=d.get('qed', 0.0),
                sa_like_score=sa,
                toxicity_risk=tox,
                metabolism_risk=met,
                synthesizability_confidence=synth,
                uncertainty=uncertainty,
                novelty_score=novelty,
                portfolio_score=portfolio,
                recommended=recommended,
                rationale=rationale,
                metadata={'docking': docking, 'source_metadata': m.metadata},
            ))
        reports.sort(key=lambda x: x.portfolio_score, reverse=True)
        return reports
