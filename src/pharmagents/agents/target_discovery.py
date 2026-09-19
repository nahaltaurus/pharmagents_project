from __future__ import annotations

from typing import List

from pharmagents.core.config import PipelineConfig
from pharmagents.core.schemas import TargetCandidate
from pharmagents.services.knowledge_base import KnowledgeBase
from pharmagents.services.uniprot_rcsb import UniProtRCSBClient


class TargetDiscoveryAgent:
    def __init__(self, kb: KnowledgeBase, config: PipelineConfig | None = None):
        self.kb = kb
        self.config = config or PipelineConfig()
        self.remote = UniProtRCSBClient(self.config)

    def run(self, disease: str, top_k: int = 3) -> List[TargetCandidate]:
        local = self._local_candidates(disease)
        if self.config.use_real_fetch and local:
            enriched = [self._enrich_remote(c) for c in local]
        else:
            enriched = local
        enriched.sort(key=lambda x: x.confidence, reverse=True)
        return enriched[:top_k]

    def _local_candidates(self, disease: str) -> List[TargetCandidate]:
        matches = self.kb.disease_candidates(disease)
        candidates: List[TargetCandidate] = []
        for row in matches:
            for t in row['targets']:
                rationale = (
                    f"Target {t['target_name']} selected because disease match '{row['name']}' shares immune/pathology context "
                    f"with '{disease}', and structure {t['pdb_id']} has a defined pocket suitable for small-molecule work."
                )
                candidates.append(
                    TargetCandidate(
                        disease=disease,
                        target_name=t['target_name'],
                        uniprot_id=t['uniprot_id'],
                        pdb_id=t['pdb_id'],
                        pocket_center=t['pocket_center'],
                        confidence=float(t['confidence']),
                        rationale=rationale,
                        source='local_kb',
                        metadata={'matched_disease': row['name'], 'tags': row.get('tags', []), 'box_size': [20, 20, 20]},
                    )
                )
        return candidates

    def _enrich_remote(self, candidate: TargetCandidate) -> TargetCandidate:
        entry = self.remote.fetch_uniprot_entry(candidate.uniprot_id)
        pdb_hits = self.remote.search_pdb_for_uniprot(candidate.uniprot_id, rows=3)
        metadata = dict(candidate.metadata)
        metadata['uniprot_entry'] = entry
        selected_pdb = candidate.pdb_id
        if pdb_hits:
            best = pdb_hits[0]
            selected_pdb = best.get('pdb_id', candidate.pdb_id)
            metadata['rcsb_best_hit'] = best
            candidate.confidence = min(0.99, candidate.confidence + 0.02)
            candidate.source = 'local_kb+remote_fetch'
            candidate.rationale += f" Remote validation found PDB {selected_pdb} ({best.get('title', 'structure')})."
        pdb_path = self.remote.download_pdb_structure(selected_pdb)
        if pdb_path:
            metadata['pdb_path'] = pdb_path
            candidate.rationale += ' Structure file downloaded locally for optional docking.'
        candidate.pdb_id = selected_pdb
        candidate.metadata = metadata
        return candidate
