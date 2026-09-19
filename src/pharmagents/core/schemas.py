from __future__ import annotations

from typing import Any, Dict, List, Optional
from pydantic import BaseModel, Field


class TargetCandidate(BaseModel):
    disease: str
    target_name: str
    uniprot_id: str
    pdb_id: str
    pocket_center: List[float]
    confidence: float
    rationale: str
    source: str = 'local_kb'
    metadata: Dict[str, Any] = Field(default_factory=dict)


class MoleculeRecord(BaseModel):
    smiles: str
    source: str
    score: float = 0.0
    descriptors: Dict[str, float] = Field(default_factory=dict)
    rationale: str = ''
    metadata: Dict[str, Any] = Field(default_factory=dict)


class EvaluationReport(BaseModel):
    smiles: str
    docking_like_score: float
    docking_backend: str = 'heuristic'
    qed_like_score: float
    sa_like_score: float
    toxicity_risk: float
    metabolism_risk: float
    synthesizability_confidence: float
    uncertainty: float
    novelty_score: float
    portfolio_score: float
    recommended: bool
    rationale: str
    metadata: Dict[str, Any] = Field(default_factory=dict)


class PipelineOutput(BaseModel):
    disease: str
    targets: List[TargetCandidate]
    initial_leads: List[MoleculeRecord]
    optimized_leads: List[MoleculeRecord]
    evaluations: List[EvaluationReport]
    final_report: str
    extension_summary: Dict[str, Any]
    run_id: Optional[str] = None
    artifacts: Dict[str, Any] = Field(default_factory=dict)
