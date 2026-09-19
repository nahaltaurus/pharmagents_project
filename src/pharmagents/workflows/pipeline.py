from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional

from pharmagents.agents.lead_identification import LeadIdentificationAgent
from pharmagents.agents.lead_optimization import LeadOptimizationAgent
from pharmagents.agents.preclinical import PreclinicalEvaluationAgent
from pharmagents.agents.target_discovery import TargetDiscoveryAgent
from pharmagents.core.config import PipelineConfig
from pharmagents.core.schemas import PipelineOutput
from pharmagents.services.docking import DockingService
from pharmagents.services.experience_db import ExperienceDB
from pharmagents.services.experiment_logger import ExperimentLogger
from pharmagents.services.knowledge_base import KnowledgeBase
from pharmagents.utils.io import save_json


class VirtualPharmaPipeline:
    def __init__(self, config: PipelineConfig | None = None):
        self.config = config or PipelineConfig()
        self.kb = KnowledgeBase()
        self.exp_db = ExperienceDB()
        self.experiment_logger = ExperimentLogger(self.config)
        self.docking_service = DockingService(self.config.docking_backend, self.config)
        self.target_agent = TargetDiscoveryAgent(self.kb, self.config)
        self.lead_agent = LeadIdentificationAgent(self.kb, self.config)
        self.opt_agent = LeadOptimizationAgent(self.exp_db, self.config)
        self.eval_agent = PreclinicalEvaluationAgent(self.kb, self.exp_db, self.config)

    def run(
        self,
        disease: str,
        top_k_targets: int = 3,
        n_leads: int = 12,
        optimization_rounds: int = 4,
        target_overrides: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        config_payload = {
            'top_k_targets': top_k_targets,
            'n_leads': n_leads,
            'optimization_rounds': optimization_rounds,
            'use_real_fetch': self.config.use_real_fetch,
            'docking_backend': self.config.docking_backend,
            'target_overrides': target_overrides or {},
        }
        run_id = self.experiment_logger.start_run(disease, config_payload)

        targets = self.target_agent.run(disease, top_k=top_k_targets)
        if target_overrides:
            for t in targets:
                t.metadata.update({k: v for k, v in target_overrides.items() if v is not None})
                if target_overrides.get('pocket_center'):
                    t.pocket_center = list(target_overrides['pocket_center'])
                if target_overrides.get('box_size'):
                    t.metadata['box_size'] = list(target_overrides['box_size'])
        if self.config.auto_prepare_receptors:
            targets = [self.docking_service.ensure_target_receptor(t) for t in targets]

        if not targets:
            out = {
                'disease': disease,
                'targets': [],
                'initial_leads': [],
                'optimized_leads': [],
                'evaluations': [],
                'final_report': f"No target candidates were found for disease '{disease}'. Add more knowledge-base entries.",
                'extension_summary': {},
                'run_id': run_id,
                'artifacts': {},
            }
            self._persist_run(run_id, disease, config_payload, out)
            return out

        self.experiment_logger.log_event(run_id, 'targets_selected', {'count': len(targets), 'targets': [t.model_dump() for t in targets]})
        primary_target = targets[0]
        leads = self.lead_agent.run(disease, primary_target, n_leads=n_leads)
        self.experiment_logger.log_event(run_id, 'lead_identification_complete', {'count': len(leads)})
        optimized = self.opt_agent.run(disease, primary_target, leads, rounds=optimization_rounds)
        self.experiment_logger.log_event(run_id, 'lead_optimization_complete', {'count': len(optimized)})
        evaluations = self.eval_agent.run(primary_target, optimized)
        top = evaluations[:5]
        accepted = [r for r in evaluations if r.recommended]
        extension_summary = {
            'paper_extension': 'uncertainty-aware active-learning prioritization and portfolio-diverse ranking',
            'top_active_learning_candidates': [r.smiles for r in top],
            'recommended_count': len(accepted),
            'docking_backend_used': top[0].docking_backend if top else 'heuristic',
            'prepared_receptor_path': primary_target.metadata.get('receptor_path'),
            'downloaded_pdb_path': primary_target.metadata.get('pdb_path'),
        }

        final_report = self._compose_report(disease, targets, evaluations, extension_summary)
        out = PipelineOutput(
            disease=disease,
            targets=targets,
            initial_leads=leads,
            optimized_leads=optimized,
            evaluations=evaluations,
            final_report=final_report,
            extension_summary=extension_summary,
            run_id=run_id,
            artifacts={},
        ).model_dump()
        artifacts = self._persist_run(run_id, disease, config_payload, out)
        out['artifacts'] = artifacts
        return out

    def _persist_run(self, run_id: str, disease: str, config_payload: Dict[str, Any], out: Dict[str, Any]) -> Dict[str, Any]:
        run_dir = Path(self.config.log_dir) / run_id
        report_path = run_dir / 'report.json'
        save_json(report_path, out)
        metrics = {
            'target_count': len(out.get('targets', [])),
            'lead_count': len(out.get('initial_leads', [])),
            'optimized_count': len(out.get('optimized_leads', [])),
            'recommended_count': out.get('extension_summary', {}).get('recommended_count', 0),
        }
        artifacts = {'report_json': str(report_path), 'run_dir': str(run_dir)}
        self.experiment_logger.finalize_run(run_id, disease, config_payload, metrics, artifacts)
        return artifacts

    def _compose_report(self, disease: str, targets, evaluations, extension_summary) -> str:
        target_text = "\n".join(
            [f"- {t.target_name} ({t.uniprot_id}, {t.pdb_id}) confidence={t.confidence:.2f} source={t.source}" for t in targets]
        )
        top = evaluations[:5]
        mol_text = "\n".join(
            [
                f"- {r.smiles} | portfolio={r.portfolio_score:.3f} | docking={r.docking_like_score:.3f} ({r.docking_backend}) | tox={r.toxicity_risk:.3f} | SA={r.sa_like_score:.3f} | novelty={r.novelty_score:.3f} | recommended={r.recommended}"
                for r in top
            ]
        )
        return f"""
PharmAgents Working Project Report
=================================
Disease: {disease}

Selected Targets
----------------
{target_text}

Top Molecules After Optimization + PCC Evaluation
-------------------------------------------------
{mol_text}

Extension Beyond Paper
----------------------
{extension_summary['paper_extension']}
Top active-learning candidates: {', '.join(extension_summary['top_active_learning_candidates'])}
Recommended count: {extension_summary['recommended_count']}
Docking backend used: {extension_summary['docking_backend_used']}
Prepared receptor path: {extension_summary.get('prepared_receptor_path')}
Downloaded PDB path: {extension_summary.get('downloaded_pdb_path')}
"""
