from __future__ import annotations

from typing import List, Optional

from fastapi import FastAPI
from pydantic import BaseModel
from pharmagents.core.config import PipelineConfig
from pharmagents.services.experience_db import ExperienceDB
from pharmagents.workflows.pipeline import VirtualPharmaPipeline

app = FastAPI(title='PharmAgents API', version='0.3.0')
pipeline = VirtualPharmaPipeline()
exp_db = ExperienceDB()


class RunRequest(BaseModel):
    disease: str
    top_k_targets: int = 3
    n_leads: int = 12
    optimization_rounds: int = 4
    real_fetch: bool = False
    docking_backend: str = 'auto'
    receptor_path: Optional[str] = None
    pdb_path: Optional[str] = None
    pocket_center: Optional[List[float]] = None
    box_size: Optional[List[float]] = None


@app.get('/health')
def health():
    return {'status': 'ok'}


@app.post('/run')
def run_pipeline(req: RunRequest):
    pipe = VirtualPharmaPipeline(config=PipelineConfig(use_real_fetch=req.real_fetch, docking_backend=req.docking_backend))
    return pipe.run(
        disease=req.disease,
        top_k_targets=req.top_k_targets,
        n_leads=req.n_leads,
        optimization_rounds=req.optimization_rounds,
        target_overrides={
            'receptor_path': req.receptor_path,
            'pdb_path': req.pdb_path,
            'pocket_center': req.pocket_center,
            'box_size': req.box_size,
        },
    )


@app.get('/history')
def history(limit: int = 20):
    return {'history': exp_db.history(limit=limit)}


@app.get('/runs')
def runs(limit: int = 20):
    return {'runs': pipeline.experiment_logger.list_runs(limit=limit)}
