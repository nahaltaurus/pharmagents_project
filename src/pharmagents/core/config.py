from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / 'data'
DB_PATH = DATA_DIR / 'experience.db'
LOG_DIR = ROOT.parent.parent / 'runs'
CACHE_DIR = ROOT.parent.parent / 'cache'
RECEPTOR_DIR = CACHE_DIR / 'receptors'


def _env_bool(name: str, default: bool = False) -> bool:
    value = os.getenv(name)
    if value is None:
        return default
    return value.strip().lower() in {'1', 'true', 'yes', 'y', 'on'}


@dataclass(slots=True)
class PipelineConfig:
    max_targets: int = int(os.getenv('PHARMAGENTS_MAX_TARGETS', '3'))
    n_leads: int = int(os.getenv('PHARMAGENTS_N_LEADS', '12'))
    optimization_rounds: int = int(os.getenv('PHARMAGENTS_OPT_ROUNDS', '4'))
    toxicity_threshold: float = float(os.getenv('PHARMAGENTS_TOX_THRESHOLD', '0.65'))
    accept_threshold: float = float(os.getenv('PHARMAGENTS_ACCEPT_THRESHOLD', '0.55'))
    use_real_fetch: bool = _env_bool('PHARMAGENTS_USE_REAL_FETCH', False)
    docking_backend: str = os.getenv('PHARMAGENTS_DOCKING_BACKEND', 'auto')
    log_dir: Path = Path(os.getenv('PHARMAGENTS_LOG_DIR', str(LOG_DIR)))
    cache_dir: Path = Path(os.getenv('PHARMAGENTS_CACHE_DIR', str(CACHE_DIR)))
    receptor_cache_dir: Path = Path(os.getenv('PHARMAGENTS_RECEPTOR_CACHE_DIR', str(RECEPTOR_DIR)))
    api_timeout_sec: int = int(os.getenv('PHARMAGENTS_API_TIMEOUT', '15'))
    max_generated_molecules: int = int(os.getenv('PHARMAGENTS_MAX_GENERATED', '128'))
    auto_prepare_receptors: bool = _env_bool('PHARMAGENTS_AUTO_PREPARE_RECEPTORS', True)
