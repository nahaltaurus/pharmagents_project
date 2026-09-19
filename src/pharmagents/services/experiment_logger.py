from __future__ import annotations

import json
import sqlite3
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

from pharmagents.core.config import PipelineConfig


class ExperimentLogger:
    def __init__(self, config: PipelineConfig | None = None):
        self.config = config or PipelineConfig()
        self.root = Path(self.config.log_dir)
        self.root.mkdir(parents=True, exist_ok=True)
        self.db_path = self.root / 'experiments.db'
        self._init_db()

    def _conn(self):
        return sqlite3.connect(self.db_path)

    def _init_db(self) -> None:
        with self._conn() as conn:
            conn.execute(
                '''
                CREATE TABLE IF NOT EXISTS runs(
                    run_id TEXT PRIMARY KEY,
                    created_at TEXT NOT NULL,
                    disease TEXT NOT NULL,
                    config_json TEXT NOT NULL,
                    metrics_json TEXT NOT NULL,
                    artifacts_json TEXT NOT NULL
                )
                '''
            )
            conn.commit()

    def start_run(self, disease: str, config: Dict[str, Any]) -> str:
        run_id = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '_' + uuid.uuid4().hex[:8]
        run_dir = self.root / run_id
        run_dir.mkdir(parents=True, exist_ok=True)
        (run_dir / 'events.jsonl').write_text('', encoding='utf-8')
        self._write_event(run_id, {'type': 'run_started', 'disease': disease, 'config': config})
        return run_id

    def log_event(self, run_id: str, event_type: str, payload: Dict[str, Any]) -> None:
        self._write_event(run_id, {'type': event_type, **payload})

    def finalize_run(self, run_id: str, disease: str, config: Dict[str, Any], metrics: Dict[str, Any], artifacts: Dict[str, Any]) -> None:
        self._write_event(run_id, {'type': 'run_finished', 'metrics': metrics, 'artifacts': artifacts})
        with self._conn() as conn:
            conn.execute(
                'INSERT OR REPLACE INTO runs(run_id, created_at, disease, config_json, metrics_json, artifacts_json) VALUES(?,?,?,?,?,?)',
                (run_id, datetime.now(timezone.utc).isoformat(), disease, json.dumps(config), json.dumps(metrics), json.dumps(artifacts)),
            )
            conn.commit()

    def list_runs(self, limit: int = 20) -> List[Dict[str, Any]]:
        with self._conn() as conn:
            rows = conn.execute('SELECT run_id, created_at, disease, metrics_json, artifacts_json FROM runs ORDER BY created_at DESC LIMIT ?', (limit,)).fetchall()
        return [
            {'run_id': r[0], 'created_at': r[1], 'disease': r[2], 'metrics': json.loads(r[3]), 'artifacts': json.loads(r[4])}
            for r in rows
        ]

    def _write_event(self, run_id: str, payload: Dict[str, Any]) -> None:
        run_dir = self.root / run_id
        run_dir.mkdir(parents=True, exist_ok=True)
        with open(run_dir / 'events.jsonl', 'a', encoding='utf-8') as f:
            f.write(json.dumps({'ts': datetime.now(timezone.utc).isoformat(), **payload}) + '\n')
