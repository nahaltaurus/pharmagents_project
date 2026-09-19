from __future__ import annotations

import json
import sqlite3
from pathlib import Path
from typing import Any, Dict, List

from pharmagents.core.config import DB_PATH


class ExperienceDB:
    def __init__(self, db_path: Path | None = None):
        self.db_path = db_path or DB_PATH
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self._init_db()

    def _conn(self):
        return sqlite3.connect(self.db_path)

    def _init_db(self) -> None:
        with self._conn() as conn:
            conn.execute(
                """
                CREATE TABLE IF NOT EXISTS experiences (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    disease TEXT NOT NULL,
                    target_name TEXT NOT NULL,
                    input_smiles TEXT NOT NULL,
                    output_smiles TEXT NOT NULL,
                    score REAL NOT NULL,
                    notes TEXT NOT NULL
                )
                """
            )
            conn.commit()

    def add(self, disease: str, target_name: str, input_smiles: str, output_smiles: str, score: float, notes: Dict[str, Any]) -> None:
        with self._conn() as conn:
            conn.execute(
                "INSERT INTO experiences(disease,target_name,input_smiles,output_smiles,score,notes) VALUES(?,?,?,?,?,?)",
                (disease, target_name, input_smiles, output_smiles, float(score), json.dumps(notes)),
            )
            conn.commit()

    def fetch_for_target(self, target_name: str, limit: int = 20) -> List[Dict[str, Any]]:
        with self._conn() as conn:
            rows = conn.execute(
                "SELECT disease,target_name,input_smiles,output_smiles,score,notes FROM experiences WHERE target_name=? ORDER BY score DESC LIMIT ?",
                (target_name, limit),
            ).fetchall()
        out = []
        for row in rows:
            out.append({
                "disease": row[0],
                "target_name": row[1],
                "input_smiles": row[2],
                "output_smiles": row[3],
                "score": row[4],
                "notes": json.loads(row[5]),
            })
        return out

    def history(self, limit: int = 50) -> List[Dict[str, Any]]:
        with self._conn() as conn:
            rows = conn.execute(
                "SELECT disease,target_name,input_smiles,output_smiles,score,notes FROM experiences ORDER BY id DESC LIMIT ?",
                (limit,),
            ).fetchall()
        return [
            {
                "disease": r[0], "target_name": r[1], "input_smiles": r[2], "output_smiles": r[3], "score": r[4], "notes": json.loads(r[5])
            }
            for r in rows
        ]
