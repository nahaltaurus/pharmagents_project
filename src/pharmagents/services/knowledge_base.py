from __future__ import annotations

from typing import Any, Dict, List
from pharmagents.core.config import DATA_DIR
from pharmagents.utils.io import load_json


class KnowledgeBase:
    def __init__(self):
        self.targets = load_json(DATA_DIR / 'disease_targets.json')
        self.library = load_json(DATA_DIR / 'molecule_library.json')
        self.tox = load_json(DATA_DIR / 'tox_reference.json')

    def disease_candidates(self, disease: str) -> List[Dict[str, Any]]:
        disease_l = disease.lower()
        all_rows = self.targets['diseases']
        scored = []
        for row in all_rows:
            name = row['name'].lower()
            aliases = [a.lower() for a in row.get('aliases', [])]
            tags = [t.lower() for t in row.get('tags', [])]
            score = 0.0
            if disease_l == name or disease_l in aliases:
                score += 1.0
            score += 0.2 * sum(tok in disease_l for tok in tags)
            score += 0.1 * sum(tok in name for tok in disease_l.split())
            scored.append((score, row))
        scored.sort(key=lambda x: x[0], reverse=True)
        return [r for s, r in scored[:8] if s > 0]

    def get_molecule_library(self) -> List[Dict[str, Any]]:
        return self.library['molecules']

    def tox_reference(self) -> List[Dict[str, Any]]:
        return self.tox['records']
