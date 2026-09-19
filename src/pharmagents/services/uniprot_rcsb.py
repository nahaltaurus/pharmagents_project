from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

import requests

from pharmagents.core.config import PipelineConfig
from pharmagents.utils.logging import get_logger

logger = get_logger(__name__)


class UniProtRCSBClient:
    def __init__(self, config: PipelineConfig | None = None):
        self.config = config or PipelineConfig()
        self.session = requests.Session()
        self.session.headers.update({'User-Agent': 'PharmAgents/0.3'})
        self.config.cache_dir.mkdir(parents=True, exist_ok=True)
        self.config.receptor_cache_dir.mkdir(parents=True, exist_ok=True)

    def _get_json(self, url: str, params: Optional[Dict[str, Any]] = None) -> Optional[Dict[str, Any]]:
        try:
            resp = self.session.get(url, params=params, timeout=self.config.api_timeout_sec)
            resp.raise_for_status()
            return resp.json()
        except Exception as exc:
            logger.warning('GET failed for %s: %s', url, exc)
            return None

    def search_uniprot(self, query: str, size: int = 5) -> List[Dict[str, Any]]:
        data = self._get_json(
            'https://rest.uniprot.org/uniprotkb/search',
            params={'query': query, 'format': 'json', 'size': size},
        )
        results = []
        for row in (data or {}).get('results', []):
            accession = row.get('primaryAccession')
            protein_desc = row.get('proteinDescription', {})
            recommended = protein_desc.get('recommendedName', {}).get('fullName', {}).get('value', '')
            genes = row.get('genes', [])
            gene = genes[0].get('geneName', {}).get('value', '') if genes else ''
            results.append({'uniprot_id': accession, 'protein_name': recommended or gene or accession, 'gene': gene})
        return results

    def fetch_uniprot_entry(self, accession: str) -> Dict[str, Any]:
        data = self._get_json(f'https://rest.uniprot.org/uniprotkb/{accession}.json') or {}
        return {
            'uniprot_id': accession,
            'protein_name': data.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', accession),
            'organism': data.get('organism', {}).get('scientificName', ''),
            'comments': data.get('comments', []),
        }

    def search_pdb_for_uniprot(self, accession: str, rows: int = 5) -> List[Dict[str, Any]]:
        payload = {
            'query': {
                'type': 'terminal',
                'service': 'text',
                'parameters': {
                    'attribute': 'rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession',
                    'operator': 'exact_match',
                    'value': accession,
                },
            },
            'request_options': {'pager': {'start': 0, 'rows': rows}, 'return_all_hits': False},
            'return_type': 'entry',
        }
        try:
            resp = self.session.post('https://search.rcsb.org/rcsbsearch/v2/query', json=payload, timeout=self.config.api_timeout_sec)
            resp.raise_for_status()
            data = resp.json()
        except Exception as exc:
            logger.warning('RCSB search failed for %s: %s', accession, exc)
            return []

        entries = []
        for item in data.get('result_set', []):
            pdb_id = item.get('identifier')
            if not pdb_id:
                continue
            meta = self.fetch_pdb_summary(pdb_id)
            ligand = self.fetch_best_ligand_info(pdb_id)
            entries.append({'pdb_id': pdb_id, **meta, 'ligand': ligand})
        return entries

    def fetch_pdb_summary(self, pdb_id: str) -> Dict[str, Any]:
        data = self._get_json(f'https://data.rcsb.org/rest/v1/core/entry/{pdb_id}') or {}
        struct = data.get('struct', {})
        title = struct.get('title', pdb_id)
        return {
            'title': title,
            'experimental_method': ', '.join(data.get('exptl', [{}])[0].get('method', '').split(';')) if data.get('exptl') else '',
            'resolution': data.get('rcsb_entry_info', {}).get('resolution_combined', [None])[0],
        }

    def fetch_best_ligand_info(self, pdb_id: str) -> Dict[str, Any]:
        data = self._get_json(f'https://data.rcsb.org/rest/v1/core/nonpolymer_entity/{pdb_id}/1') or {}
        chem = data.get('chem_comp', {})
        return {
            'id': chem.get('id', ''),
            'name': chem.get('name', ''),
            'formula': chem.get('formula', ''),
        }

    def download_pdb_structure(self, pdb_id: str, overwrite: bool = False) -> Optional[str]:
        out_path = self.config.receptor_cache_dir / f'{pdb_id}.pdb'
        if out_path.exists() and not overwrite:
            return str(out_path)
        url = f'https://files.rcsb.org/download/{pdb_id}.pdb'
        try:
            resp = self.session.get(url, timeout=self.config.api_timeout_sec)
            resp.raise_for_status()
            out_path.parent.mkdir(parents=True, exist_ok=True)
            out_path.write_bytes(resp.content)
            return str(out_path)
        except Exception as exc:
            logger.warning('Could not download PDB %s: %s', pdb_id, exc)
            return None
