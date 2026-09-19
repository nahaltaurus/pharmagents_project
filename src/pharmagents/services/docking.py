from __future__ import annotations

import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Dict

from pharmagents.core.config import PipelineConfig
from pharmagents.core.schemas import TargetCandidate
from pharmagents.utils.chem import basic_descriptors, make_3d_mol_file
from pharmagents.utils.logging import get_logger

logger = get_logger(__name__)


class DockingService:
    def __init__(self, backend: str = 'auto', config: PipelineConfig | None = None):
        self.config = config or PipelineConfig(docking_backend=backend)
        self.backend = backend
        self.vina_bin = shutil.which('vina')
        self.gnina_bin = shutil.which('gnina')
        self.obabel_bin = shutil.which('obabel')
        self.prepare_receptor_bin = shutil.which('prepare_receptor4.py') or shutil.which('mk_prepare_receptor.py')
        self.config.receptor_cache_dir.mkdir(parents=True, exist_ok=True)

    def available_backend(self) -> str:
        if self.backend == 'gnina' and self.gnina_bin:
            return 'gnina'
        if self.backend == 'vina' and self.vina_bin:
            return 'vina'
        if self.backend == 'auto':
            if self.gnina_bin:
                return 'gnina'
            if self.vina_bin:
                return 'vina'
        return 'heuristic'

    def ensure_target_receptor(self, target: TargetCandidate) -> TargetCandidate:
        metadata = dict(target.metadata)
        receptor_path = metadata.get('receptor_path')
        if receptor_path and Path(receptor_path).exists():
            return target
        pdb_path = metadata.get('pdb_path')
        backend = self.available_backend()
        if not pdb_path or not Path(pdb_path).exists():
            return target
        if backend == 'heuristic':
            metadata['receptor_path'] = str(pdb_path)
            target.metadata = metadata
            return target
        prepared = self._prepare_receptor_file(Path(pdb_path), backend)
        if prepared:
            metadata['receptor_path'] = prepared
            target.metadata = metadata
        return target

    def score(self, smiles: str, target: TargetCandidate) -> Dict[str, Any]:
        target = self.ensure_target_receptor(target)
        backend = self.available_backend()
        if backend == 'heuristic':
            return self._heuristic_score(smiles, target)
        try:
            return self._score_with_binary(smiles, target, backend)
        except Exception as exc:
            logger.warning('Docking with %s failed, falling back to heuristic: %s', backend, exc)
            return self._heuristic_score(smiles, target)

    def _heuristic_score(self, smiles: str, target: TargetCandidate) -> Dict[str, Any]:
        d = basic_descriptors(smiles)
        target_bias = {'JAK1': 1.2, 'JAK2': 1.0, 'JAK3': 1.1}.get(target.target_name, 0.8)
        score = float(-(6.0 + 2.0 * d.get('qed', 0.0) + 0.35 * d.get('fsp3', 0.0) + 0.2 * target_bias))
        return {'backend': 'heuristic', 'affinity': score, 'details': {'reason': 'Local heuristic estimate'}}

    def _prepare_receptor_file(self, pdb_path: Path, backend: str) -> str | None:
        if backend == 'gnina':
            return str(pdb_path)
        out_path = self.config.receptor_cache_dir / f'{pdb_path.stem}.pdbqt'
        if out_path.exists():
            return str(out_path)
        if self.prepare_receptor_bin:
            cmd = [self.prepare_receptor_bin, '-r', str(pdb_path), '-o', str(out_path)]
        elif self.obabel_bin:
            cmd = [self.obabel_bin, str(pdb_path), '-O', str(out_path)]
        else:
            logger.warning('No receptor preparation binary found for %s', pdb_path)
            return None
        proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if proc.returncode != 0 or not out_path.exists():
            logger.warning('Receptor preparation failed: %s %s', proc.stdout[-500:], proc.stderr[-500:])
            return None
        return str(out_path)

    def _prepare_ligand_file(self, smiles: str, backend: str, workdir: Path) -> str:
        lig_sdf = workdir / 'ligand.sdf'
        if not make_3d_mol_file(smiles, str(lig_sdf)):
            raise RuntimeError('Failed to generate 3D ligand')
        if backend == 'gnina':
            return str(lig_sdf)
        if not self.obabel_bin:
            raise RuntimeError('Open Babel (obabel) is required to convert ligand to PDBQT for Vina')
        lig_pdbqt = workdir / 'ligand.pdbqt'
        proc = subprocess.run([self.obabel_bin, str(lig_sdf), '-O', str(lig_pdbqt)], capture_output=True, text=True, check=False)
        if proc.returncode != 0 or not lig_pdbqt.exists():
            raise RuntimeError(proc.stderr.strip() or proc.stdout.strip() or 'Ligand conversion to PDBQT failed')
        return str(lig_pdbqt)

    def _score_with_binary(self, smiles: str, target: TargetCandidate, backend: str) -> Dict[str, Any]:
        receptor = target.metadata.get('receptor_path')
        if not receptor or not Path(receptor).exists():
            raise RuntimeError('No prepared receptor path available for docking')
        with tempfile.TemporaryDirectory(prefix='pharmagents_docking_') as tmp:
            tmp_path = Path(tmp)
            ligand = self._prepare_ligand_file(smiles, backend, tmp_path)
            out_path = tmp_path / ('out.sdf' if backend == 'gnina' else 'out.pdbqt')
            center = target.pocket_center or [0.0, 0.0, 0.0]
            size = target.metadata.get('box_size', [20, 20, 20])
            cmd = [
                self.gnina_bin if backend == 'gnina' else self.vina_bin,
                '--receptor', str(receptor),
                '--ligand', str(ligand),
                '--center_x', str(center[0]), '--center_y', str(center[1]), '--center_z', str(center[2]),
                '--size_x', str(size[0]), '--size_y', str(size[1]), '--size_z', str(size[2]),
                '--out', str(out_path),
            ]
            if backend == 'gnina':
                cmd.extend(['--score_only'])
            proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
            if proc.returncode != 0:
                raise RuntimeError(proc.stderr.strip() or proc.stdout.strip())
            affinity = self._parse_affinity(proc.stdout + '\n' + proc.stderr)
            return {
                'backend': backend,
                'affinity': affinity,
                'details': {
                    'stdout': proc.stdout[-4000:],
                    'stderr': proc.stderr[-4000:],
                    'receptor_path': str(receptor),
                    'ligand_path': str(ligand),
                },
            }

    @staticmethod
    def _parse_affinity(text: str) -> float:
        patterns = [
            r'affinity[^\-\d]*(-?\d+\.\d+)',
            r'cnnscore[^\-\d]*(-?\d+\.\d+)',
            r'\b1\s+(-?\d+\.\d+)\s+\d',
        ]
        lower = text.lower()
        for pattern in patterns:
            match = re.search(pattern, lower)
            if match:
                return float(match.group(1))
        for line in text.splitlines():
            parts = [p for p in line.replace(':', ' ').split() if re.fullmatch(r'-?\d+(?:\.\d+)?', p)]
            if parts:
                return float(parts[0])
        raise RuntimeError('Could not parse docking score')
