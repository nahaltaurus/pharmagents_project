from __future__ import annotations

import math
import random
from itertools import islice
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, BRICS, Crippen, Descriptors, Lipinski, QED, rdMolDescriptors
from rdkit.DataStructs import TanimotoSimilarity

RDLogger.DisableLog('rdApp.*')


COMMON_DECORATIONS = ['F', 'Cl', 'Br', 'C#N', 'OC', 'NC', 'N(C)C', 'C(=O)N']


def mol_from_smiles(smiles: str):
    return Chem.MolFromSmiles(smiles)


def canonicalize(smiles: str) -> str:
    mol = mol_from_smiles(smiles)
    if mol is None:
        return smiles
    return Chem.MolToSmiles(mol)


def sanitize_smiles(smiles: str) -> Optional[str]:
    mol = mol_from_smiles(smiles)
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(mol)
        return Chem.MolToSmiles(mol)
    except Exception:
        return None


def basic_descriptors(smiles: str) -> Dict[str, float]:
    mol = mol_from_smiles(smiles)
    if mol is None:
        return {'valid': 0.0}
    return {
        'valid': 1.0,
        'mw': float(Descriptors.MolWt(mol)),
        'logp': float(Crippen.MolLogP(mol)),
        'hbd': float(Lipinski.NumHDonors(mol)),
        'hba': float(Lipinski.NumHAcceptors(mol)),
        'rot_bonds': float(Lipinski.NumRotatableBonds(mol)),
        'tpsa': float(rdMolDescriptors.CalcTPSA(mol)),
        'rings': float(rdMolDescriptors.CalcNumRings(mol)),
        'qed': float(QED.qed(mol)),
        'heavy_atoms': float(mol.GetNumHeavyAtoms()),
        'fsp3': float(rdMolDescriptors.CalcFractionCSP3(mol)),
    }


def sa_like_score(smiles: str) -> float:
    d = basic_descriptors(smiles)
    if d.get('valid', 0.0) == 0.0:
        return 0.0
    penalties = 0.0
    penalties += max(0.0, (d['mw'] - 550) / 550)
    penalties += max(0.0, (d['rot_bonds'] - 8) / 8)
    penalties += max(0.0, (d['rings'] - 4) / 4)
    penalties += max(0.0, (d['tpsa'] - 140) / 140)
    penalties += max(0.0, (0.1 - d['fsp3']) / 0.1) * 0.15
    return max(0.0, 1.0 - 0.35 * penalties)


def fingerprint(smiles: str):
    mol = mol_from_smiles(smiles)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)


def tanimoto(smiles_a: str, smiles_b: str) -> float:
    fa, fb = fingerprint(smiles_a), fingerprint(smiles_b)
    if fa is None or fb is None:
        return 0.0
    return float(TanimotoSimilarity(fa, fb))


def diversity_select(smiles_list: List[str], scores: Dict[str, float], k: int, max_similarity: float = 0.78) -> List[str]:
    ordered = sorted(smiles_list, key=lambda s: scores.get(s, 0.0), reverse=True)
    selected: List[str] = []
    for smi in ordered:
        if not selected or all(tanimoto(smi, x) < max_similarity for x in selected):
            selected.append(smi)
        if len(selected) >= k:
            break
    return selected


def nearest_neighbors(query: str, library: Iterable[str], top_k: int = 5) -> List[Tuple[str, float]]:
    vals = [(s, tanimoto(query, s)) for s in library]
    vals.sort(key=lambda x: x[1], reverse=True)
    return vals[:top_k]


def score_druglikeness(smiles: str) -> float:
    d = basic_descriptors(smiles)
    if d.get('valid', 0.0) == 0.0:
        return 0.0
    lipinski_bonus = 1.0 if d['mw'] < 550 and d['hba'] <= 10 and d['hbd'] <= 5 and d['logp'] < 5.5 else 0.6
    return 0.55 * d['qed'] + 0.25 * sa_like_score(smiles) + 0.20 * lipinski_bonus


def enumerate_analogs(smiles: str, limit: int = 24) -> List[str]:
    base = sanitize_smiles(smiles)
    if not base:
        return []
    mol = mol_from_smiles(base)
    candidates = {base}
    replacements = [
        ('Cl', 'F'), ('F', 'Cl'), ('Br', 'Cl'), ('OC', 'O'), ('N', 'NC'), ('c1ccccc1', 'c1ccncc1'),
        ('C(=O)N', 'C(=O)NC'), ('C', 'CC'),
    ]
    for a, b in replacements:
        if a in base:
            s = sanitize_smiles(base.replace(a, b, 1))
            if s:
                candidates.add(s)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetTotalNumHs() > 0:
            rw = Chem.RWMol(mol)
            new_idx = rw.AddAtom(Chem.Atom(9))
            rw.AddBond(atom.GetIdx(), new_idx, Chem.BondType.SINGLE)
            s = sanitize_smiles(Chem.MolToSmiles(rw))
            if s:
                candidates.add(s)
            if len(candidates) >= limit:
                break
    out = [s for s in candidates if basic_descriptors(s).get('valid', 0.0) == 1.0]
    out.sort(key=score_druglikeness, reverse=True)
    return out[:limit]


def brics_recombine(smiles_list: Sequence[str], limit: int = 48) -> List[str]:
    mols = [mol_from_smiles(s) for s in smiles_list if mol_from_smiles(s) is not None]
    frags = []
    for mol in mols[:10]:
        try:
            frags.extend(BRICS.BRICSDecompose(mol))
        except Exception:
            continue
    if not frags:
        return []
    try:
        built = list(islice(BRICS.BRICSBuild([Chem.MolFromSmiles(f) for f in set(frags) if Chem.MolFromSmiles(f) is not None]), limit * 3))
    except Exception:
        return []
    out = []
    for mol in built:
        if mol is None:
            continue
        smi = sanitize_smiles(Chem.MolToSmiles(mol))
        if smi:
            out.append(smi)
    uniq = list(dict.fromkeys(out))
    uniq.sort(key=score_druglikeness, reverse=True)
    return uniq[:limit]


def make_3d_mol_file(smiles: str, out_path: str) -> bool:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    mol = Chem.AddHs(mol)
    try:
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        AllChem.MMFFOptimizeMolecule(mol)
        writer = Chem.SDWriter(out_path)
        writer.write(mol)
        writer.close()
        return True
    except Exception:
        return False
