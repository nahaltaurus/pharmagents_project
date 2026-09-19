# PharmAgents Production-Ready Working Project

This repository upgrades the earlier runnable demo into a more production-ready implementation inspired by the PharmAgents paper.

## Added in this version

- real **UniProt / RCSB** fetching support
- automatic **PDB download** for selected structures
- **AutoDock Vina / GNINA** integration hooks with automatic fallback to heuristic scoring
- optional **receptor preparation** from `.pdb` to `.pdbqt` when helper binaries are available
- stronger local molecule generation using **analog enumeration + BRICS recombination**
- persistent **experiment logging** with JSONL event logs and SQLite run registry
- **Streamlit web UI** with receptor upload and docking box controls
- better config management, run artifacts, and API/CLI controls

## Important reality check

This repository is runnable end to end, but heavyweight bioinformatics steps still depend on your machine and installed tools:

- live UniProt / RCSB fetch requires internet access
- real Vina / GNINA docking requires the binary installed and available on PATH
- Vina also needs receptor and ligand conversion support, usually via `obabel` or `prepare_receptor4.py`
- without those tools, the system automatically falls back to local heuristic scoring so the full pipeline still runs

So this is a **real working project** with a clean path to stronger deployment, not a fake placeholder.

## Setup

```bash
cd pharmagents_project
python -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scriptsctivate
pip install -r requirements.txt
export PYTHONPATH=src       # Windows PowerShell: $env:PYTHONPATH='src'
```

## CLI

Basic run:

```bash
python -m pharmagents.main run --disease "atopic dermatitis" --real-fetch --docking-backend auto
```

Use your own receptor and docking box:

```bash
python -m pharmagents.main run   --disease "atopic dermatitis"   --real-fetch   --docking-backend vina   --pdb-path /absolute/path/to/receptor.pdb   --center 10 10 10   --box-size 20 20 20
```

If you already have a prepared receptor:

```bash
python -m pharmagents.main run   --disease "atopic dermatitis"   --docking-backend gnina   --receptor-path /absolute/path/to/receptor.pdbqt   --center 10 10 10   --box-size 20 20 20
```

History:

```bash
python -m pharmagents.main history
```

## API

```bash
uvicorn pharmagents.api.app:app --reload
```

Endpoints:
- `GET /health`
- `POST /run`
- `GET /history`
- `GET /runs`

## Streamlit UI

```bash
streamlit run streamlit_app.py
```

The UI now supports:
- disease input
- live UniProt/RCSB fetch toggle
- backend selection
- receptor upload (`.pdb` or `.pdbqt`)
- custom pocket center
- custom docking box size

## Real docking setup

### Option 1: GNINA

1. Install `gnina` and ensure `gnina` is on PATH.
2. Provide either:
   - a prepared receptor `.pdbqt`, or
   - a `.pdb` file if your GNINA install accepts it directly.
3. Run with `--docking-backend gnina`.

### Option 2: AutoDock Vina

1. Install `vina` and ensure `vina` is on PATH.
2. Install one of these helpers:
   - `obabel`
   - `prepare_receptor4.py`
3. Provide either:
   - `--receptor-path receptor.pdbqt`, or
   - `--pdb-path receptor.pdb` so the project can auto-prepare it.
4. Pass a pocket center and box size.

If those are missing, PharmAgents will automatically fall back to heuristic docking-like scoring.

## Live structure fetching

When `--real-fetch` is enabled, the pipeline:
1. enriches local targets with real UniProt metadata
2. queries RCSB by UniProt accession
3. downloads the selected `.pdb` structure into the receptor cache
4. tries to prepare a receptor if docking tools are installed

## Experiment logging

Each run stores:
- a `runs/<run_id>/events.jsonl` log
- a `runs/<run_id>/report.json` artifact
- a global `runs/experiments.db` registry

## Architecture

```text
src/pharmagents/
  agents/
  api/
  core/
  data/
  services/
    docking.py
    experiment_logger.py
    generation.py
    uniprot_rcsb.py
  utils/
  workflows/
streamlit_app.py
```

## Suggested next upgrade after this

For a paper-quality version 3, the next steps would be:
- automatic pocket detection from co-crystal ligands
- batch docking queue
- Chemprop/DeepChem affinity prediction
- better ADMET models
- containerization with Docker + CI/CD
