from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import streamlit as st

from pharmagents.core.config import PipelineConfig
from pharmagents.workflows.pipeline import VirtualPharmaPipeline

st.set_page_config(page_title='PharmAgents', page_icon='💊', layout='wide')
st.title('💊 PharmAgents Production Demo')
st.caption('Virtual pharma pipeline with live UniProt/RCSB fetch, optional receptor preparation, docking integration hooks, experiment logging, and explainable reports.')

upload_dir = Path('.streamlit_uploads')
upload_dir.mkdir(parents=True, exist_ok=True)

with st.sidebar:
    st.header('Run settings')
    disease = st.text_input('Disease', value='atopic dermatitis')
    top_k_targets = st.slider('Top targets', 1, 5, 3)
    n_leads = st.slider('Lead count', 4, 30, 12)
    optimization_rounds = st.slider('Optimization rounds', 1, 8, 4)
    real_fetch = st.toggle('Use live UniProt/RCSB fetch', value=False)
    docking_backend = st.selectbox('Docking backend', ['auto', 'heuristic', 'vina', 'gnina'], index=0)
    st.markdown('**Optional docking inputs**')
    receptor_upload = st.file_uploader('Upload receptor file (.pdb or .pdbqt)', type=['pdb', 'pdbqt'])
    center_x = st.number_input('Pocket center X', value=0.0, step=0.5)
    center_y = st.number_input('Pocket center Y', value=0.0, step=0.5)
    center_z = st.number_input('Pocket center Z', value=0.0, step=0.5)
    size_x = st.number_input('Box size X', value=20.0, step=1.0)
    size_y = st.number_input('Box size Y', value=20.0, step=1.0)
    size_z = st.number_input('Box size Z', value=20.0, step=1.0)
    run_clicked = st.button('Run PharmAgents', type='primary', use_container_width=True)

uploaded_path = None
if receptor_upload is not None:
    uploaded_path = upload_dir / receptor_upload.name
    uploaded_path.write_bytes(receptor_upload.getbuffer())

if run_clicked:
    with st.spinner('Running virtual pharma pipeline...'):
        pipe = VirtualPharmaPipeline(config=PipelineConfig(use_real_fetch=real_fetch, docking_backend=docking_backend))
        overrides = {
            'receptor_path': str(uploaded_path) if uploaded_path and uploaded_path.suffix.lower() == '.pdbqt' else None,
            'pdb_path': str(uploaded_path) if uploaded_path and uploaded_path.suffix.lower() == '.pdb' else None,
            'pocket_center': [center_x, center_y, center_z],
            'box_size': [size_x, size_y, size_z],
        }
        result = pipe.run(disease, top_k_targets=top_k_targets, n_leads=n_leads, optimization_rounds=optimization_rounds, target_overrides=overrides)
    st.success(f"Run complete. Run ID: {result.get('run_id')}")

    col1, col2 = st.columns(2)
    with col1:
        st.subheader('Selected targets')
        st.dataframe(pd.DataFrame(result['targets']), use_container_width=True)
    with col2:
        st.subheader('Extension summary')
        st.json(result['extension_summary'])

    st.subheader('Top evaluations')
    eval_df = pd.DataFrame(result['evaluations'])
    cols = ['smiles', 'portfolio_score', 'docking_like_score', 'docking_backend', 'toxicity_risk', 'sa_like_score', 'recommended']
    st.dataframe(eval_df[cols], use_container_width=True)

    st.subheader('Final report')
    st.text(result['final_report'])

    st.subheader('Artifacts')
    st.json(result.get('artifacts', {}))

    st.download_button('Download report JSON', data=json.dumps(result, indent=2), file_name=f"pharmagents_{result.get('run_id', 'run')}.json", mime='application/json')
else:
    st.info('Set the parameters in the sidebar and click Run PharmAgents.')

st.divider()
st.subheader('Recent experiment runs')
pipe = VirtualPharmaPipeline()
runs = pipe.experiment_logger.list_runs(limit=10)
if runs:
    st.dataframe(pd.DataFrame(runs), use_container_width=True)
else:
    st.write('No runs logged yet.')
