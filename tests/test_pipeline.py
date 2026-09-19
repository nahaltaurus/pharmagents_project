from pharmagents.core.config import PipelineConfig
from pharmagents.workflows.pipeline import VirtualPharmaPipeline


def test_pipeline_runs():
    pipe = VirtualPharmaPipeline(config=PipelineConfig(use_real_fetch=False, docking_backend='heuristic'))
    result = pipe.run('atopic dermatitis', top_k_targets=2, n_leads=6, optimization_rounds=2)
    assert result['disease'] == 'atopic dermatitis'
    assert len(result['targets']) >= 1
    assert len(result['initial_leads']) >= 1
    assert len(result['optimized_leads']) >= 1
    assert 'run_id' in result and result['run_id']
    assert 'artifacts' in result
