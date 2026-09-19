from pharmagents.workflows.pipeline import VirtualPharmaPipeline

pipe = VirtualPharmaPipeline()
result = pipe.run("atopic dermatitis", top_k_targets=3, n_leads=8, optimization_rounds=3)
print(result["final_report"])
