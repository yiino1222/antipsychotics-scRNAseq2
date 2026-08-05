# Release manifest

| Path | Purpose |
| --- | --- |
| `README.md` | Reviewer-facing overview and run instructions. |
| `LICENSE` | Repository license copied for public release. |
| `Dockerfile` | Container recipe from the analysis repository. |
| `requirements.txt` | Python dependency list. |
| `code/preprocess.py` | Standardizes raw scRNA-seq inputs into AnnData. |
| `code/run_antipsychotic_pipeline.py` | Main antipsychotic-response and clozapine-selectivity workflow. |
| `code/run_camp_pattern_analysis.py` | Standalone inhibitor-cocktail cAMP-pattern analysis. |
| `code/run_human_sz_pfc_clz_response_deg.py` | Human Sz/control PFC clozapine-responsive DEG workflow. |
| `code/clz_response_expression_analysis.py` | Helper module for clozapine-response expression comparisons. |
| `code/calculation_tool.py` | Shared receptor/drug-response calculation utilities. |
| `code/rapids_scanpy_funcs.py` | RAPIDS/Scanpy helper functions used by GPU notebooks/workflows. |
| `code/utils.py` | Shared utility functions. |
| `config/preprocess_config.json` | Example preprocessing configuration. |
| `config/antipsychotic_pipeline_config.json` | Example main pipeline configuration. |
| `config/camp_analysis_config.json` | Example cAMP-pattern analysis configuration. |
| `docs/manuscript_methods.md` | Methods description matching the implemented workflow. |
