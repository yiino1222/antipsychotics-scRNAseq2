# Manuscript code release package

This folder contains a self-contained code package prepared for public repository upload so that reviewers can inspect and rerun the analyses used to generate manuscript data for the bioRxiv preprint:

- DOI/preprint URL: https://www.biorxiv.org/content/10.64898/2026.06.14.732113v1.abstract

The package intentionally includes code, configuration templates, environment files, and methods documentation only. Raw or processed large single-cell data files are not bundled and should be placed locally by the reviewer according to the paths in `config/*.json`.

## Folder layout

```text
public_repository_upload/
├── README.md
├── LICENSE
├── Dockerfile
├── requirements.txt
├── code/
│   ├── preprocess.py
│   ├── run_antipsychotic_pipeline.py
│   ├── run_camp_pattern_analysis.py
│   ├── run_human_sz_pfc_clz_response_deg.py
│   ├── clz_response_expression_analysis.py
│   ├── calculation_tool.py
│   ├── rapids_scanpy_funcs.py
│   └── utils.py
├── config/
│   ├── preprocess_config.json
│   ├── antipsychotic_pipeline_config.json
│   └── camp_analysis_config.json
└── docs/
    └── manuscript_methods.md
```

## What each entry point does

1. `code/preprocess.py`
   - Converts one or more raw count matrices into a standardized AnnData file.
   - Harmonizes gene symbols to uppercase and records dataset/species/region metadata.

2. `code/run_antipsychotic_pipeline.py`
   - Runs the main receptor-expression and antipsychotic-response workflow.
   - Calculates modeled cAMP/Ca response scores, clozapine-selective cell labels, and optional inhibitor-cocktail cAMP summaries.

3. `code/run_camp_pattern_analysis.py`
   - Runs the inhibitor-cocktail cAMP-pattern analysis as a standalone workflow.

4. `code/run_human_sz_pfc_clz_response_deg.py`
   - Reproduces the human schizophrenia/control PFC clozapine-responsive DEG workflow.
   - Writes CSV outputs and publication-editable PDF figures.

5. `code/clz_response_expression_analysis.py`
   - Notebook-friendly helper functions for clozapine-responsive expression comparisons and receptor-focused plots.

6. `code/calculation_tool.py`, `code/rapids_scanpy_funcs.py`, and `code/utils.py`
   - Shared functions used by the workflows above.

## Minimal reviewer workflow

From inside `public_repository_upload/`, install dependencies, edit config paths to point to local data, and run:

```bash
python code/preprocess.py --config config/preprocess_config.json
python code/run_antipsychotic_pipeline.py --config config/antipsychotic_pipeline_config.json
```

Optional analyses:

```bash
python code/run_camp_pattern_analysis.py --config config/camp_analysis_config.json
python code/run_human_sz_pfc_clz_response_deg.py --help
```

## Configuration notes

- `config/preprocess_config.json` is a template for raw-data-to-standardized-AnnData conversion.
- `config/antipsychotic_pipeline_config.json` assumes the standardized output from preprocessing is available at `outputs/preprocess/standardized_matrix.h5ad`.
- `config/camp_analysis_config.json` is a standalone cAMP-pattern template and requires precomputed GPCR-normalized and receptor-type tables.
- GPU usage is optional. Set `preprocess.is_gpu` and `preprocess.gpu_device` in `config/antipsychotic_pipeline_config.json` if a CUDA/RAPIDS environment is available.

## Data policy

No controlled, private, or large binary data are included in this upload folder. Reviewers should obtain public datasets from the manuscript's data-availability statement or place authorized local data at the configured paths.
