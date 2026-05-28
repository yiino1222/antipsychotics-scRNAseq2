from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc


def _read_input(path: str, file_type: str):
    p = Path(path)
    if file_type == "h5ad":
        return ad.read_h5ad(p)
    if file_type == "10x_mtx":
        return sc.read_10x_mtx(p, var_names="gene_symbols", cache=False)
    if file_type == "10x_h5":
        return sc.read_10x_h5(p)
    if file_type == "csv":
        df = pd.read_csv(p, index_col=0)
        return ad.AnnData(X=df.values, obs=pd.DataFrame(index=df.index), var=pd.DataFrame(index=df.columns))
    raise ValueError(f"unsupported file_type: {file_type}")


def _harmonize(adata, species: str | None = None):
    adata.var_names = pd.Index([str(x).upper() for x in adata.var_names])
    adata.var_names_make_unique()
    if species is not None:
        adata.obs["species"] = species
    return adata


def preprocess_from_config(cfg_path: Path):
    cfg = json.loads(cfg_path.read_text())
    out_dir = Path(cfg["output"]["out_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)

    matrices = []
    for i, ds in enumerate(cfg["datasets"], start=1):
        print(f"[{i}/{len(cfg['datasets'])}] loading {ds['name']}")
        adata = _read_input(ds["path"], ds["file_type"])
        adata = _harmonize(adata, ds.get("species"))
        adata.obs["dataset"] = ds["name"]
        adata.obs["region"] = ds.get("region", "unknown")
        matrices.append(adata)

    if len(matrices) == 1:
        merged = matrices[0]
    else:
        merged = ad.concat(matrices, axis=0, join="outer", label="batch", fill_value=0)

    merged.write_h5ad(out_dir / cfg["output"].get("matrix_name", "standardized_matrix.h5ad"))
    merged.obs.to_csv(out_dir / "standardized_obs.csv")
    pd.DataFrame(index=merged.var_names).to_csv(out_dir / "standardized_var.csv")
    print(f"done: cells={merged.n_obs}, genes={merged.n_vars}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="preprocess_config.json")
    args = parser.parse_args()
    preprocess_from_config(Path(args.config))
