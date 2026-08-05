from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu, ttest_ind

def compute_camp_response_for_pattern(
    adata,
    gpcr_norm_df,
    gpcr_type_df,
    drug_conc,
    pattern,
    group_col="is_clz_selective",
    selected_label=True,
    ki_inhibited=0.01,
    ki_not_inhibited=10000,
):
    ki_df = gpcr_type_df.copy()
    ki_df["Ki"] = ki_not_inhibited
    for receptor_col, inhibit in pattern.items():
        if receptor_col in ki_df.index and inhibit:
            ki_df.loc[receptor_col, "Ki"] = ki_inhibited

    camp = ct.calc_camp(adata, gpcr_norm_df, ki_df, drug_conc)
    df_plot = pd.DataFrame(
        {
            "cAMP_response": camp,
            "group": np.where(adata.obs[group_col] == selected_label, "clz_selective", "non_clz_selective"),
            "leiden": adata.obs.get("leiden", pd.Series(index=adata.obs_names, dtype=str)).astype(str).values,
        },
        index=adata.obs_names,
    )
    summary = (
        df_plot.groupby("group")["cAMP_response"]
        .agg(["count", "mean", "median", "std"])
        .rename(columns={"count": "n"})
    )
    return df_plot, summary


def run(config_path: Path):
    cfg = json.loads(config_path.read_text())
    preprocess_cfg = cfg.get("preprocess", {})
    gpu_device = preprocess_cfg.get("gpu_device")
    if gpu_device is not None and str(gpu_device) != "":
        os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_device)
        print(f"[env] CUDA_VISIBLE_DEVICES={os.environ['CUDA_VISIBLE_DEVICES']}")

    import calculation_tool as ct

    out_dir = Path(cfg["output"]["out_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)

    print("[1/5] load reference parameters")
    D_R_mtx, GPCR_type_df, drug_list, _ = ct.load_parameters()

    print("[2/5] preprocess snRNAseq")
    adata, GPCR_df = ct.preprocess_adata_in_bulk(
        adata_path=cfg["input"]["adata_h5ad"],
        label=cfg["input"].get("label_filter"),
        add_markers=cfg["input"].get("add_markers"),
        is_gpu=cfg.get("preprocess", {}).get("is_gpu", False),
    )

    print("[3/5] predict drug response + clz-selective cells")
    drug_conc = cfg["drug_response"]["drug_conc"]
    adata = ct.calc_drug_response(adata, GPCR_df, GPCR_type_df, drug_list, D_R_mtx, drug_conc)
    adata, n_clz = ct.calc_clz_selective_cell(adata, drug_list, cfg["drug_response"].get("selectivity_threshold", 1.5))
    pd.Series(adata.obs["is_clz_selective"]).value_counts().to_csv(out_dir / "is_clz_selective_counts.csv")

    print("[4/5] optional cAMP pattern analysis")
    if cfg.get("camp_pattern", {}).get("enabled", True):
        GPCR_adata_norm_df = pd.DataFrame(index=adata.obs_names)
        gpcr_cols = [c for c in GPCR_df.columns]
        # rebuild normalized GPCR expression table (columns without _raw suffix)
        import scanpy as sc
        gpcr_adata = ct.anndata.AnnData(X=GPCR_df)
        gpcr_norm = sc.pp.normalize_total(gpcr_adata, target_sum=1e4, inplace=False)["X"]
        GPCR_adata_norm_df = pd.DataFrame(gpcr_norm)
        GPCR_adata_norm_df.columns = [str(c)[:-4] for c in gpcr_cols]

        df_plot, summary = compute_camp_response_for_pattern(
            adata=adata,
            gpcr_norm_df=GPCR_adata_norm_df,
            gpcr_type_df=GPCR_type_df,
            drug_conc=drug_conc,
            pattern=cfg["camp_pattern"]["pattern"],
            group_col=cfg["camp_pattern"].get("group_col", "is_clz_selective"),
            selected_label=cfg["camp_pattern"].get("selected_label", True),
            ki_inhibited=cfg["camp_pattern"].get("ki_inhibited", 0.01),
            ki_not_inhibited=cfg["camp_pattern"].get("ki_not_inhibited", 10000),
        )
        df_cmp = df_plot.dropna(subset=["cAMP_response", "group"]).copy()
        clz_vals = df_cmp.loc[df_cmp["group"] == "clz_selective", "cAMP_response"].values
        non_vals = df_cmp.loc[df_cmp["group"] == "non_clz_selective", "cAMP_response"].values
        stats_df = pd.DataFrame({
            "test": ["Welch_ttest", "Mann_Whitney_U"],
            "statistic": [ttest_ind(clz_vals, non_vals, equal_var=False, nan_policy="omit").statistic,
                          mannwhitneyu(clz_vals, non_vals, alternative="two-sided").statistic],
            "p_value": [ttest_ind(clz_vals, non_vals, equal_var=False, nan_policy="omit").pvalue,
                        mannwhitneyu(clz_vals, non_vals, alternative="two-sided").pvalue],
        })
        df_plot.to_csv(out_dir / "camp_df_plot.csv")
        summary.to_csv(out_dir / "camp_summary.csv")
        stats_df.to_csv(out_dir / "camp_stats.csv", index=False)

    print("[5/5] save processed outputs")
    if cfg["output"].get("save_processed_h5ad", True):
        adata.write(out_dir / cfg["output"].get("processed_h5ad_name", "adata_processed.h5ad"))

    print(f"done: n_clz_selective={n_clz}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="antipsychotic_pipeline_config.json")
    args = parser.parse_args()
    run(Path(args.config))
