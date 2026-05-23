from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu, ttest_ind

import calculation_tool as ct


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
    out_dir = Path(cfg["output"]["out_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)

    print("[1/6] load inputs")
    adata = anndata.read_h5ad(cfg["input"]["adata_h5ad"])
    gpcr_norm_df = pd.read_csv(cfg["input"]["gpcr_norm_csv"], index_col=0)
    gpcr_type_df = pd.read_csv(cfg["input"]["gpcr_type_csv"], index_col=0)

    a = cfg["analysis"]
    pattern = a["pattern"]
    print("[2/6] validate pattern")
    missing = [k for k in pattern if k not in gpcr_type_df.index]
    if missing:
        print(f"[WARN] missing receptors in gpcr_type_df index: {missing}")

    print("[3/6] compute cAMP")
    df_plot, summary = compute_camp_response_for_pattern(
        adata=adata,
        gpcr_norm_df=gpcr_norm_df,
        gpcr_type_df=gpcr_type_df,
        drug_conc=a["drug_conc"],
        pattern=pattern,
        group_col=a.get("group_col", "is_clz_selective"),
        selected_label=a.get("selected_label", True),
        ki_inhibited=a.get("ki_inhibited", 0.01),
        ki_not_inhibited=a.get("ki_not_inhibited", 10000),
    )

    print("[4/6] statistical tests")
    df_cmp = df_plot.dropna(subset=["cAMP_response", "group"]).copy()
    clz_vals = df_cmp.loc[df_cmp["group"] == "clz_selective", "cAMP_response"].values
    non_vals = df_cmp.loc[df_cmp["group"] == "non_clz_selective", "cAMP_response"].values

    welch_t = ttest_ind(clz_vals, non_vals, equal_var=False, nan_policy="omit")
    mann_u = mannwhitneyu(clz_vals, non_vals, alternative="two-sided")
    stats_df = pd.DataFrame(
        {
            "test": ["Welch_ttest", "Mann_Whitney_U"],
            "statistic": [welch_t.statistic, mann_u.statistic],
            "p_value": [welch_t.pvalue, mann_u.pvalue],
            "n_clz": [len(clz_vals), len(clz_vals)],
            "n_non": [len(non_vals), len(non_vals)],
            "mean_clz": [np.mean(clz_vals), np.mean(clz_vals)],
            "mean_non": [np.mean(non_vals), np.mean(non_vals)],
            "median_clz": [np.median(clz_vals), np.median(clz_vals)],
            "median_non": [np.median(non_vals), np.median(non_vals)],
        }
    )

    print("[5/6] write intermediate outputs")
    df_plot.to_csv(out_dir / "camp_df_plot.csv")
    summary.to_csv(out_dir / "camp_summary.csv")
    stats_df.to_csv(out_dir / "camp_stats.csv", index=False)

    print("[6/6] save plots")
    plot_df = df_cmp.copy()
    plot_df["group"] = pd.Categorical(plot_df["group"], ["non_clz_selective", "clz_selective"])

    plt.figure(figsize=(6, 5))
    sns.violinplot(data=plot_df, x="group", y="cAMP_response", inner=None, cut=0, color="lightgray")
    sns.boxplot(data=plot_df, x="group", y="cAMP_response", width=0.35, showfliers=False, boxprops={"facecolor": "white"})
    sns.swarmplot(data=plot_df, x="group", y="cAMP_response", size=2.0, alpha=0.5, color="black")
    plt.title("Predicted cAMP under 5-HT1A/M2/H3 inhibition")
    plt.tight_layout()
    plt.savefig(out_dir / "camp_group_comparison.png", dpi=200)
    plt.close()

    if cfg.get("plot", {}).get("make_cluster_plot", True) and cfg.get("plot", {}).get("cluster_col", "leiden") in plot_df.columns:
        plt.figure(figsize=(12, 4))
        sns.boxplot(data=plot_df, x=cfg.get("plot", {}).get("cluster_col", "leiden"), y="cAMP_response", hue="group", showfliers=False)
        plt.tight_layout()
        plt.savefig(out_dir / "camp_cluster_comparison.png", dpi=200)
        plt.close()

    print("done")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="camp_analysis_config.json")
    args = parser.parse_args()
    run(Path(args.config))
