"""Predict clozapine-responsive human PFC cells and run Sz-vs-control DEG.

This script follows the workflow in ``human_Sz_PFC_20260624.ipynb``:
preprocess separate schizophrenia (Sz) and healthy-control AnnData files with
``calculation_tool.preprocess_adata_in_bulk``, calculate GPCR-based drug
responses, label clozapine-selective cells, and compare gene expression between
Sz and healthy-control clozapine-responsive cells.

Outputs are written as CSV tables plus Illustrator-editable PDF figures. The PDF
backend keeps text editable by embedding TrueType fonts (fonttype 42).
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import sparse

import calculation_tool as ct


DEFAULT_SZ_H5AD = "/data/human_Sz_PFC_each/merged_SZ_adata.h5ad"
DEFAULT_CONTROL_H5AD = "/data/human_Sz_PFC_each/merged_control_adata.h5ad"
DEFAULT_OUTPUT_DIR = "/data/human_Sz_PFC_each/clz_response_deg"


def configure_plotting() -> None:
    """Configure matplotlib PDFs so labels remain editable in Illustrator."""
    matplotlib.rcParams["pdf.fonttype"] = 42
    matplotlib.rcParams["ps.fonttype"] = 42
    matplotlib.rcParams["font.family"] = "DejaVu Sans"
    matplotlib.rcParams["svg.fonttype"] = "none"
    sns.set_theme(style="whitegrid", context="talk")


def preprocess_and_predict(
    adata_path: str | Path,
    group_label: str,
    drug_conc: float,
    selectivity_threshold: float,
    is_gpu: bool,
    label_filter: str | None = None,
    add_markers: list[str] | None = None,
) -> tuple[anndata.AnnData, pd.DataFrame]:
    """Preprocess one cohort and add clozapine-response annotations."""
    d_r_mtx, gpcr_type_df, drug_list, _ = ct.load_parameters()
    adata, gpcr_df = ct.preprocess_adata_in_bulk(
        adata_path=str(adata_path),
        label=label_filter,
        add_markers=add_markers,
        is_gpu=is_gpu,
    )
    adata = ct.calc_drug_response(adata, gpcr_df, gpcr_type_df, drug_list, d_r_mtx, drug_conc)
    adata, _ = ct.calc_clz_selective_cell(adata, drug_list, selectivity_threshold)
    adata.obs["diagnosis_group"] = group_label
    adata.obs["source_h5ad"] = str(adata_path)
    return adata, gpcr_df


def make_unique_obs_names(adata: anndata.AnnData, prefix: str) -> anndata.AnnData:
    """Prefix cell IDs before concatenating cohorts."""
    adata = adata.copy()
    adata.obs_names = pd.Index([f"{prefix}_{cell}" for cell in adata.obs_names.astype(str)])
    return adata


def concatenate_clz_responsive(sz_adata: anndata.AnnData, control_adata: anndata.AnnData) -> anndata.AnnData:
    """Return a combined AnnData containing only predicted clozapine-responsive cells."""
    sz_clz = sz_adata[sz_adata.obs["is_clz_selective"].astype(str).isin(["True", "1", "true"])].copy()
    control_clz = control_adata[control_adata.obs["is_clz_selective"].astype(str).isin(["True", "1", "true"])].copy()
    if sz_clz.n_obs == 0 or control_clz.n_obs == 0:
        raise ValueError(
            "Both Sz and healthy-control cohorts must contain at least one predicted clozapine-responsive cell. "
            f"Observed Sz={sz_clz.n_obs}, healthy_control={control_clz.n_obs}."
        )
    sz_clz = make_unique_obs_names(sz_clz, "Sz")
    control_clz = make_unique_obs_names(control_clz, "HC")
    combined = anndata.concat(
        {"Sz": sz_clz, "healthy_control": control_clz},
        label="cohort_key",
        join="inner",
        merge="same",
        index_unique=None,
    )
    combined.obs["diagnosis_group"] = pd.Categorical(
        combined.obs["diagnosis_group"], categories=["healthy_control", "Sz"]
    )
    return combined


def run_deg(clz_adata: anndata.AnnData, output_dir: Path, method: str) -> pd.DataFrame:
    """Run Scanpy DEG for Sz vs healthy control among clozapine-responsive cells."""
    deg_adata = clz_adata.copy()
    if sparse.issparse(deg_adata.X):
        deg_adata.X = deg_adata.X.tocsr()
    sc.tl.rank_genes_groups(
        deg_adata,
        groupby="diagnosis_group",
        groups=["Sz"],
        reference="healthy_control",
        method=method,
        corr_method="benjamini-hochberg",
    )
    deg_df = sc.get.rank_genes_groups_df(deg_adata, group="Sz")
    deg_df = deg_df.rename(
        columns={
            "names": "gene",
            "scores": "score",
            "logfoldchanges": "log2FC_Sz_vs_healthy_control",
            "pvals": "p_value",
            "pvals_adj": "p_adj_bh",
        }
    )
    deg_df.to_csv(output_dir / "clz_responsive_cells_DEG_Sz_vs_healthy_control.csv", index=False)
    return deg_df


def save_prediction_summaries(sz_adata: anndata.AnnData, control_adata: anndata.AnnData, output_dir: Path) -> None:
    """Save cell-level annotations and clozapine-responsive counts."""
    summary = pd.DataFrame(
        [
            {
                "diagnosis_group": "Sz",
                "n_cells_total": int(sz_adata.n_obs),
                "n_clz_responsive": int(sz_adata.obs["is_clz_selective"].astype(str).isin(["True", "1", "true"]).sum()),
            },
            {
                "diagnosis_group": "healthy_control",
                "n_cells_total": int(control_adata.n_obs),
                "n_clz_responsive": int(control_adata.obs["is_clz_selective"].astype(str).isin(["True", "1", "true"]).sum()),
            },
        ]
    )
    summary["pct_clz_responsive"] = summary["n_clz_responsive"] / summary["n_cells_total"] * 100
    summary.to_csv(output_dir / "clz_response_prediction_summary.csv", index=False)

    response_cols = [
        col for col in sz_adata.obs.columns.union(control_adata.obs.columns) if col.startswith("cAMP_") or col.startswith("Ca_")
    ]
    obs_cols = ["diagnosis_group", "is_clz_selective", "is_clz_activated", "is_clz_inhibited", *response_cols]
    pd.concat(
        [sz_adata.obs.reindex(columns=obs_cols), control_adata.obs.reindex(columns=obs_cols)],
        axis=0,
    ).to_csv(output_dir / "clz_response_cell_annotations.csv")


def plot_response_counts(output_dir: Path) -> None:
    """Plot number and percentage of predicted clozapine-responsive cells."""
    summary = pd.read_csv(output_dir / "clz_response_prediction_summary.csv")
    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.6))
    sns.barplot(data=summary, x="diagnosis_group", y="n_clz_responsive", ax=axes[0], palette=["#4c78a8", "#e45756"])
    sns.barplot(data=summary, x="diagnosis_group", y="pct_clz_responsive", ax=axes[1], palette=["#4c78a8", "#e45756"])
    axes[0].set_xlabel("")
    axes[0].set_ylabel("Clozapine-responsive cells")
    axes[1].set_xlabel("")
    axes[1].set_ylabel("Clozapine-responsive cells (%)")
    for ax in axes:
        ax.tick_params(axis="x", rotation=20)
    fig.tight_layout()
    fig.savefig(output_dir / "clz_response_prediction_summary.pdf")
    plt.close(fig)


def plot_deg_volcano(deg_df: pd.DataFrame, output_dir: Path, fdr_cutoff: float, logfc_cutoff: float) -> None:
    """Save an Illustrator-editable DEG volcano plot."""
    plot_df = deg_df.copy()
    plot_df["neg_log10_fdr"] = -np.log10(plot_df["p_adj_bh"].fillna(1.0).clip(lower=1e-300))
    plot_df["significant"] = (plot_df["p_adj_bh"] < fdr_cutoff) & (plot_df["log2FC_Sz_vs_healthy_control"].abs() >= logfc_cutoff)
    fig, ax = plt.subplots(figsize=(6.2, 5.0))
    sns.scatterplot(
        data=plot_df,
        x="log2FC_Sz_vs_healthy_control",
        y="neg_log10_fdr",
        hue="significant",
        palette={True: "#d62728", False: "#7f7f7f"},
        s=18,
        linewidth=0,
        ax=ax,
        legend=False,
    )
    for _, row in plot_df.sort_values("p_adj_bh", na_position="last").head(15).iterrows():
        ax.text(row["log2FC_Sz_vs_healthy_control"], row["neg_log10_fdr"], row["gene"], fontsize=7)
    ax.axhline(-np.log10(fdr_cutoff), color="black", linewidth=0.8, linestyle="--")
    ax.axvline(logfc_cutoff, color="black", linewidth=0.8, linestyle="--")
    ax.axvline(-logfc_cutoff, color="black", linewidth=0.8, linestyle="--")
    ax.set_xlabel("log2FC (Sz / healthy control), clozapine-responsive cells")
    ax.set_ylabel("-log10(FDR)")
    ax.set_title("DEG in predicted clozapine-responsive cells")
    fig.tight_layout()
    fig.savefig(output_dir / "clz_responsive_cells_DEG_volcano.pdf")
    plt.close(fig)


def plot_top_gene_heatmap(clz_adata: anndata.AnnData, deg_df: pd.DataFrame, output_dir: Path, top_n: int) -> None:
    """Plot mean expression z-scores for top DE genes by diagnosis group."""
    top_genes = deg_df.dropna(subset=["p_adj_bh"]).head(top_n)["gene"].astype(str).tolist()
    top_genes = [gene for gene in top_genes if gene in clz_adata.var_names]
    if not top_genes:
        return
    expr = sc.get.obs_df(clz_adata, keys=[*top_genes, "diagnosis_group"])
    mean_expr = expr.groupby("diagnosis_group", observed=True)[top_genes].mean().T
    z = mean_expr.sub(mean_expr.mean(axis=1), axis=0).div(mean_expr.std(axis=1).replace(0, np.nan), axis=0).fillna(0)
    fig_height = max(4.0, 0.22 * len(z) + 1.4)
    fig, ax = plt.subplots(figsize=(4.2, fig_height))
    sns.heatmap(z, cmap="vlag", center=0, ax=ax, cbar_kws={"label": "row z-score"})
    ax.set_xlabel("")
    ax.set_ylabel("Top DEG")
    fig.tight_layout()
    fig.savefig(output_dir / "clz_responsive_cells_top_DEG_heatmap.pdf")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sz-h5ad", default=DEFAULT_SZ_H5AD, help="Schizophrenia AnnData h5ad path.")
    parser.add_argument("--control-h5ad", default=DEFAULT_CONTROL_H5AD, help="Healthy-control AnnData h5ad path.")
    parser.add_argument("--out-dir", default=DEFAULT_OUTPUT_DIR, help="Directory for CSV/PDF outputs.")
    parser.add_argument("--drug-conc", type=float, default=1.0, help="Drug concentration passed to ct.calc_drug_response.")
    parser.add_argument("--selectivity-threshold", type=float, default=1.5, help="Threshold for ct.calc_clz_selective_cell.")
    parser.add_argument("--deg-method", default="wilcoxon", choices=["wilcoxon", "t-test", "t-test_overestim_var", "logreg"])
    parser.add_argument("--is-gpu", action="store_true", help="Use GPU preprocessing path in calculation_tool.")
    parser.add_argument("--label-filter", default=None, help="Optional adata.obs['label'] value passed to preprocessing.")
    parser.add_argument("--add-marker", action="append", default=[], help="Additional marker/gene to retain during preprocessing; can be repeated.")
    parser.add_argument("--top-n-heatmap", type=int, default=40, help="Number of top DE genes for heatmap.")
    parser.add_argument("--fdr-cutoff", type=float, default=0.05, help="FDR cutoff drawn on volcano plot.")
    parser.add_argument("--logfc-cutoff", type=float, default=0.25, help="Absolute log2FC cutoff drawn on volcano plot.")
    parser.add_argument("--save-clz-h5ad", action="store_true", help="Save concatenated clozapine-responsive AnnData.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    configure_plotting()
    output_dir = Path(args.out_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    sz_adata, sz_gpcr_df = preprocess_and_predict(
        args.sz_h5ad, "Sz", args.drug_conc, args.selectivity_threshold, args.is_gpu, args.label_filter, args.add_marker
    )
    control_adata, control_gpcr_df = preprocess_and_predict(
        args.control_h5ad,
        "healthy_control",
        args.drug_conc,
        args.selectivity_threshold,
        args.is_gpu,
        args.label_filter,
        args.add_marker,
    )
    sz_gpcr_df.to_csv(output_dir / "Sz_GPCR_raw_expression.csv")
    control_gpcr_df.to_csv(output_dir / "healthy_control_GPCR_raw_expression.csv")

    save_prediction_summaries(sz_adata, control_adata, output_dir)
    plot_response_counts(output_dir)

    clz_adata = concatenate_clz_responsive(sz_adata, control_adata)
    deg_df = run_deg(clz_adata, output_dir, args.deg_method)
    plot_deg_volcano(deg_df, output_dir, args.fdr_cutoff, args.logfc_cutoff)
    plot_top_gene_heatmap(clz_adata, deg_df, output_dir, args.top_n_heatmap)

    if args.save_clz_h5ad:
        clz_adata.write(output_dir / "clz_responsive_cells_Sz_and_healthy_control.h5ad")

    print(f"Saved clozapine-response DEG outputs to: {output_dir}")


if __name__ == "__main__":
    main()
