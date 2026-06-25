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
import shutil
import subprocess
from pathlib import Path

import anndata
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import sparse, stats
from statsmodels.stats.multitest import multipletests

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


def add_clz_selective_labels(
    adata: anndata.AnnData,
    drug_list: list[str],
    selectivity_threshold: float,
) -> tuple[anndata.AnnData, int]:
    """Add clozapine-selective labels without opening Scanpy plot windows.

    ``calculation_tool.calc_clz_selective_cell`` also calls ``sc.pl.umap`` at the
    end, which can open an interactive window and block this batch script. This
    function keeps the same label definitions but intentionally does not plot.
    """
    adata = adata.copy()
    clz_col = "cAMP_CLOZAPINE"
    if clz_col not in adata.obs.columns:
        raise KeyError(f"{clz_col} is required before clozapine-selective cell labeling.")

    adata.obs["is_clz_activated"] = (adata.obs[clz_col] > 10).astype("category")
    adata.obs["is_clz_inhibited"] = (adata.obs[clz_col] < -10).astype("category")

    drug_cols = [f"cAMP_{drug}" for drug in drug_list if drug != "CLOZAPINE"]
    missing_drug_cols = [col for col in drug_cols if col not in adata.obs.columns]
    if missing_drug_cols:
        raise KeyError(f"Missing drug-response columns: {missing_drug_cols}")

    for col in [*drug_cols, clz_col]:
        adata.obs[col] = adata.obs[col].astype("float32")

    epsilon = 1e-9
    adata.obs["cAMP_mean_other_than_czp"] = adata.obs[drug_cols].mean(axis=1) + epsilon
    adata.obs["cAMP_clz_selectivity"] = (
        adata.obs[clz_col] ** 2 / adata.obs["cAMP_mean_other_than_czp"] ** 2
    )
    adata.obs["is_clz_selective"] = (
        (adata.obs["cAMP_clz_selectivity"] > selectivity_threshold) & (adata.obs[clz_col] > 0)
    ).astype("category")

    num_clz_selective = int(adata.obs["is_clz_selective"].astype(bool).sum())
    print("clz selective cells")
    print("# of clz selective cells:", adata.obs["is_clz_selective"].value_counts())
    return adata, num_clz_selective


def check_hvg_did_not_subset_genes(
    adata: anndata.AnnData,
    gpcr_df: pd.DataFrame,
    expected_gpcrs: list[str],
    group_label: str,
    output_dir: Path,
) -> None:
    """Validate that preprocessing kept non-HVG genes and retained GPCR inputs."""
    hvg_count = int(adata.var["highly_variable"].sum()) if "highly_variable" in adata.var.columns else np.nan
    has_non_hvg_genes = bool("highly_variable" not in adata.var.columns or (~adata.var["highly_variable"]).any())
    gpcr_raw_genes = {str(col).removesuffix("_raw").upper() for col in gpcr_df.columns}
    var_genes = {str(gene).upper() for gene in adata.var_names}
    rows = []
    for gpcr in expected_gpcrs:
        gpcr_upper = str(gpcr).upper()
        rows.append(
            {
                "diagnosis_group": group_label,
                "gpcr": gpcr_upper,
                "in_processed_var_names": gpcr_upper in var_genes,
                "in_gpcr_raw_expression": gpcr_upper in gpcr_raw_genes,
            }
        )
    report_df = pd.DataFrame(rows)
    report_df["retained_for_drug_response"] = report_df["in_gpcr_raw_expression"]
    report_df.to_csv(output_dir / f"{group_label}_gpcr_retention_check.csv", index=False)

    summary_df = pd.DataFrame(
        [
            {
                "diagnosis_group": group_label,
                "n_vars_after_preprocess": int(adata.n_vars),
                "n_highly_variable_genes": hvg_count,
                "has_non_hvg_genes_after_preprocess": has_non_hvg_genes,
                "n_expected_gpcrs": len(expected_gpcrs),
                "n_gpcrs_in_processed_var_names": int(report_df["in_processed_var_names"].sum()),
                "n_gpcrs_in_gpcr_raw_expression": int(report_df["in_gpcr_raw_expression"].sum()),
            }
        ]
    )
    summary_path = output_dir / "hvg_and_gpcr_retention_summary.csv"
    summary_df.to_csv(summary_path, mode="a", header=not summary_path.exists(), index=False)

    if not has_non_hvg_genes:
        raise RuntimeError(
            f"{group_label}: preprocessing appears to contain only highly variable genes. "
            "GPCRs may have been dropped by HVG subsetting."
        )
    missing_for_response = report_df.loc[~report_df["retained_for_drug_response"], "gpcr"].tolist()
    if len(missing_for_response) == len(expected_gpcrs):
        raise RuntimeError(
            f"{group_label}: no GPCRs were retained in GPCR_df; drug-response calculation cannot proceed."
        )
    if missing_for_response:
        print(
            f"[WARN] {group_label}: these expected GPCRs were not retained in GPCR_df "
            f"(absent before/after filtering, not HVG-subsetted): {missing_for_response}"
        )


def preprocess_and_predict(
    adata_path: str | Path,
    group_label: str,
    drug_conc: float,
    selectivity_threshold: float,
    is_gpu: bool,
    output_dir: Path,
    label_filter: str | None = None,
    add_markers: list[str] | None = None,
) -> tuple[anndata.AnnData, pd.DataFrame]:
    """Preprocess one cohort and add clozapine-response annotations."""
    d_r_mtx, gpcr_type_df, drug_list, gpcr_list = ct.load_parameters()
    adata, gpcr_df = ct.preprocess_adata_in_bulk(
        adata_path=str(adata_path),
        label=label_filter,
        add_markers=add_markers,
        is_gpu=is_gpu,
        run_drug_response=False,
    )
    check_hvg_did_not_subset_genes(adata, gpcr_df, gpcr_list, group_label, output_dir)
    adata = ct.calc_drug_response(adata, gpcr_df, gpcr_type_df, drug_list, d_r_mtx, drug_conc)
    adata, _ = add_clz_selective_labels(adata, drug_list, selectivity_threshold)
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


def _obs_marker_score(adata: anndata.AnnData, genes: list[str]) -> pd.Series:
    """Score marker raw-count columns stored in ``adata.obs``."""
    score = pd.Series(0.0, index=adata.obs_names)
    for gene in genes:
        raw_col = f"{gene}_raw"
        if raw_col in adata.obs.columns:
            score = score + pd.to_numeric(adata.obs[raw_col], errors="coerce").fillna(0.0)
    return score


def annotate_broad_cell_types(clz_adata: anndata.AnnData, output_dir: Path) -> anndata.AnnData:
    """Annotate broad clozapine-responsive cell classes from retained marker raw counts."""
    adata = clz_adata.copy()
    inhibitory_score = _obs_marker_score(adata, ["GAD1", "GAD2"])
    excitatory_score = _obs_marker_score(adata, ["SLC17A7", "SLC17A6", "SATB2"])
    non_neuronal_score = _obs_marker_score(adata, ["CX3CR1", "CLDN5", "GLUL", "NDRG2", "PCDH15", "PLP1", "MBP"])
    neuronal_score = _obs_marker_score(adata, ["SNAP25"])

    broad_type = pd.Series("unclassified", index=adata.obs_names, dtype="object")
    broad_type.loc[(non_neuronal_score > 0) & (neuronal_score <= 0)] = "non_neuronal"
    broad_type.loc[excitatory_score > 0] = "excitatory"
    broad_type.loc[inhibitory_score > 0] = "inhibitory"
    broad_type.loc[
        (broad_type == "unclassified") & (non_neuronal_score > excitatory_score) & (non_neuronal_score > inhibitory_score)
    ] = "non_neuronal"

    adata.obs["broad_cell_type"] = pd.Categorical(
        broad_type,
        categories=["excitatory", "inhibitory", "non_neuronal", "unclassified"],
    )
    marker_scores = pd.DataFrame(
        {
            "excitatory_marker_score": excitatory_score,
            "inhibitory_marker_score": inhibitory_score,
            "non_neuronal_marker_score": non_neuronal_score,
            "neuronal_marker_score": neuronal_score,
            "broad_cell_type": adata.obs["broad_cell_type"].astype(str),
            "diagnosis_group": adata.obs["diagnosis_group"].astype(str),
        },
        index=adata.obs_names,
    )
    marker_scores.to_csv(output_dir / "clz_responsive_cells_broad_cell_type_annotations.csv")
    (
        adata.obs.groupby(["broad_cell_type", "diagnosis_group"], observed=False)
        .size()
        .rename("n_cells")
        .reset_index()
        .to_csv(output_dir / "clz_responsive_cells_broad_cell_type_counts.csv", index=False)
    )
    return adata


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


def _matrix_to_dense_array(matrix) -> np.ndarray:
    """Return a dense numpy array from a dense or sparse expression matrix."""
    if sparse.issparse(matrix):
        return matrix.toarray()
    return np.asarray(matrix)


def make_pseudobulk_tables(
    clz_adata: anndata.AnnData,
    donor_col: str,
    output_dir: Path,
    layer: str = "raw_counts",
    target_sum: float = 1_000_000.0,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Pool clozapine-responsive cells by donor so each donor is one sample."""
    if donor_col not in clz_adata.obs.columns:
        raise KeyError(
            f"Donor column '{donor_col}' is not present in clz_adata.obs. "
            "Specify the correct donor/sample column with --donor-col."
        )
    if "diagnosis_group" not in clz_adata.obs.columns:
        raise KeyError("diagnosis_group is required for pseudobulk analysis.")

    source_matrix = clz_adata.layers[layer] if layer in clz_adata.layers else clz_adata.X
    rows = []
    sample_ids = []
    metadata_rows = []
    grouped = clz_adata.obs.groupby(["diagnosis_group", donor_col], observed=True)
    for (diagnosis_group, donor_id), obs_subset in grouped:
        cell_positions = clz_adata.obs_names.get_indexer(obs_subset.index)
        if (cell_positions < 0).any():
            raise RuntimeError("Could not map pseudobulk cell names back to AnnData positions.")
        pooled = source_matrix[cell_positions, :].sum(axis=0)
        pooled_array = _matrix_to_dense_array(pooled).ravel().astype(float)
        sample_id = f"{diagnosis_group}__{donor_id}"
        rows.append(pooled_array)
        sample_ids.append(sample_id)
        metadata_rows.append(
            {
                "pseudobulk_sample": sample_id,
                "diagnosis_group": diagnosis_group,
                "donor": donor_id,
                "n_clz_responsive_cells": int(len(obs_subset)),
                "library_size": float(pooled_array.sum()),
                "expression_layer": layer if layer in clz_adata.layers else "X",
            }
        )

    pseudobulk_counts = pd.DataFrame(rows, index=sample_ids, columns=clz_adata.var_names)
    pseudobulk_metadata = pd.DataFrame(metadata_rows).set_index("pseudobulk_sample")
    library_size = pseudobulk_counts.sum(axis=1).replace(0, np.nan)
    pseudobulk_cpm = pseudobulk_counts.div(library_size, axis=0) * target_sum
    pseudobulk_log_cpm = np.log1p(pseudobulk_cpm)

    pseudobulk_counts.to_csv(output_dir / "clz_responsive_cells_pseudobulk_counts_by_donor.csv")
    pseudobulk_log_cpm.to_csv(output_dir / "clz_responsive_cells_pseudobulk_logCPM_by_donor.csv")
    pseudobulk_metadata.to_csv(output_dir / "clz_responsive_cells_pseudobulk_metadata.csv")
    pseudobulk_counts.round().astype(int).T.to_csv(
        output_dir / "clz_responsive_cells_pseudobulk_counts_for_deseq2_edger.csv"
    )
    pseudobulk_metadata.to_csv(output_dir / "clz_responsive_cells_pseudobulk_coldata_for_deseq2_edger.csv")
    (output_dir / "clz_responsive_cells_pseudobulk_deseq2_edger_notes.txt").write_text(
        "Counts file: clz_responsive_cells_pseudobulk_counts_for_deseq2_edger.csv\n"
        "  - rows are genes, columns are donor-level pseudobulk samples.\n"
        "Metadata file: clz_responsive_cells_pseudobulk_coldata_for_deseq2_edger.csv\n"
        "  - rows are donor-level pseudobulk samples.\n"
        "Suggested design for edgeR/DESeq2: ~ diagnosis_group\n",
        encoding="utf-8",
    )
    return pseudobulk_counts, pseudobulk_log_cpm, pseudobulk_metadata


def run_pseudobulk_deg(
    pseudobulk_counts: pd.DataFrame,
    pseudobulk_log_cpm: pd.DataFrame,
    pseudobulk_metadata: pd.DataFrame,
    output_dir: Path,
    min_cpm: float,
    min_samples: int,
) -> pd.DataFrame:
    """Run donor-level Welch t-tests on pseudobulked clozapine-responsive cells."""
    sz_samples = pseudobulk_metadata.index[pseudobulk_metadata["diagnosis_group"].astype(str) == "Sz"]
    control_samples = pseudobulk_metadata.index[pseudobulk_metadata["diagnosis_group"].astype(str) == "healthy_control"]
    if len(sz_samples) == 0 or len(control_samples) == 0:
        raise ValueError(
            "Pseudobulk DEG requires at least one Sz donor and one healthy-control donor "
            f"(observed Sz={len(sz_samples)}, healthy_control={len(control_samples)})."
        )

    pseudobulk_cpm = np.expm1(pseudobulk_log_cpm)
    if min_cpm > 0 and min_samples > 0:
        keep_genes = (pseudobulk_cpm >= min_cpm).sum(axis=0) >= min_samples
    else:
        keep_genes = pd.Series(True, index=pseudobulk_log_cpm.columns)
    gene_filter_df = pd.DataFrame(
        {
            "gene": pseudobulk_log_cpm.columns,
            "n_samples_ge_min_cpm": (pseudobulk_cpm >= min_cpm).sum(axis=0).to_numpy(),
            "min_cpm": min_cpm,
            "min_samples": min_samples,
            "kept_for_pseudobulk_deg": keep_genes.to_numpy(dtype=bool),
        }
    )
    gene_filter_df.to_csv(output_dir / "clz_responsive_cells_pseudobulk_gene_filter.csv", index=False)

    rows = []
    for gene in pseudobulk_log_cpm.columns[keep_genes]:
        sz_log = pseudobulk_log_cpm.loc[sz_samples, gene].astype(float).to_numpy()
        control_log = pseudobulk_log_cpm.loc[control_samples, gene].astype(float).to_numpy()
        sz_cpm = np.expm1(sz_log)
        control_cpm = np.expm1(control_log)

        if len(sz_log) >= 2 and len(control_log) >= 2:
            statistic, p_value = stats.ttest_ind(sz_log, control_log, equal_var=False, nan_policy="omit")
        else:
            statistic, p_value = np.nan, np.nan

        mean_sz_cpm = float(np.nanmean(sz_cpm)) if len(sz_cpm) else np.nan
        mean_control_cpm = float(np.nanmean(control_cpm)) if len(control_cpm) else np.nan
        rows.append(
            {
                "gene": gene,
                "n_donors_Sz": int(len(sz_samples)),
                "n_donors_healthy_control": int(len(control_samples)),
                "mean_logCPM_Sz": float(np.nanmean(sz_log)) if len(sz_log) else np.nan,
                "mean_logCPM_healthy_control": float(np.nanmean(control_log)) if len(control_log) else np.nan,
                "mean_CPM_Sz": mean_sz_cpm,
                "mean_CPM_healthy_control": mean_control_cpm,
                "log2FC_Sz_vs_healthy_control": float(
                    np.log2((mean_sz_cpm + 0.5) / (mean_control_cpm + 0.5))
                ),
                "welch_t_statistic": float(statistic) if not np.isnan(statistic) else np.nan,
                "p_value": float(p_value) if not np.isnan(p_value) else np.nan,
            }
        )

    result_columns = [
        "gene",
        "n_donors_Sz",
        "n_donors_healthy_control",
        "mean_logCPM_Sz",
        "mean_logCPM_healthy_control",
        "mean_CPM_Sz",
        "mean_CPM_healthy_control",
        "log2FC_Sz_vs_healthy_control",
        "welch_t_statistic",
        "p_value",
    ]
    deg_df = pd.DataFrame(rows, columns=result_columns)
    deg_df["p_adj_bh"] = np.nan
    valid = deg_df["p_value"].notna()
    if valid.any():
        deg_df.loc[valid, "p_adj_bh"] = multipletests(deg_df.loc[valid, "p_value"], method="fdr_bh")[1]
    deg_df["abs_log2FC"] = deg_df["log2FC_Sz_vs_healthy_control"].abs()
    deg_df["neg_log10_p_value"] = -np.log10(deg_df["p_value"].fillna(1.0).clip(lower=1e-300))
    deg_df["neg_log10_p_adj"] = -np.log10(deg_df["p_adj_bh"].fillna(1.0).clip(lower=1e-300))
    deg_df["signed_neg_log10_p_value"] = np.sign(deg_df["log2FC_Sz_vs_healthy_control"]) * deg_df["neg_log10_p_value"]
    deg_df["signed_neg_log10_p_adj"] = np.sign(deg_df["log2FC_Sz_vs_healthy_control"]) * deg_df["neg_log10_p_adj"]
    deg_df = deg_df.sort_values(["p_adj_bh", "p_value", "gene"], na_position="last").reset_index(drop=True)
    deg_df.to_csv(output_dir / "clz_responsive_cells_pseudobulk_DEG_Sz_vs_healthy_control.csv", index=False)
    return deg_df


def write_edger_deseq2_scripts(output_dir: Path) -> tuple[Path, Path]:
    """Write runnable R scripts for edgeR and DESeq2 pseudobulk DEG."""
    edge_r_script = output_dir / "run_pseudobulk_edgeR.R"
    deseq2_script = output_dir / "run_pseudobulk_DESeq2.R"
    edge_r_script.write_text(
        r'''
suppressPackageStartupMessages({
  library(edgeR)
})

out_dir <- normalizePath(".", mustWork = TRUE)
counts_path <- file.path(out_dir, "clz_responsive_cells_pseudobulk_counts_for_deseq2_edger.csv")
coldata_path <- file.path(out_dir, "clz_responsive_cells_pseudobulk_coldata_for_deseq2_edger.csv")

counts <- read.csv(counts_path, row.names = 1, check.names = FALSE)
coldata <- read.csv(coldata_path, row.names = 1, check.names = FALSE)
counts <- round(as.matrix(counts))
counts <- counts[, rownames(coldata), drop = FALSE]
coldata$diagnosis_group <- relevel(factor(coldata$diagnosis_group), ref = "healthy_control")

y <- DGEList(counts = counts, samples = coldata, group = coldata$diagnosis_group)
design <- model.matrix(~ diagnosis_group, data = coldata)
keep <- filterByExpr(y, design = design)
write.csv(
  data.frame(gene = rownames(y), kept_by_filterByExpr = keep),
  file.path(out_dir, "clz_responsive_cells_pseudobulk_edgeR_gene_filter.csv"),
  row.names = FALSE
)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = "diagnosis_groupSz")
res <- topTags(qlf, n = Inf, sort.by = "PValue")$table
res$gene <- rownames(res)
res <- res[, c("gene", setdiff(colnames(res), "gene"))]
write.csv(res, file.path(out_dir, "clz_responsive_cells_pseudobulk_edgeR_DEG_Sz_vs_healthy_control.csv"), row.names = FALSE)
'''.lstrip(),
        encoding="utf-8",
    )
    deseq2_script.write_text(
        r'''
suppressPackageStartupMessages({
  library(DESeq2)
})

out_dir <- normalizePath(".", mustWork = TRUE)
counts_path <- file.path(out_dir, "clz_responsive_cells_pseudobulk_counts_for_deseq2_edger.csv")
coldata_path <- file.path(out_dir, "clz_responsive_cells_pseudobulk_coldata_for_deseq2_edger.csv")

counts <- read.csv(counts_path, row.names = 1, check.names = FALSE)
coldata <- read.csv(coldata_path, row.names = 1, check.names = FALSE)
counts <- round(as.matrix(counts))
counts <- counts[, rownames(coldata), drop = FALSE]
coldata$diagnosis_group <- relevel(factor(coldata$diagnosis_group), ref = "healthy_control")

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ diagnosis_group
)
keep <- rowSums(counts(dds)) >= 10
write.csv(
  data.frame(gene = rownames(dds), kept_by_rowSum10 = keep),
  file.path(out_dir, "clz_responsive_cells_pseudobulk_DESeq2_gene_filter.csv"),
  row.names = FALSE
)
dds <- dds[keep, ]
dds <- DESeq(dds)
res <- results(dds, contrast = c("diagnosis_group", "Sz", "healthy_control"))
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df <- res_df[, c("gene", setdiff(colnames(res_df), "gene"))]
res_df <- res_df[order(res_df$padj, res_df$pvalue, na.last = TRUE), ]
write.csv(res_df, file.path(out_dir, "clz_responsive_cells_pseudobulk_DESeq2_DEG_Sz_vs_healthy_control.csv"), row.names = FALSE)
'''.lstrip(),
        encoding="utf-8",
    )
    return edge_r_script, deseq2_script


def run_r_script(script_path: Path) -> None:
    """Run an R script with Rscript in the script directory."""
    rscript = shutil.which("Rscript")
    if rscript is None:
        raise RuntimeError("Rscript was not found on PATH; install R or run the generated R script manually.")
    subprocess.run([rscript, str(script_path.name)], cwd=script_path.parent, check=True)


def parse_gmt(gmt_path: str | Path) -> dict[str, set[str]]:
    """Read a GMT gene-set file."""
    gene_sets: dict[str, set[str]] = {}
    with open(gmt_path, encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3:
                gene_sets[fields[0]] = {gene.upper() for gene in fields[2:] if gene}
    return gene_sets


def _weighted_enrichment_score(ranked_genes: list[str], rank_values: np.ndarray, gene_set: set[str]) -> float:
    """Compute a weighted preranked GSEA enrichment score."""
    hits = np.array([gene in gene_set for gene in ranked_genes], dtype=bool)
    n_hits = int(hits.sum())
    n_genes = len(ranked_genes)
    if n_hits == 0 or n_hits == n_genes:
        return np.nan
    hit_weights = np.abs(rank_values[hits])
    hit_weights_sum = hit_weights.sum()
    if hit_weights_sum == 0:
        hit_weights = np.ones(n_hits) / n_hits
    else:
        hit_weights = hit_weights / hit_weights_sum
    miss_weight = 1.0 / (n_genes - n_hits)
    running = np.empty(n_genes, dtype=float)
    running[hits] = hit_weights
    running[~hits] = -miss_weight
    cumulative = np.cumsum(running)
    max_es = cumulative.max()
    min_es = cumulative.min()
    return float(max_es if abs(max_es) >= abs(min_es) else min_es)


def run_preranked_gsea(
    deg_df: pd.DataFrame,
    gmt_path: str | Path,
    output_path: Path,
    min_size: int,
    max_size: int,
    permutations: int,
    seed: int,
) -> pd.DataFrame:
    """Run a lightweight preranked GSEA using DEG direction and p-value rank."""
    if deg_df.empty:
        result = pd.DataFrame()
        result.to_csv(output_path, index=False)
        return result

    rank_df = deg_df.dropna(subset=["gene", "log2FC_Sz_vs_healthy_control"]).copy()
    p_col = "p_value" if "p_value" in rank_df.columns else "pvals"
    rank_df[p_col] = pd.to_numeric(rank_df[p_col], errors="coerce").fillna(1.0)
    rank_df["rank_metric"] = (
        np.sign(pd.to_numeric(rank_df["log2FC_Sz_vs_healthy_control"], errors="coerce"))
        * -np.log10(rank_df[p_col].clip(lower=1e-300))
    )
    rank_df = rank_df.replace([np.inf, -np.inf], np.nan).dropna(subset=["rank_metric"])
    rank_df["gene_upper"] = rank_df["gene"].astype(str).str.upper()
    rank_df = rank_df.sort_values("rank_metric", ascending=False).drop_duplicates("gene_upper")
    ranked_genes = rank_df["gene_upper"].tolist()
    rank_values = rank_df["rank_metric"].to_numpy(dtype=float)
    gene_sets = parse_gmt(gmt_path)
    rng = np.random.default_rng(seed)

    rows = []
    for pathway, genes in gene_sets.items():
        overlap = set(ranked_genes).intersection(genes)
        if len(overlap) < min_size or len(overlap) > max_size:
            continue
        es = _weighted_enrichment_score(ranked_genes, rank_values, overlap)
        if np.isnan(es):
            continue
        null_es = np.array(
            [
                _weighted_enrichment_score(ranked_genes, rng.permutation(rank_values), overlap)
                for _ in range(permutations)
            ],
            dtype=float,
        )
        if es >= 0:
            same_tail = null_es[null_es >= 0]
            p_value = (np.sum(same_tail >= es) + 1) / (len(same_tail) + 1) if len(same_tail) else 1.0
            nes = es / np.mean(same_tail) if len(same_tail) and np.mean(same_tail) != 0 else np.nan
        else:
            same_tail = null_es[null_es < 0]
            p_value = (np.sum(same_tail <= es) + 1) / (len(same_tail) + 1) if len(same_tail) else 1.0
            nes = es / abs(np.mean(same_tail)) if len(same_tail) and np.mean(same_tail) != 0 else np.nan
        rows.append(
            {
                "pathway": pathway,
                "size": len(overlap),
                "ES": es,
                "NES": nes,
                "p_value": p_value,
                "leading_edge_genes": ";".join([gene for gene in ranked_genes if gene in overlap][:50]),
            }
        )

    result = pd.DataFrame(rows)
    if not result.empty:
        result["p_adj_bh"] = multipletests(result["p_value"], method="fdr_bh")[1]
        result = result.sort_values(["p_adj_bh", "p_value", "pathway"]).reset_index(drop=True)
    result.to_csv(output_path, index=False)
    return result


def plot_gsea_dotplot(gsea_df: pd.DataFrame, output_path: Path, title: str, top_n: int = 20) -> None:
    """Plot top GSEA pathways as an Illustrator-editable PDF."""
    if gsea_df.empty:
        return
    plot_df = gsea_df.sort_values(["p_adj_bh", "p_value"]).head(top_n).copy()
    plot_df["neg_log10_fdr"] = -np.log10(plot_df["p_adj_bh"].fillna(1.0).clip(lower=1e-300))
    fig_height = max(4.0, 0.28 * len(plot_df) + 1.5)
    fig, ax = plt.subplots(figsize=(7.0, fig_height))
    sns.scatterplot(
        data=plot_df,
        x="NES",
        y="pathway",
        size="size",
        hue="neg_log10_fdr",
        palette="viridis",
        ax=ax,
    )
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_title(title)
    ax.set_xlabel("Normalized enrichment score (NES)")
    ax.set_ylabel("")
    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)


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


def plot_pseudobulk_deg_volcano(
    deg_df: pd.DataFrame,
    output_dir: Path,
    fdr_cutoff: float,
    logfc_cutoff: float,
) -> None:
    """Save an Illustrator-editable pseudobulk DEG volcano plot."""
    plot_df = deg_df.copy()
    plot_df["neg_log10_fdr"] = -np.log10(plot_df["p_adj_bh"].fillna(1.0).clip(lower=1e-300))
    plot_df["significant"] = (plot_df["p_adj_bh"] < fdr_cutoff) & (
        plot_df["log2FC_Sz_vs_healthy_control"].abs() >= logfc_cutoff
    )
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
        if np.isfinite(row["neg_log10_fdr"]):
            ax.text(row["log2FC_Sz_vs_healthy_control"], row["neg_log10_fdr"], row["gene"], fontsize=7)
    ax.axhline(-np.log10(fdr_cutoff), color="black", linewidth=0.8, linestyle="--")
    ax.axvline(logfc_cutoff, color="black", linewidth=0.8, linestyle="--")
    ax.axvline(-logfc_cutoff, color="black", linewidth=0.8, linestyle="--")
    ax.set_xlabel("log2FC (Sz / healthy control), donor pseudobulk")
    ax.set_ylabel("-log10(FDR)")
    ax.set_title("Pseudobulk DEG in clozapine-responsive cells")
    fig.tight_layout()
    fig.savefig(output_dir / "clz_responsive_cells_pseudobulk_DEG_volcano.pdf")
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


def run_stratified_cell_type_analyses(
    clz_adata: anndata.AnnData,
    output_dir: Path,
    args: argparse.Namespace,
) -> None:
    """Run DEG, pseudobulk DEG, optional R runners, and optional GSEA per broad cell type."""
    celltype_root = output_dir / "broad_cell_type_stratified"
    celltype_root.mkdir(parents=True, exist_ok=True)
    summary_rows = []
    for cell_type in ["excitatory", "inhibitory", "non_neuronal"]:
        subset = clz_adata[clz_adata.obs["broad_cell_type"].astype(str) == cell_type].copy()
        group_counts = subset.obs["diagnosis_group"].astype(str).value_counts()
        celltype_dir = celltype_root / cell_type
        celltype_dir.mkdir(parents=True, exist_ok=True)
        summary_rows.append(
            {
                "broad_cell_type": cell_type,
                "n_cells": int(subset.n_obs),
                "n_cells_Sz": int(group_counts.get("Sz", 0)),
                "n_cells_healthy_control": int(group_counts.get("healthy_control", 0)),
                "ran_analysis": bool(group_counts.get("Sz", 0) > 0 and group_counts.get("healthy_control", 0) > 0),
            }
        )
        if group_counts.get("Sz", 0) == 0 or group_counts.get("healthy_control", 0) == 0:
            continue

        deg_df = run_deg(subset, celltype_dir, args.deg_method)
        plot_deg_volcano(deg_df, celltype_dir, args.fdr_cutoff, args.logfc_cutoff)
        plot_top_gene_heatmap(subset, deg_df, celltype_dir, args.top_n_heatmap)
        if args.gmt:
            gsea_df = run_preranked_gsea(
                deg_df=deg_df,
                gmt_path=args.gmt,
                output_path=celltype_dir / "cell_level_preranked_GSEA.csv",
                min_size=args.gsea_min_size,
                max_size=args.gsea_max_size,
                permutations=args.gsea_permutations,
                seed=args.gsea_seed,
            )
            plot_gsea_dotplot(gsea_df, celltype_dir / "cell_level_preranked_GSEA_dotplot.pdf", f"{cell_type} cell-level GSEA")

        if args.skip_pseudobulk:
            continue
        pseudobulk_counts, pseudobulk_log_cpm, pseudobulk_metadata = make_pseudobulk_tables(
            clz_adata=subset,
            donor_col=args.donor_col,
            output_dir=celltype_dir,
        )
        pseudobulk_deg_df = run_pseudobulk_deg(
            pseudobulk_counts=pseudobulk_counts,
            pseudobulk_log_cpm=pseudobulk_log_cpm,
            pseudobulk_metadata=pseudobulk_metadata,
            output_dir=celltype_dir,
            min_cpm=args.pseudobulk_min_cpm,
            min_samples=args.pseudobulk_min_samples,
        )
        plot_pseudobulk_deg_volcano(pseudobulk_deg_df, celltype_dir, args.fdr_cutoff, args.logfc_cutoff)
        if args.gmt:
            gsea_df = run_preranked_gsea(
                deg_df=pseudobulk_deg_df,
                gmt_path=args.gmt,
                output_path=celltype_dir / "pseudobulk_preranked_GSEA.csv",
                min_size=args.gsea_min_size,
                max_size=args.gsea_max_size,
                permutations=args.gsea_permutations,
                seed=args.gsea_seed,
            )
            plot_gsea_dotplot(gsea_df, celltype_dir / "pseudobulk_preranked_GSEA_dotplot.pdf", f"{cell_type} pseudobulk GSEA")
        edge_r_script, deseq2_script = write_edger_deseq2_scripts(celltype_dir)
        if args.run_edger:
            run_r_script(edge_r_script)
        if args.run_deseq2:
            run_r_script(deseq2_script)

    pd.DataFrame(summary_rows).to_csv(celltype_root / "broad_cell_type_stratified_analysis_summary.csv", index=False)


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
    parser.add_argument("--donor-col", default="sample", help="obs column used as donor ID for pseudobulk; one donor becomes n=1.")
    parser.add_argument("--pseudobulk-min-cpm", type=float, default=1.0, help="Keep genes with CPM >= this value before pseudobulk DEG.")
    parser.add_argument("--pseudobulk-min-samples", type=int, default=2, help="Keep genes detected above min CPM in at least this many donor pseudobulk samples.")
    parser.add_argument("--run-edger", action="store_true", help="Run generated edgeR pseudobulk DEG script with Rscript.")
    parser.add_argument("--run-deseq2", action="store_true", help="Run generated DESeq2 pseudobulk DEG script with Rscript.")
    parser.add_argument("--skip-cell-type-stratified", action="store_true", help="Skip excitatory/inhibitory/non-neuronal stratified DEG.")
    parser.add_argument("--gmt", default=None, help="GMT gene-set file for preranked GSEA.")
    parser.add_argument("--gsea-min-size", type=int, default=10, help="Minimum overlap size for preranked GSEA.")
    parser.add_argument("--gsea-max-size", type=int, default=500, help="Maximum overlap size for preranked GSEA.")
    parser.add_argument("--gsea-permutations", type=int, default=1000, help="Number of rank permutations for preranked GSEA p-values.")
    parser.add_argument("--gsea-seed", type=int, default=0, help="Random seed for preranked GSEA permutations.")
    parser.add_argument("--skip-pseudobulk", action="store_true", help="Skip donor-level pseudobulk DEG analysis.")
    parser.add_argument("--save-clz-h5ad", action="store_true", help="Save concatenated clozapine-responsive AnnData.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    configure_plotting()
    output_dir = Path(args.out_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    hvg_summary_path = output_dir / "hvg_and_gpcr_retention_summary.csv"
    if hvg_summary_path.exists():
        hvg_summary_path.unlink()

    sz_adata, sz_gpcr_df = preprocess_and_predict(
        args.sz_h5ad,
        "Sz",
        args.drug_conc,
        args.selectivity_threshold,
        args.is_gpu,
        output_dir,
        args.label_filter,
        args.add_marker,
    )
    control_adata, control_gpcr_df = preprocess_and_predict(
        args.control_h5ad,
        "healthy_control",
        args.drug_conc,
        args.selectivity_threshold,
        args.is_gpu,
        output_dir,
        args.label_filter,
        args.add_marker,
    )
    sz_gpcr_df.to_csv(output_dir / "Sz_GPCR_raw_expression.csv")
    control_gpcr_df.to_csv(output_dir / "healthy_control_GPCR_raw_expression.csv")

    save_prediction_summaries(sz_adata, control_adata, output_dir)
    plot_response_counts(output_dir)

    clz_adata = annotate_broad_cell_types(concatenate_clz_responsive(sz_adata, control_adata), output_dir)
    deg_df = run_deg(clz_adata, output_dir, args.deg_method)
    plot_deg_volcano(deg_df, output_dir, args.fdr_cutoff, args.logfc_cutoff)
    plot_top_gene_heatmap(clz_adata, deg_df, output_dir, args.top_n_heatmap)
    if args.gmt:
        gsea_df = run_preranked_gsea(
            deg_df=deg_df,
            gmt_path=args.gmt,
            output_path=output_dir / "cell_level_preranked_GSEA.csv",
            min_size=args.gsea_min_size,
            max_size=args.gsea_max_size,
            permutations=args.gsea_permutations,
            seed=args.gsea_seed,
        )
        plot_gsea_dotplot(gsea_df, output_dir / "cell_level_preranked_GSEA_dotplot.pdf", "All clz-responsive cell-level GSEA")

    if not args.skip_pseudobulk:
        pseudobulk_counts, pseudobulk_log_cpm, pseudobulk_metadata = make_pseudobulk_tables(
            clz_adata=clz_adata,
            donor_col=args.donor_col,
            output_dir=output_dir,
        )
        pseudobulk_deg_df = run_pseudobulk_deg(
            pseudobulk_counts=pseudobulk_counts,
            pseudobulk_log_cpm=pseudobulk_log_cpm,
            pseudobulk_metadata=pseudobulk_metadata,
            output_dir=output_dir,
            min_cpm=args.pseudobulk_min_cpm,
            min_samples=args.pseudobulk_min_samples,
        )
        plot_pseudobulk_deg_volcano(
            pseudobulk_deg_df,
            output_dir,
            args.fdr_cutoff,
            args.logfc_cutoff,
        )
        edge_r_script, deseq2_script = write_edger_deseq2_scripts(output_dir)
        if args.run_edger:
            run_r_script(edge_r_script)
        if args.run_deseq2:
            run_r_script(deseq2_script)
        if args.gmt:
            gsea_df = run_preranked_gsea(
                deg_df=pseudobulk_deg_df,
                gmt_path=args.gmt,
                output_path=output_dir / "pseudobulk_preranked_GSEA.csv",
                min_size=args.gsea_min_size,
                max_size=args.gsea_max_size,
                permutations=args.gsea_permutations,
                seed=args.gsea_seed,
            )
            plot_gsea_dotplot(gsea_df, output_dir / "pseudobulk_preranked_GSEA_dotplot.pdf", "All clz-responsive pseudobulk GSEA")

    if not args.skip_cell_type_stratified:
        run_stratified_cell_type_analyses(clz_adata, output_dir, args)

    if args.save_clz_h5ad:
        clz_adata.write(output_dir / "clz_responsive_cells_Sz_and_healthy_control.h5ad")

    print(f"Saved clozapine-response DEG outputs to: {output_dir}")


if __name__ == "__main__":
    main()
