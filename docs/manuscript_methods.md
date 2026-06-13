# Manuscript Methods Draft

## Single-cell RNA-seq preprocessing

Single-cell RNA-seq count matrices were analyzed using an AnnData/Scanpy-based workflow. For each dataset, gene symbols were converted to uppercase and duplicate gene names were made unique before downstream analysis. When multiple matrices were analyzed together, cells were concatenated across datasets using an outer join on genes, and dataset-level metadata, including dataset identifier, species, and brain region, were retained for downstream stratification. Cells with fewer than 200 detected genes or more than 6,000 detected genes were excluded, and genes detected in at least one cell were retained. Counts were normalized to 10,000 counts per cell and log-transformed. The log-normalized expression matrix was stored as the primary layer for differential-expression analyses. Where indicated, total counts were regressed out and expression values were scaled with a maximum value of 10. Highly variable genes were annotated using the top 4,000 genes but were not used to remove genes from receptor-focused or differential-expression analyses.

Dimensionality reduction was performed using principal component analysis (PCA; 50 components). A nearest-neighbor graph was constructed using 15 neighbors and 50 principal components, followed by UMAP embedding for visualization. Louvain and Leiden community-detection algorithms were applied for unsupervised clustering. GPU-accelerated routines were used when available; otherwise, the same preprocessing parameters and statistical definitions were applied using CPU-based Scanpy functions.

## Mouse scRNA-seq analysis and pharmacological response modeling

For mouse scRNA-seq data, receptor-expression features were extracted for a curated panel of antipsychotic-relevant G-protein-coupled receptors (GPCRs), including serotonergic, dopaminergic, histaminergic, muscarinic, adrenergic, and adenosine receptors. Raw receptor counts were retained in per-cell metadata fields when available. For pharmacological modeling, the receptor-expression matrix was normalized to 10,000 counts per cell. A drug-by-receptor affinity matrix was then used to calculate model-derived cAMP- and Ca-associated responses for each cell and each drug. These response values were used as computational phenotypes for identifying clozapine-responsive and clozapine-selective cells; they should be interpreted as receptor-expression- and affinity-weighted response scores rather than direct biochemical measurements.

## cAMP-response model

The receptor panel was stratified by canonical G-protein coupling using the receptor-type annotation table. Let \(E_{i,r}\) denote the normalized expression of receptor \(r\) in cell \(i\), \(C\) denote the modeled drug concentration, and \(K_{d,r}\) denote the affinity value for drug \(d\) and receptor \(r\) from the drug-by-receptor matrix. For Gs- and Gi-coupled receptors, the drug-adjusted residual receptor signal was calculated as

\[
E^{\mathrm{drug}}_{i,r,d} = \frac{E_{i,r}}{1 + C/K_{d,r}}.
\]

Basal cAMP-related signaling was defined as the difference between total Gs-coupled and total Gi-coupled receptor expression:

\[
B^{\mathrm{cAMP}}_i = \sum_{r \in Gs} E_{i,r} - \sum_{r \in Gi} E_{i,r}.
\]

The modeled cAMP response for cell \(i\) and drug \(d\), stored as `cAMP_<DRUG>`, was calculated as the drug-adjusted Gs-minus-Gi signal minus the basal Gs-minus-Gi signal:

\[
R^{\mathrm{cAMP}}_{i,d} =
\left(\sum_{r \in Gs} \frac{E_{i,r}}{1 + C/K_{d,r}} -
      \sum_{r \in Gi} \frac{E_{i,r}}{1 + C/K_{d,r}}\right)
- B^{\mathrm{cAMP}}_i.
\]

Thus, the cAMP score represents the modeled drug-induced change from each cell's basal Gs-minus-Gi receptor-expression balance. Under this formulation, reduced contribution from Gi-coupled receptors shifts the score upward, whereas reduced contribution from Gs-coupled receptors shifts the score downward.

## Ca-response model

For Gq-coupled receptors, basal Ca-related signaling was defined as the total normalized expression of Gq-coupled receptors:

\[
B^{\mathrm{Ca}}_i = \sum_{r \in Gq} E_{i,r}.
\]

The drug-weighted Gq signal was calculated directly from the affinity matrix as \(\sum_{r \in Gq} E_{i,r}/K_{d,r}\). The modeled Ca response for cell \(i\) and drug \(d\), stored as `Ca_<DRUG>`, was therefore

\[
R^{\mathrm{Ca}}_{i,d} = \sum_{r \in Gq} \frac{E_{i,r}}{K_{d,r}} - B^{\mathrm{Ca}}_i + 10^{-4}.
\]

The \(10^{-4}\) offset was added to all Ca-response values after calculation. Ca and cAMP response scores were calculated for clozapine and for each comparator antipsychotic in the drug-by-receptor matrix, enabling cell-level comparison of clozapine-associated response patterns with those predicted for other antipsychotic drugs.

## Definition of clozapine-responsive and clozapine-selective cells

Clozapine-responsive cells were defined from the modeled cAMP response to clozapine. Cells with `cAMP_CLOZAPINE > 10` were annotated as clozapine-activated, and cells with `cAMP_CLOZAPINE < -10` were annotated as clozapine-inhibited. A separate clozapine-selective label was used to identify cells with a cAMP response preferentially associated with clozapine compared with other antipsychotic drugs.

For the selectivity calculation, all cAMP-response columns corresponding to non-clozapine drugs were first collected as comparator responses. For each cell, the mean comparator response was calculated as

\[
\overline{R}_{\mathrm{other},i} = \frac{1}{N_{\mathrm{other}}} \sum_{d \ne \mathrm{CLOZAPINE}} R^{\mathrm{cAMP}}_{i,d} + \epsilon,
\]

where \(N_{\mathrm{other}}\) is the number of non-clozapine drugs and \(\epsilon = 10^{-9}\) was added to avoid division by zero. Clozapine selectivity was then calculated as the squared ratio

\[
S_i = \frac{\left(R^{\mathrm{cAMP}}_{i,\mathrm{CLOZAPINE}}\right)^2}{\left(\overline{R}_{\mathrm{other},i}\right)^2}.
\]

Cells were classified as clozapine-selective when \(S_i\) exceeded the prespecified selectivity threshold and the clozapine cAMP response was positive:

\[
S_i > \theta \quad \mathrm{and} \quad R^{\mathrm{cAMP}}_{i,\mathrm{CLOZAPINE}} > 0.
\]

In the standardized pipeline configuration, the modeled drug concentration was set to 1,000 and the selectivity threshold \(\theta\) was set to 1.5. This definition excludes cells with large-magnitude but negative clozapine cAMP responses from the clozapine-selective class. The resulting binary clozapine-selective annotation was used for downstream differential-expression analysis, enrichment analysis, and comparison with other antipsychotic-response profiles.

## Inhibitor-cocktail cAMP-pattern analysis

Candidate inhibitor-cocktail patterns were evaluated by specifying a binary inhibition pattern across receptor columns. Receptors included in the cocktail were assigned a low effective Ki value (`Ki_inhibited = 0.01` by default), whereas receptors not included in the cocktail were assigned a high effective Ki value (`Ki_not_inhibited = 10000` by default). Pattern keys corresponded to GPCR expression columns, such as `HTR1A_raw`, and receptors absent from the pattern were treated as not inhibited. For each receptor column \(r\), an effective \(K_r\) was defined from this binary pattern, and the cAMP response was recalculated using the same Gs/Gi framework:

\[
R^{\mathrm{cocktail}}_{i} =
\left(\sum_{r \in Gs} \frac{E_{i,r}}{1 + C/K_r} -
      \sum_{r \in Gi} \frac{E_{i,r}}{1 + C/K_r}\right)
- \left(\sum_{r \in Gs} E_{i,r} - \sum_{r \in Gi} E_{i,r}\right).
\]

Only Gs- and Gi-coupled receptors present in the normalized GPCR expression matrix were included in this cocktail cAMP-pattern calculation. Cells were grouped by clozapine-selective status, and cAMP-response distributions were summarized by cell count, mean, median, and standard deviation. Differences between clozapine-selective and nonselective cells were evaluated using both Welch's two-sample t-test and a two-sided Mann-Whitney U test.

## Human scRNA-seq analysis

Human scRNA-seq datasets were processed using the same AnnData/Scanpy framework and pharmacological modeling definitions used for mouse data. Gene symbols were harmonized to uppercase, quality-control filters were applied, and expression matrices were normalized, log-transformed, embedded by PCA/UMAP, and clustered. Clozapine-selective response scores were inferred with the same receptor-affinity and receptor-coupling models after restricting calculations to receptors detected in the human expression matrix. When curated cell-type annotations were unavailable, broad cell classes were inferred from marker-gene module scores, including excitatory neurons, inhibitory neurons, astrocytes, oligodendrocyte-lineage cells, endothelial cells, and microglia. Annotated or inferred cell-type labels were used for enrichment analyses and for stratified interpretation of clozapine-selective populations.

## Differential-expression analysis

Differentially expressed genes (DEGs) between clozapine-selective and nonselective cells were identified using Scanpy's `rank_genes_groups` function with the Wilcoxon rank-sum test applied to the log-normalized expression layer. The clozapine-selective annotation was encoded as a categorical grouping variable, and gene rankings, test statistics, log fold changes, adjusted p values, and detection fractions were extracted with `scanpy.get.rank_genes_groups_df`. Unless otherwise specified, genes were considered significant at a Benjamini-Hochberg false-discovery-rate (FDR) threshold of 0.05. An additional absolute log2 fold-change threshold of 0.25 was used for volcano-plot visualization and interpretation. For receptor-focused analyses, raw receptor-count columns were preferred when present; otherwise, expression values were extracted from the AnnData expression matrix. Receptor and marker expression differences between clozapine-selective and nonselective cells were tested using two-sided Mann-Whitney U tests, followed by Benjamini-Hochberg correction across tested genes.

## Excitatory-neuron-only DEG analysis

To determine whether clozapine-selective transcriptional differences were present within excitatory neurons, a subset analysis was performed after restricting the AnnData object to cells annotated or inferred as excitatory neurons. Excitatory cells were identified from available cell-type annotations or marker-based labels using glutamatergic/excitatory markers such as SLC17A7, SLC17A6, SATB2, and related cortical excitatory-neuron markers. Within this excitatory-only subset, DEGs were recalculated between clozapine-selective and nonselective cells using the same Wilcoxon rank-sum framework, log-normalized expression layer, and multiple-testing correction described above. Excitatory-only volcano plots were generated using an FDR threshold of 0.05 and an absolute log2 fold-change threshold of 0.25 unless otherwise indicated. Selected disease-risk or mechanistic genes were labeled only by exact gene-symbol matching to avoid partial-match annotation artifacts.

## Cell-type and cluster enrichment analysis

Enrichment of clozapine-selective cells within annotated cell types, inferred cell types, or Leiden clusters was tested using two-by-two contingency tables. For each category, the table compared the number of clozapine-selective and nonselective cells inside the category with the corresponding numbers outside the category. Two-sided Fisher's exact tests were used to compute enrichment P values. Odds ratios and 95% confidence intervals were estimated using a Haldane-Anscombe correction by adding 0.5 to each contingency-table cell before calculating the log odds ratio and its standard error. P values across tested categories were corrected using the Benjamini-Hochberg FDR procedure. Enrichment results were visualized as forest plots on a log odds-ratio scale.

## Gene-set enrichment analysis

Gene-set enrichment was performed from DEG results using both over-representation and ranked-list approaches. For over-representation analysis, significant upregulated or downregulated genes were defined using the FDR and log2 fold-change thresholds described above, and each gene set was tested against the expressed-gene background with Fisher's exact test followed by Benjamini-Hochberg correction. For ranked enrichment analysis, genes were ranked using the Scanpy differential-expression score or log-fold-change statistic, and preranked enrichment was performed for Gene Ontology, Reactome, and Hallmark gene-set collections when enabled. Enrichment outputs were interpreted together with cell-type enrichment results to distinguish receptor-model-associated cellular enrichment from transcriptional pathway enrichment.

## Curated pathway enrichment analysis for clozapine-selective DEGs

A targeted enrichment analysis was additionally performed for manually curated gene sets selected a priori to support mechanistic interpretation of DEGs observed in clozapine-selective cells. All gene symbols were converted to uppercase before matching, and mouse symbols were mapped to their human-style uppercase ortholog symbols when needed so that the same curated terms could be applied consistently across mouse and human analyses. The tested universe was defined as the set of genes retained in the DEG result table after expression filtering and gene-symbol harmonization, rather than the whole genome. For each curated term, genes not present in the tested universe were retained in the reported definition of the term but were excluded from the contingency-table counts.

The curated terms and member genes were defined as follows:

| Curated term | Included genes | Rationale |
| --- | --- | --- |
| Synapse and synaptic transmission | `SYN1`, `SYN2`, `SYN3`, `SYP`, `SNAP25`, `STX1A`, `STXBP1`, `VAMP2`, `SYT1`, `DLG4`, `SHANK1`, `SHANK2`, `SHANK3`, `HOMER1`, `NRXN1`, `NRXN2`, `NRXN3`, `NLGN1`, `NLGN2`, `NLGN3`, `GRIN1`, `GRIN2A`, `GRIN2B`, `GRIA1`, `GRIA2`, `CAMK2A`, `CAMK2B` | Presynaptic vesicle, release-machinery, postsynaptic-density, adhesion, and glutamatergic synapse genes used to evaluate whether clozapine-selective DEGs were enriched for neuronal synaptic biology. |
| cAMP/PKA signaling | `ADCY1`, `ADCY2`, `ADCY3`, `ADCY5`, `ADCY6`, `ADCY7`, `ADCY8`, `ADCY9`, `PRKACA`, `PRKACB`, `PRKACG`, `PRKAR1A`, `PRKAR1B`, `PRKAR2A`, `PRKAR2B`, `PDE1A`, `PDE1B`, `PDE1C`, `PDE2A`, `PDE3A`, `PDE3B`, `PDE4A`, `PDE4B`, `PDE4C`, `PDE4D`, `CREB1`, `ATF1`, `RAPGEF3`, `RAPGEF4`, `AKAP5`, `AKAP7`, `AKAP9`, `AKAP11`, `AKAP12` | Adenylyl cyclases, PKA catalytic and regulatory subunits, phosphodiesterases, cAMP-responsive transcriptional mediators, EPAC proteins, and anchoring proteins used to evaluate transcriptional changes related to the receptor-derived cAMP model. |
| GPCR and monoamine receptor signaling | `HTR1A`, `HTR2A`, `HTR2C`, `HTR6`, `HTR7`, `DRD1`, `DRD2`, `DRD3`, `DRD4`, `HRH1`, `HRH3`, `CHRM1`, `CHRM2`, `CHRM3`, `CHRM4`, `CHRM5`, `ADRA1A`, `ADRA2A`, `ADORA2A`, `GNAI1`, `GNAI2`, `GNAI3`, `GNAS`, `GNAQ`, `GNA11`, `ARRB1`, `ARRB2` | Antipsychotic-relevant receptors and proximal G-protein or arrestin signaling components used to compare DEG enrichment with the pharmacological receptor-response model. |
| Calcium signaling and neuronal excitability | `CACNA1A`, `CACNA1B`, `CACNA1C`, `CACNA1D`, `CACNA1E`, `CACNA1G`, `CACNA1H`, `CACNA1I`, `CAMK2A`, `CAMK2B`, `CALM1`, `CALM2`, `CALM3`, `PPP3CA`, `PPP3CB`, `PPP3CC`, `RYR1`, `RYR2`, `RYR3`, `ITPR1`, `ITPR2`, `ITPR3` | Voltage-gated calcium channels, calmodulin/calcineurin components, and intracellular calcium-release genes used to assess Ca-associated signaling and excitability-related DEG patterns. |
| Immediate-early and activity-regulated transcription | `FOS`, `FOSB`, `JUN`, `JUNB`, `JUND`, `EGR1`, `EGR2`, `EGR3`, `EGR4`, `NPAS4`, `ARC`, `IER2`, `IER3`, `NR4A1`, `NR4A2`, `NR4A3`, `BDNF` | Activity-dependent transcriptional regulators and plasticity-associated genes used to evaluate whether clozapine-selective cells showed an activated or activity-regulated expression signature. |

For the curated over-representation test, significant DEGs were separated into upregulated and downregulated lists using the same FDR and absolute log2 fold-change thresholds used for DEG interpretation. For each direction and each curated term, a two-by-two contingency table was constructed from: significant DEG genes inside the curated term, significant DEG genes outside the term, nonsignificant background genes inside the term, and nonsignificant background genes outside the term. Enrichment P values were calculated using a two-sided Fisher's exact test. Odds ratios were reported with Haldane-Anscombe correction when needed for zero cells, and P values across curated terms and DEG directions were adjusted using the Benjamini-Hochberg FDR procedure. Curated terms were considered enriched at FDR < 0.05; nominal P values, adjusted P values, odds ratios, overlap counts, term sizes, and overlapping gene symbols were reported for transparency.

As a sensitivity analysis, ranked enrichment was performed for the same curated terms using all genes in the DEG table ranked by the signed differential-expression statistic. The primary ranking metric was the Scanpy Wilcoxon score multiplied by the sign of the log fold change when both fields were available; otherwise, genes were ranked by log2 fold change. Enrichment was evaluated by comparing the observed concentration of curated-term genes near the top or bottom of the ranked DEG list with label-permuted or preranked enrichment statistics, depending on the implementation available in the analysis environment. Curated enrichment results were interpreted only when supported by both the direction of the DEG overlaps and the underlying gene-level effects, to avoid overinterpreting small gene sets driven by one or two highly ranked genes.

## Statistical analysis

Unless otherwise specified, statistical tests were two-sided. Differential-expression testing used Wilcoxon rank-sum tests as implemented in Scanpy, with Benjamini-Hochberg adjustment for multiple comparisons. Receptor-level and marker-level comparisons between clozapine-selective and nonselective cells used Mann-Whitney U tests with Benjamini-Hochberg correction. Cell-type and cluster enrichment analyses used Fisher's exact tests with Benjamini-Hochberg correction across categories. Inhibitor-cocktail cAMP-pattern comparisons were summarized descriptively and tested using Welch's two-sample t-test and the Mann-Whitney U test. For visualization, adjusted P values equal to zero were replaced by the smallest positive adjusted P value in the corresponding result table or by a numerical floor to permit log-scale plotting. Data are reported as cell counts, means, medians, standard deviations, odds ratios with 95% confidence intervals, log2 fold changes, and FDR-adjusted P values as appropriate.
