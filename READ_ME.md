# antipsychotic pipeline (cross-region / cross-species)

はい、その設計が良いです。  
**`preprocess.py` で raw_data → 統一フォーマット(cell x gene) を作成し、下流は同一 pipeline に流す**構成にしました。

## 追加したもの

- `preprocess.py` : 複数データセットを読み込み、遺伝子名を大文字へ統一、metadata付与後に結合して標準化 matrix を出力
- `preprocess_config.json` : raw 入力群と出力先を定義
- `run_antipsychotic_pipeline.py` : 標準化済み matrix を前提に、薬剤反応予測〜clz selective 判定〜任意 cAMP 比較を実施
- `antipsychotic_pipeline_config.json` : 下流解析の設定
- `run_camp_pattern_analysis.py` / `camp_analysis_config.json` : cAMP 比較単体実行

## 推奨ワークフロー

### 1) 前処理（raw → standardized matrix）

```bash
python preprocess.py --config preprocess_config.json
```

出力:
- `outputs/preprocess/standardized_matrix.h5ad`
- `outputs/preprocess/standardized_obs.csv`
- `outputs/preprocess/standardized_var.csv`

### 2) 下流パイプライン（統一処理）

`antipsychotic_pipeline_config.json` の `input.adata_h5ad` を
`outputs/preprocess/standardized_matrix.h5ad` に向けて実行:

```bash
python run_antipsychotic_pipeline.py --config antipsychotic_pipeline_config.json
```

### GPU を明示指定する

`antipsychotic_pipeline_config.json` の `preprocess.gpu_device` を設定すると、
実行時に `CUDA_VISIBLE_DEVICES` を自動設定します。

例:

```json
"preprocess": {
  "is_gpu": true,
  "gpu_device": "1"
}
```

※ `nvidia-smi` で空いている GPU 番号を選んでください。

## preprocess 設計ポイント

- `file_type` は `h5ad`, `10x_mtx`, `10x_h5`, `csv` をサポート
- 遺伝子名は uppercase 化 + unique 化
- `obs` に `dataset`, `species`, `region` を付与
- 複数 dataset は `anndata.concat(..., join="outer")` で結合

## 将来拡張（推奨）

- species 別 ortholog マッピング（human/mouse統合時）
- 低品質 cell / doublet QC の前処理段階での明示化
- `preprocess` と `pipeline` の schema バリデーション追加

## Notebookでclozapine-responsive cellのFigure案を対話的に出す

`clz_responsive_figure_analysis.py` を使うと、既に `is_clz_selective` を持つ `adata` から、Figure構成案と主要な補助テーブル/プロットをNotebook上で対話的に確認できます。

```python
import clz_responsive_figure_analysis as clzfig

result = clzfig.run_clz_responsive_figure_story(
    adata,
    group_col="is_clz_selective",
    cell_type_col="cell_type",  # データセットに合わせて annotation / subclass / label などに変更
    output_dir="outputs/clz_responsive_figures",
)
```

出力される主なファイル:

- `figure_plan.csv`: 推奨Figure panel構成
- `umap_clz_responsive.png`: clozapine-responsive cells / response score のUMAP
- `cell_type_enrichment.csv` / `cell_type_enrichment_odds_ratio.png`: cell type enrichment
- `layer_marker_dotplot_table.csv` / `layer_marker_dotplot.png`: excitatory/layer marker dot plot用データ
- `receptor_signature_summary.csv`: HTR1A/CHRM2/HRH3 と cAMP関連 signature
- `functional_gene_set_summary.csv`: synapse/plasticity/cAMP gene set score
