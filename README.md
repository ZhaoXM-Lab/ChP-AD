# ChP-AD

Single-nucleus transcriptomic and chromatin accessibility analysis of the human choroid plexus in Alzheimer disease, with mouse transcriptomic comparisons.

## Analysis workflows

`RNA/` contains sequencing preprocessing, species-specific quality control, Harmony integration, cell-type annotation, CPEC subclustering, differential expression, functional enrichment, cell-cell communication, hormone signaling and drug-target analysis.

`ATAC/` contains sequencing preprocessing, Signac and ArchR quality control, LSI/Harmony integration, RNA label transfer, cCRE analysis, enhancer-promoter interactions, regulatory networks, MCL clustering, CellOracle perturbations, stratified LD score regression and SuSiE fine-mapping.

Each workflow uses numbered directories and scripts. Reference tables are under `00_reference/`; figure scripts are under `11_figures/`. LDSC is under `ATAC/09_ldsc/`.

## Environment

The analysis scripts run in a Linux environment. Python dependencies are listed in `conda-environment.yml`, R dependencies in `R-packages.tsv`, and the CellOracle environment in `ATAC/00_reference/celloracle-environment.yml`. External command-line tools include Cell Ranger 9.0.0, Cell Ranger ATAC 2.2.0, MACS2, BEDTools, MCL and LDSC.

```bash
conda env create -f conda-environment.yml
conda activate chp-ad-rna
export PROJ_DIR="$PWD"
export THREADS=16
```

`PROJ_DIR` is the repository root. Analysis inputs and outputs are stored under `data/`; figures are written under `figures/`. Table paths are relative to `PROJ_DIR`. Input schemas and reference files are described in [INPUTS.md](INPUTS.md).

## RNA

```bash
export SPECIES=human
python RNA/01_data_preprocessing/01_cellranger.py
python RNA/01_data_preprocessing/02_load_data.py
python RNA/02_qc/01_quality_control.py
python RNA/03_integration/01_integrate_harmony.py
python RNA/04_umap/01_run_umap.py
python RNA/05_annotation/01_annotate_cell_types.py
python RNA/05_annotation/02_cpec_subclustering.py
python RNA/05_annotation/03_export_seurat_input.py
Rscript RNA/05_annotation/04_build_seurat.R
```

Preprocessing and QC use the species column in the RNA sample manifest. Integration and annotation select the species with `SPECIES=human` or `SPECIES=mouse`. CPEC subclustering uses `CPEC_RESOLUTION`, with a default of 0.1.

| QC criterion | Human | Mouse |
| --- | --- | --- |
| Detected genes per nucleus | >1000 and <8000 | >500 and <4000 |
| UMIs per nucleus | >1000 and <30000 | >500 and <10000 |
| Mitochondrial percentage | <10% | <10% |
| Scrublet score | <0.1 | <0.15 |
| Minimum cells per gene | 3 | 2 |

Differential expression is under `RNA/06_differential_expression/`. Communication and drug-target analyses use the annotated RNA atlas and differential-expression results.

## ATAC and downstream analyses

Execute `ATAC/00_reference/01_prepare_annotations.R`, followed by ATAC stages 01 through 05 in numeric order. RNA annotation supplies the reference for label transfer. Stage 06 performs peak calling, cCRE analysis, NMF and GREAT enrichment. Functional concordance also uses the RNA enrichment results. Stage 07 infers and validates enhancer-promoter interactions.

Stage 08 builds regulatory networks and scans MCL inflation values. Selected inflation values are stored in `ATAC/00_reference/mcl_inflation.tsv`. CellOracle uses `CELL_TYPE` to select a lineage; example perturbations are TCF4 in CPEC and AHR in macrophages.

For stage 09, set `LDSC_DIR` and `LDSC_BUILD` (`hg19` or `hg38`); `LDSC_PYTHON` selects the LDSC interpreter. Stage 10 uses `FINEMAP_BUILD` (`hg19` or `hg38`). The matching liftOver chain is used for hg19 reference panels.

Run the figure scripts after their corresponding analysis stages.
