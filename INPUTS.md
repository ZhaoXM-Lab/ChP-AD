# Analysis inputs

Sample tables, cluster annotations and external reference paths define the analysis inputs. Sequencing data and reference datasets are stored outside version control under `data/`.

## Metadata tables

| Table | Columns |
| --- | --- |
| `RNA/00_reference/samples.tsv` | sample, species, disease, source, include, matrix_dir, fastq_dir, fastq_sample, transcriptome |
| `ATAC/00_reference/samples.tsv` | sample, disease, include, outs_dir, fastq_dir, fastq_sample, reference |
| `RNA/00_reference/cluster_annotations.tsv` | species, cluster, cell_type |
| `RNA/00_reference/cpec_annotations.tsv` | species, cluster, cell_type, include |
| `ATAC/00_reference/epi_validation.tsv` | cell_type, assay, file |
| `ATAC/00_reference/mcl_inflation.tsv` | cell_type, inflation |
| `ATAC/00_reference/gwas.tsv` | trait, sumstats |
| `ATAC/00_reference/loci.tsv` | gene, cell_type, lead_snp |

`species` is `human` or `mouse`; `disease` is `AD` or `WT`; `include` is `1` or `0`. Sample identifiers are unique and contain no path separators. `matrix_dir` points to a 10X gene-expression matrix. Human nuclear calling uses raw-feature matrices. ATAC `outs_dir` contains `filtered_peak_bc_matrix.h5`, `singlecell.csv`, `fragments.tsv.gz` and its index.

Broad annotations use `cell_type`; `cell_subtype` distinguishes CPEC-1 and CPEC-2. Cell names and marker panels are defined in `RNA/00_reference/cell_types.json`. Cluster identifiers map to the corresponding species and clustering run. The CPEC include column controls subcluster retention.

## Genomic references

| Reference | Location or format |
| --- | --- |
| Gene annotation | `data/reference/genes.gtf.gz`, GRCh38 |
| TF motifs | `data/reference/Vierstra2.1_Motif_individual.rds`, named PWM list |
| TF families | `data/reference/tf_families.tsv`, columns TF and family |
| PLINK panel | `data/reference/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.{1..22}.{bed,bim,fam}` |
| HapMap3 SNP list | `data/reference/ldsc/hm3.snplist` |
| Baseline annotations | `data/reference/ldsc/baselineLD_v2.1/` |
| Regression weights | `data/reference/ldsc/weights/weights.hm3_noMHC.{chr}` |
| Allele frequencies | `data/reference/ldsc/frequencies/1000G.EUR.QC.{chr}.frq` |
| GWAS summary statistics | Paths in `ATAC/00_reference/gwas.tsv` |
| Fine-mapping summary statistics | `data/reference/gwas/AD.finemap.tsv`, columns SNP, A1, A2, Z and N |
| Coordinate conversion | `data/reference/hg38ToHg19.over.chain` and `hg19ToHg38.over.chain` |

Genomic intervals use UCSC chromosome names. BED coordinates are zero-based and half-open; GRanges coordinates are one-based and closed. SNP rows follow the reference BIM order. `LDSC_BUILD` and `FINEMAP_BUILD` identify the reference-panel coordinate system.

`RNA/05_annotation/05_mouse_orthologs.py` queries Ensembl BioMart for human-mouse orthologues and saves the table under `data/RNA/output/annotation/`.

## EPI validation

The validation table links each cell type to an assay and input file. eQTL tables contain Chromosome, Start, End and gold_gene_id. pcHiC tables contain chr1_hg38, start1_hg38, end1_hg38, chr2_hg38, start2_hg38 and end2_hg38. Both use GRCh38 one-based closed coordinates; the first pcHiC anchor is the promoter and the second is the enhancer.

## Expression matrices

RNA AnnData objects store raw counts in `layers['counts']` and log-normalized full-gene expression in `X` and `raw`. CellOracle and hormone2cell receive raw-count inputs. Saved analysis objects and figure tables are written by the corresponding numbered scripts.
