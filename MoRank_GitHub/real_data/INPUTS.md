# Input files for the lung adenocarcinoma experiments

Pass the root of the frozen input directory with `--input-dir`. The path names below match the analysis snapshots used in the manuscript reproduction. Supply the TCGA working RData with `--tcga-rdata` or place it at `TCGA-LUAD_list.RData` in the input root. It contains `list_data$puried_data` and the aligned sample group labels. The analysis checks the 541 tumour and 59 normal samples and the final 29,196-gene working matrix.

| Stage | Files relative to the input root |
| --- | --- |
| Candidate set and identifiers | `round33_model_optimization/results/TCGA_exploratory_ranks.csv`, `round9/results/gene_mapping_audit.csv` |
| Discovery and subsampling | `round7/results/all_gene_complement_audit.csv`, `round13/data/TCGA_aux_Wilcoxon_rawP.csv`, `round6/raw/real/real_plans_first12.rds`, `round6/raw/real/real_scores_rep1.rds` through `real_scores_rep12.rds` |
| GSE81089 reranking and paired bulk analysis | `round13/sources/GSE81089_readcounts_featurecounts.tsv.gz`, `round13/data/external_frozen_features.rds`, `round13/data/ensembl_symbol_map.csv`, `round13/results/external_full_cohort_ranking.csv`, `round13/results/external_paired_all_tested.csv.gz` |
| GSE131907 epithelial expression and depth | `round10/data/raw_gene_mapping.csv`, `round10/data/classified_epithelial_counts.rds`, `round10/data/pseudobulk_full.rds`, `round10/data/primary_lung_input_cache.rds`, `round10/data/primary_lung_cell_metadata.csv.gz`, `round10/results/pseudobulk_all_tested_genes.csv` |
| STRING interactions | Under `round20_network/sources/`, `9606.protein.info.v12.5.txt.gz`, `9606.protein.links.v12.5.txt.gz` and `9606.protein.physical.links.v12.5.txt.gz` |
| Pathway enrichment and activity | `round9/results/mapped_universe.csv`, `round9/results/term_gene.csv`, `round9/results/database_backgrounds.csv`, `round20_network/data/full_term_gene.csv`, `round13/data/external_analysis_manifest.rds`, plus the bulk counts, pseudobulk counts and mappings above |
| IRS2 modules and perturbation | `round10/results/module_IRS2_paired_changes.csv`, `round10/results/module_tests_all6.csv`, `round10/results/reference_KO/perturbations_all102.rds`, `round10/data/IRS2_reference_subnetwork.rds` |

The optional `tcga_aux` check also reads `round3/inputs/gene_scores.csv`. The `singlecell-input` check reads `round10/sources/GSE131907_Lung_Cancer_raw_UMI_matrix.rds` and `round10/sources/cell_annotation.tsv`, and compares the reconstructed count subsets with the frozen files. For `ko`, provide the exact archived `scTenifoldNet_manifoldAlignment.R` function file via `--alignment-r`.

Public sources are [TCGA-LUAD](https://portal.gdc.cancer.gov/projects/TCGA-LUAD), [GSE81089](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81089), [GSE131907](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131907) and [STRING](https://string-db.org/). Public source downloads and the study-specific frozen intermediates serve different roles. Use the exact inputs for numerical verification against the paper.
