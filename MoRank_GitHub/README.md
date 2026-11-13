# MoRank

Experiment code for *MoRank: Moment and multiscale rank features for gene expression fluctuations* by Shunjie Chen, Baokuan Zhang, Yunxian Wang and Li Dai.

MoRank ranks genes by changes in their expression distributions. This repository contains the simulation and lung adenocarcinoma analysis scripts, together with the frozen XGBoost models used for evaluation. The scripts and weights correspond to the results checked in the manuscript reproduction run. The final manuscript revision uses the same experiments.

## Setup

Use R 4.5 and Python 3.10 or later. The reference run uses R 4.5.3, xgboost 3.2.1.1, matrixStats 1.5.0, data.table 1.18.4, DESeq2 1.50.2 and edgeR 4.8.2. Keep these versions when checking numerical equality. The common-gene benchmark additionally uses waddR and the bundled MMDAgg implementation.

```sh
Rscript simulation/install_dependencies.R
python -m pip install -r requirements.txt
python run.py smoke
```

If Rscript is outside PATH, pass its location with `python run.py --rscript PATH_TO_RSCRIPT smoke`. All experiment outputs and logs default to `runs/`. Choose another location with `--output` and adjust concurrency with `--workers` before the command. The runner verifies frozen file checksums before execution.

## Simulation experiments

| Command | Experiment |
| --- | --- |
| `python run.py simulation primary` | 240 negative-binomial datasets, the reference learner and the original comparators |
| `python run.py simulation ablations` | Frozen ablation models on the 240 primary datasets. Run `primary` first |
| `python run.py simulation stress` | 80 parameter-stress datasets and the final-rule replay. Run `primary` first |
| `python run.py simulation robustness` | First 192-dataset held-out count-generator evaluation |
| `python run.py simulation replication` | Independent 192-dataset replication plus 48 smaller negative-binomial datasets |
| `python run.py simulation null` | 52 global-null datasets with 99 whole-sample label permutations per dataset |
| `python run.py simulation common` | waddR and MMDAgg comparison on the same 400 genes per dataset across 24 datasets |

Run every stage with `python run.py simulation all`. Full runs can take hours. For model development, `python run.py train` regenerates fitting and validation features, refits the reference and ablation models, and fits the mixed models, specialists and routing classifier under `runs/training/`. Released model files are loaded from the repository, and training outputs go to the run directory.

Run `python run.py simulation routing` to compare the original routing classifier against the two alternative classifiers on the separate profile-validation set.

The source IDs S6–S9 in the original simulation implementation correspond to scenarios S5–S8 in manuscript Table 2. The internal model name `WSKR_A` refers to the reference MoRank learner. The final rule also uses the mixed-generator model, negative-binomial specialist and routing classifier.

## Lung adenocarcinoma experiments

```sh
python run.py real all --input-dir PATH_TO_FROZEN_INPUTS --tcga-rdata PATH_TO_TCGA_LUAD_LIST_RDATA
```

The real-data workflow uses the analysis-specific count snapshots, mapping tables, paired metadata and reference networks described in [the input manifest](real_data/INPUTS.md). Obtain the public cohorts from TCGA-LUAD, GSE81089 and GSE131907. The exact frozen and derived inputs listed in the manifest are needed to replay the paper results. The repository contains the code and frozen model weights. Obtain the study input snapshots separately.

The stages are `discovery`, `bulk`, `singlecell`, `string` and `context`. Use the same output directory to pass intermediate results between stages. The optional `tcga_deseq2`, `tcga_edger`, `tcga_aux`, `bulk-comparators` and `singlecell-input` stages reproduce additional comparator or input checks. The `ko` stage reruns manifold alignment when the archived reference alignment function is provided with `--alignment-r`. The main `context` stage recomputes the four matched-reference endpoints from the frozen perturbation response vectors.

## Third-party code

The MMDAgg implementation in `simulation/common/vendor/` retains its accompanying MIT license. The benchmark calls DESeq2, edgeR and waddR as external R packages.
