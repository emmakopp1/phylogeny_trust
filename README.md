# Phylogenetic Reconstruction Accuracy in Linguistic Evolution

> Evaluating the performance of MCC, consensus tree, and HIPSTR methods for phylogenetic reconstruction in linguistic evolution, combining simulation studies with real-world analyses of Indo-European and Sino-Tibetan language families.

---

## Table of Contents

- [Overview](#overview)
- [Project Structure](#project-structure)
- [Methodology](#methodology)
- [Quick Start](#quick-start)
- [Script Reference](#script-reference)
- [Dependencies](#dependencies)
- [Data Availability](#data-availability)
- [Citation](#citation)

---

## Overview

This repository provides the full reproducible pipeline for a study quantifying the temporal limits of phylogenetic inference in historical linguistics. Using 850 simulations with root ages from 1 to 17 ka BP alongside empirical analyses of Sino-Tibetan and Indo-European, it shows that the probability of correctly recovering tree topology drops sharply beyond 8 ka BP and becomes indistinguishable from chance at 11–12 ka BP.

The study combines:
- **Simulation experiments** — birth-death trees at 17 ages, with 1500–12 000 binary traits
- **Real-data analyses** — Bayesian phylogenetic inference on Indo-European and Sino-Tibetan using BEAST 2
- **Ancestral state reconstruction** — tracing the evolution of semantic meanings (cognate presence/absence) across posterior tree samples

---

## Project Structure

```
phylogeny_trust/
├── data/
│   ├── real/
│   │   ├── iecor_ctmc-strict-M1/          # Indo-European BEAST results
│   │   └── st_ctmc-strict-fbd-uni/        # Sino-Tibetan BEAST results
│   ├── simulated-2025-07-22-1500/         # Simulation — 1 500 traits
│   ├── simulated-2025-07-22-6000/         # Simulation — 6 000 traits
│   ├── simulated-2025-07-22-12000/        # Simulation — 12 000 traits
│   ├── simulated-2025-07-28/              # Main simulation — 3 000 traits
│   ├── beast-data-sim.xml                 # BEAST template (sequence simulation)
│   └── ctmc-strict-bd-template.xml        # BEAUti template (CTMC birth-death)
│
├── src/
│   ├── simulations/                       # Steps 01–04 (see Script Reference)
│   ├── ancestral_reconstruction/          # Steps 00–02 (see Script Reference)
│   ├── shared_cognate.R                   # Theoretical shared-cognate probability
│   └── plots.R                            # All publication figures
│
└── output/
    ├── results/                           # Processed datasets and statistics
```

---

## Methodology

### Simulation Study

50-taxon birth-death trees (calibrated on Sino-Tibetan clades) are generated at **17 ages** and evaluated across **four trait-count conditions** (1 500, 3 000, 6 000, 12 000). Binary traits evolve under a CTMC model; posterior tree samples are summarised with three methods and compared to the true tree.

| Dimension | Values |
|---|---|
| Tree size | 50 taxa |
| Ages evaluated | 17 |
| Trait counts | 1 500 · 3 000 · 6 000 · 12 000 |
| Summary methods | MCC · 50% MRC · HIPSTR |
| Evaluation metrics | Node accuracy · First split · Marginal prob. · Robinson-Foulds |

### Real Data Analysis

| Family | Dataset | Tree prior |
|---|---|---|
| Indo-European | IECor (Heggarty et al. 2023) | Birth-death |
| Sino-Tibetan | ST (Sagart et al. 2019) | Fossilised birth-death |

---

## Quick Start

> **Prerequisites**: R ≥ 4.0, BEAST 2, and the R packages listed in [Dependencies](#dependencies).

### Reproduce the simulation analyses

```bash
# 1. Generate trees and BEAST XMLs (optional — pre-computed data provided)
Rscript src/simulations/01_tree_simulation_main.R

# 2. Run BEAST 2 on the generated XMLs

# 3. Compute summary trees
#    ⚠ Edit the `simulation_folder` variable in each script before running
Rscript src/simulations/02_compute_consensus.R
Rscript src/simulations/02_compute_mcc.R
bash   src/simulations/03_compute_hipstr.sh

# 4. Run analyses
Rscript src/simulations/03_compute_rf.R
Rscript src/simulations/03_marginal_prob_first_split_ic.R   # ⚠ Costly — see note below
Rscript src/simulations/03_marginal_prob_first_split_resumed.R
Rscript src/simulations/03_reconstruction_difficulty.R
Rscript src/simulations/03_resume_to_true_TF.R
Rscript src/simulations/03_true_false_uncertain.R

# 5. Consolidate and model
Rscript src/simulations/04_simulation_analyses_main.R
Rscript src/simulations/04_simulation_analyses_trait_influence.R
```

> **Note** — `03_marginal_prob_first_split_ic.R` is computationally intensive (originally run on a cluster). For testing, reduce the number of sampled posterior files to ~50 in the script, or reuse the existing output CSV.

> **Note** — Scripts in steps 3–5 rely on **hardcoded paths** pointing to a specific `data/simulated-*` folder. Edit the `simulation_folder` / `cluster_directory` variable at the top of each script before running.

### Reproduce the ancestral reconstruction

```bash
Rscript src/ancestral_reconstruction/00_compute_meaning_set.R
Rscript src/ancestral_reconstruction/00_tracelogs.R
Rscript src/ancestral_reconstruction/01_ancestral_state_reconstruction_ie.R   # ⚠ Costly
Rscript src/ancestral_reconstruction/01_ancestral_state_reconstruction_st.R   # ⚠ Costly
Rscript src/ancestral_reconstruction/02_post_process.R
```

### Shared cognate analysis

```bash
Rscript src/shared_cognate.R
```

### Generate all figures

```bash
Rscript src/plots.R
```

---

## Script Reference

### `src/simulations/`

| Script | Purpose | Key inputs | Key outputs |
|---|---|---|---|
| `01_tree_simulation_*.R` | Simulate birth-death trees and BEAST XMLs | ST trees, BEAST templates | `tree-sim-{n}-{k}.tree`, BEAST XMLs |
| `02_compute_consensus.R` | 50% majority-rule consensus trees | Posterior `.trees` | `consensus-{age}.tree` |
| `02_compute_mcc.R` | MCC trees | Posterior `.trees` | `mcc-{age}.tree` |
| `03_compute_hipstr.sh` | HIPSTR summary trees via TreeAnnotator | Posterior `.trees` | `hipstr-{age}.tree` |
| `03_compute_rf.R` | Robinson-Foulds distances (true vs. posterior) | `.trees`, `tree-sim-*` | `rf_values[_{N}].csv` |
| `03_marginal_prob_first_split_ic.R` | Posterior prob. of recovering the true first bipartition | Posterior `.trees` | `marginal_prob_first_split_ic.csv` |
| `03_marginal_prob_first_split_resumed.R` | First-split clade support in MCC/consensus/HIPSTR | Summary trees | `marginal_prob_first_split_mcc_consensus_hipstr.csv` |
| `03_reconstruction_difficulty.R` | Root age and first-split age from summary trees | Summary trees | `first_split_age_{mcc,cs,hipstr}.csv` |
| `03_resume_to_true_TF.R` | Strict true/false monophyly check per node | True + summary trees | `resume_to_true_TF[_{N}].csv` |
| `03_true_false_uncertain.R` | 3-way node classification (true/false/uncertain) | True + consensus trees | `true_false_uncertain_nodes[_{N}].csv` |
| `04_simulation_analyses_main.R` | HDIs, proportions, logistic regressions (3 000-trait) | Multiple CSVs | Model RDS, summary CSVs |
| `04_simulation_analyses_trait_influence.R` | Accuracy comparison across trait counts | CSVs from step 03 | Long-format comparison CSVs |

### `src/ancestral_reconstruction/`

| Script | Purpose | Key outputs |
|---|---|---|
| `00_compute_meaning_set.R` | Extract character index ranges from NEXUS blocks | `meanings_sets_{ie,st}.csv` |
| `00_tracelogs.R` | Post-burnin trace log summaries (clock rate, tree height) | `ntipschars.csv`, `tracelog_summary.csv` |
| `01_ancestral_state_reconstruction_{ie,st}.R` | Mk-model ancestral state reconstruction across 200 trees | `ancestral_reconstruction_{ie,st}.csv` |
| `02_post_process.R` | Per-meaning summary statistics and per-trait IE breakdown | `ancestral_reconstruction_summary_*.csv` |

---

## Dependencies

### R packages

| Category | Packages |
|---|---|
| Phylogenetics | `ape` · `phangorn` · `phytools` · `TreeSim` · `castor` · `adephylo` · `TreeTools` · `treeio` |
| BEAST interface | `beastier` · `tracerer` |
| Data wrangling | `tidyverse` · `dplyr` · `tidyr` · `purrr` · `reshape2` · `stringr` · `readr` · `tibble` · `magrittr` |
| Modelling | `broom` · `stats` · `Matrix` · `parallel` |
| Visualisation | `ggplot2` · `patchwork` · `ggeffects` · `dotwhisker` · `pheatmap` · `gridExtra` |
| Utilities | `here` · `xml2` |

### External software

| Tool | Version | Purpose |
|---|---|---|
| [BEAST 2](https://www.beast2.org/) | 2.6.7 | Bayesian phylogenetic inference and HIPSTR (via TreeAnnotator) |

---

## Data Availability

- **Real data**: phylogenetic analyses of Indo-European and Sino-Tibetan language families (sources in [Citation](#citation))
- **Processed results**: all intermediate CSVs in `output/results/`
- **Figures**: publication-ready PDFs in `output/figs/`

---

## Citation

If you use this code or data, please cite:

> Sagart, L., Jacques, G., Lai, Y., Ryder, R. J., Thouzeau, V., Greenhill, S. J., & List, J. M. (2019). Dated language phylogenies shed light on the ancestry of Sino-Tibetan. *PNAS*, 116(21), 10317–10322.

> Heggarty, P., Anderson, C., Scarborough, M., King, B., Bouckaert, R., Jocz, L., … & others. (2023). Language trees with sampled ancestors support a hybrid model for the origin of Indo-European languages. *Science*, 381(6656), eabg0818.
