# Phylogenetic Reconstruction Accuracy in Linguistic Evolution

## Abstract

This repository contains code and data for evaluating phylogenetic reconstruction accuracy in linguistic evolution. The study examines the performance of Maximum Clade Credibility (MCC), consensus tree, and HIPSTR methods across varying data quantities and linguistic families, combining simulation studies with real-world phylogenetic analyses.

## Project Structure

```
phylogeny_trust/
├── data/                                  # Raw and processed data files
│   ├── real/                              # Real linguistic phylogenetic data
│   │   ├── iecor_ctmc-strict-M1/          # Indo-European (IE) BEAST results
│   │   └── st_ctmc-strict-fbd-uni/        # Sino-Tibetan (ST) BEAST results
│   ├── shared_cognates/                   # (currently empty; created by simulate_phylogenies.R)
│   ├── simulated-2025-07-22-1500/         # Simulation dataset (1500 traits)
│   ├── simulated-2025-07-22-6000/         # Simulation dataset (6000 traits)
│   ├── simulated-2025-07-22-12000/        # Simulation dataset (12000 traits)
│   ├── simulated-2025-07-28/              # Main simulation dataset (3000 traits)
│   ├── beast-data-sim.xml                 # BEAST template (sequence simulation)
│   ├── ctmc-strict-bd-template.xml        # BEAUti template (CTMC birth-death)
│   └── bounds_real_tb_by_sens.csv
├── src/                                   # Source code organized by analysis type
│   ├── simulations/                       # Phylogenetic simulation analyses
│   │   ├── 01_tree_simulation_main.R
│   │   ├── 01_tree_simulation_1500.R
│   │   ├── 01_tree_simulation_6000.R
│   │   ├── 01_tree_simulation_12000.R
│   │   ├── 02_compute_consensus.R
│   │   ├── 02_compute_mcc.R
│   │   ├── 03_compute_number_of_nodes_resumed.R
│   │   ├── 03_compute_rf.R                # currently broken, see note below
│   │   ├── 03_marginal_prob_first_split_ic.R
│   │   ├── 03_marginal_prob_first_split_resumed.R
│   │   ├── 03_reconstruction_difficulty.R
│   │   ├── 03_resume_to_true_TF.R
│   │   ├── 03_true_false_uncertain.R
│   │   ├── 04_simulation_analyses_main.R
│   │   └── 04_simulation_analyses_trait_influence.R
│   ├── shared_cognates/                   # Shared cognate analysis
│   │   ├── shared_cognate_theorique.R     # Theoretical (closed-form) estimate
│   │   └── simulate_phylogenies.R         # Empirical estimate via trait simulation
│   ├── ancestral_reconstruction/          # Ancestral state reconstruction
│   │   ├── 00_compute_meaning_set.R
│   │   ├── 00_tracelogs.R
│   │   ├── 01_ancestral_state_reconstruction_ie.R
│   │   ├── 01_ancestral_state_reconstruction_st.R
│   │   ├── 02_tree_pruned_recuperation.R
│   │   └── 03_post_process.R
│   └── plots.R                            # All final figures for the whole project
├── output/                                # Generated results and figures
│   ├── results/                           # Processed datasets and statistics
│   ├── figs/                              # Publication-ready figures
│   └── trees/                             # Pruned tree objects (ancestral reconstruction)
└── archive/                                # Superseded/legacy scripts (not on GitHub)
```

## Methodology

### Simulation Study

The simulation study evaluates phylogenetic reconstruction accuracy using birth-death trees with varying evolutionary ages and trait counts:

- **Tree Simulation**: 50-taxon birth-death trees (calibrated on Sino-Tibetan clades) scaled to 17 ages
- **Trait Evolution**: Binary trait evolution under CTMC models
- **Data Sizes**: 1500, 3000 (main), 6000, and 12000 traits
- **Reconstruction Methods**: MCC trees, 50% majority-rule consensus trees, and HIPSTR summary trees
- **Evaluation Metrics**: Node accuracy, first split identification, marginal probabilities, Robinson-Foulds distance

### Real Data Analysis

Phylogenetic reconstruction of multiple language families using:
- **Families**: Indo-European (data outputs files taken from [2]), Sino-Tibetan (input data taken from [1])
- **Methods**: BEAST with birth-death/fossilized birth-death tree priors
- **Ancestral Reconstruction**: Semantic meaning (cognate presence/absence) evolution analysis

## Script Documentation

### 1. Simulation Analysis (`src/simulations/`)

**`01_tree_simulation_main.R`, `01_tree_simulation_1500.R`, `01_tree_simulation_6000.R`, `01_tree_simulation_12000.R`**
- **Purpose**: Simulate 50-taxon birth-death trees at 17 ages (main) or a single age (variants), then generate and run BEAST XMLs to simulate trait sequences (3000/1500/6000/12000 traits respectively)
- **Input**: `data/real/st_ctmc-strict-fbd-uni/st_ctmc-strict-fbd-uniform.trees`, `data/beast-data-sim.xml`, `data/ctmc-strict-bd-template.xml`
- **Output**: `data/simulated-{date}[-{1500,6000,12000}]/beast-data-sim-{n}/.../`:
  - `tree-sim-{n}-{k}.tree` — true phylogenetic trees
  - `beast-data-sim-{n}-{k}.xml`, `beast-simulated-seq-{n}-{k}.xml`, `ctmc-strict-bd-{n}-{k}.xml` — BEAST config/output files

**`02_compute_consensus.R` / `02_compute_mcc.R`**
- **Purpose**: Compute 50% majority-rule consensus trees / MCC trees from posterior `.trees` distributions
- **Input**: `*.trees` files under a **hardcoded** `simulation_folder` variable (must be edited by hand to point at the desired `data/simulated-*` dataset — no CLI parameter)
- **Output**: `consensus-{age}.tree` / `mcc-{age}.tree`, written alongside each input `.trees` file

**`03_compute_number_of_nodes_resumed.R`**
- **Purpose**: Counts internal nodes in MCC and consensus trees
- **Input**: `data/simulated-*/{tree-sim,consensus,mcc}-*.tree` (hardcoded dataset path, edited by hand per trait count)
- **Output**: `output/results/number_nodes_mcc_cs[_{N_traits}].csv`

**`03_compute_rf.R`**
- **Purpose**: Intended to compute Robinson-Foulds distance summary statistics between true trees and posterior tree samples
- **Status**: **Currently broken** — contains invalid `data.frame()` syntax and references an undefined `cl` object. Needs a fix before it can be run; `plots.R` still expects its output (`rf_values*.csv`) to exist in `output/results/`.

**`03_marginal_prob_first_split_ic.R`**
- **Purpose**: Computes the probability (with 95% CI) that posterior tree samples recover the true first bipartition (outgroup split); parallelized, computationally costly
- **Input**: `data/simulated-2025-07-28/{tree-sim,*.trees}` (hardcoded)
- **Output**: `data/simulated-2025-07-28/marginal_prob_first_split_ic.csv` (written into the data folder, not `output/results/`)

**`03_marginal_prob_first_split_resumed.R`**
- **Purpose**: Identifies the first split in MCC, consensus, and HIPSTR summary trees and extracts posterior clade support for it
- **Input**: `data/simulated-*/{mcc,consensus,hipstr}-*` (hardcoded, edited per trait count)
- **Output**: `output/results/marginal_prob_first_split_mcc_consensus_hipstr.csv`

**`03_reconstruction_difficulty.R`**
- **Purpose**: Extracts root age and first-split age from MCC and consensus trees
- **Input**: `data/simulated-2025-07-28/{mcc,consensus}-*`
- **Output**: `output/results/first_split_age_mcc.csv`, `output/results/first_split_age_cs.csv`

**`03_resume_to_true_TF.R`**
- **Purpose**: Strict true/false check of whether each node in consensus/MCC/HIPSTR trees is monophyletic relative to the true tree
- **Input**: `data/simulated-*/{tree-sim,consensus,mcc,hipstr}-*` (hardcoded, edited per trait count)
- **Output**: `output/results/resume_to_true_TF[_{N_traits}].csv`

**`03_true_false_uncertain.R`**
- **Purpose**: 3-way classification (true/false/uncertain) of true-tree nodes as reconstructed in the consensus tree
- **Input**: `data/simulated-*/{tree-sim,consensus}-*` (hardcoded, edited per trait count)
- **Output**: `output/results/true_false_uncertain_nodes[_{N_traits}].csv`

**`04_simulation_analyses_main.R`**
- **Purpose**: Aggregates all per-node/per-split results for the main (3000-trait) study, computes HDIs and proportions, fits logistic regressions of reconstruction success (MCC/consensus/HIPSTR)
- **Input**: Multiple CSVs from `output/results/` (see script)
- **Output**: `output/results/{marginal_prob_first_split_hdi,prop_true_to_cs,prop_mcc_to_true,prop_cs_to_true,prop_hipstr_to_true,rf_hdi,mcc_reconstruction_proba_first_split,cs_reconstruction_proba_first_split,hipstr_reconstruction_proba_first_split,...}.csv`, `output/results/model_{mcc,cs,hipstr}.rds`

**`04_simulation_analyses_trait_influence.R`**
- **Purpose**: Compares reconstruction accuracy across the 1500/3000/6000/12000-trait studies
- **Input**: `output/results/{true_false_uncertain_nodes,resume_to_true_TF}[_{N_traits}].csv`, `prop_true_to_{cs,mcc}_*`
- **Output**: `output/results/prop_true_to_cs_{1500,3000,6000,12000}.rds`, `output/results/prop_true_to_mcc_{1500,3000,6000,12000}.csv` (⚠ actually written with `saveRDS`, despite the `.csv` extension), `output/results/prop_cs_to_true_{1500,3000,6000,12000}.csv`, `output/results/prop_true_to_cs_data_long.csv`, `output/results/prop_true_to_mcc_data_long.csv`

### 2. Shared Cognate Analysis (`src/shared_cognates/`)

**`shared_cognate_theorique.R`**
- **Purpose**: Computes the theoretical (closed-form, no-homoplasy) probability that a cognate present at the root is retained and shared between the two root-level subgroups, as a function of tree age
- **Input**: `data/real/st_ctmc-strict-fbd-uni/st_ctmc-strict-fbd-uniform.trees`, `output/results/meanings_sets_st.csv`, `output/results/tracelog_summary.csv`
- **Output**: `output/results/shared_cognate_thq_no_homoplasie.csv`

**`simulate_phylogenies.R`**
- **Purpose**: Simulates trait-presence histories (birth-death model) along the real ST tree at 17 scaled ages to empirically estimate the proportion of shared cognates (with/without homoplasy). Despite the name, it does not simulate tree topologies — only trait histories on the fixed real tree.
- **Input**: `data/real/st_ctmc-strict-fbd-uni/st_ctmc-strict-fbd-uniform.log`, `data/simulated-2025-07-28/beast-data-sim-1/.../tree-sim-1-{1..17}.tree`
- **Output**: `output/results/shared_cognates.rds`, `output/results/shared_cognate_tip_pair.rds`

### 3. Ancestral State Reconstruction (`src/ancestral_reconstruction/`)

**`00_compute_meaning_set.R`**
- **Purpose**: Extracts semantic meaning boundaries (character index ranges) from NEXUS character-label blocks
- **Input**: `data/real/iecor_ctmc-strict-M1/iecor.nex`, `data/real/st_ctmc-strict-fbd-uni/st.nex`
- **Output**: `output/results/meanings_sets_ie.csv`, `output/results/meanings_sets_st.csv`

**`00_tracelogs.R`**
- **Purpose**: Parses BEAST trace logs and tree files to compute post-burnin summary statistics (state frequencies, clock/mutation rate, tree height) and taxon/character counts, for IE and ST
- **Input**: BEAST tree/log/NEXUS files from `data/real/`
- **Output**: `output/results/ntipschars.csv`, `output/results/tracelog_summary.csv`

**`01_ancestral_state_reconstruction_ie.R` / `01_ancestral_state_reconstruction_st.R`**
- **Purpose**: Performs ancestral state (cognate presence/absence) reconstruction across 200 thinned posterior trees using a fitted Mk model, recording the maximum reconstruction depth per meaning/tree/trait
- **Input**: BEAST trees + NEXUS alignment for the family, `output/results/meanings_sets_{ie,st}.csv`, `output/results/tracelog_summary.csv`
- **Output**: `output/results/ancestral_reconstruction_ie.csv` / `output/results/ancestral_reconstruction_st.csv`
- **Note**: Computationally costly — historically run on a cluster. For testing, reduce the number of posterior trees sampled in the script; for full results use the full posterior sample.

**`02_tree_pruned_recuperation.R`**
- **Purpose**: Re-derives and saves the pruned tree object corresponding to each row of the ancestral-reconstruction result CSVs, for downstream reuse without recomputation
- **Input**: BEAST trees/NEXUS for ST and IE, `meanings_sets_*.csv`, `tracelog_summary.csv`, `ancestral_reconstruction_{st,ie}.csv`
- **Output**: `output/trees/ancr/tree_pruned_{st,ie}/tree_{meaning}_t{t}_trait{trait}.rds`

**`03_post_process.R`**
- **Purpose**: Post-processes raw ancestral-reconstruction CSVs into per-meaning summary statistics (mean max reconstruction depth, mean presence in an early-diverging outgroup lineage), plus a per-trait IE breakdown with root-language annotation
- **Input**: `output/results/ancestral_reconstruction_st.csv`, `output/results/meanings_sets_ie.csv` (⚠ script currently reads this from `data/real/meanings_sets_ie.csv`, which does not exist — likely a stale path, needs fixing), `output/results/ancestral_reconstruction_ie.csv`, `data/real/iecor_ctmc-strict-M1/{iecor_roots.csv,iecor.nex,...}`
- **Output**: `output/results/ancestral_reconstruction_summary_st.csv`, `output/results/ancestral_reconstruction_summary_ie.csv`, `output/results/ancestral_reconstruction_summary_ie_by_sens.csv`

### 4. Figures (`src/plots.R`)

**`plots.R`** (repo root, not project-specific — covers the entire pipeline)
- **Purpose**: Generates all final publication figures: theoretical/empirical shared-cognate plots, node-reconstruction-accuracy barplots, ancestral-reconstruction-depth-vs-outgroup scatterplots, marginal first-split-probability plots, trait-count-influence plots, illustrative tree comparisons, and first-split-support calibration plots
- **Input**: Most CSV/RDS files in `output/results/`, plus example tree files from `data/simulated-2025-07-28/`
- **Output**: All figures in `output/figs/` (PDF, via `ggsave`/`cairo_pdf`)

> Historical note: earlier versions of this pipeline had separate `shared_cognates.R`/`visualization.R` and an ancestral-reconstruction-specific `03_visualization.R`. These have been superseded by the scripts above; the old versions are kept in `archive/` (not tracked on GitHub) for reference only.

## Key Output Files

### Statistical Results (`output/results/`)

- **Simulation Accuracy**: `marginal_prob_first_split_*.csv`
- **Node Classification**: `true_false_uncertain_nodes*.csv`
- **Tree Comparison**: `resume_to_true_TF*.csv`, `rf_values*.csv`
- **Regression Analysis**: `model_{mcc,cs,hipstr}.rds`
- **Ancestral States**: `ancestral_reconstruction*.csv`

### Figures (`output/figs/`)

- **Reconstruction Accuracy**: `marginal_probability_first_split_hdi.pdf`, `{mcc,cs,hipstr}_reconstruction_proba_first_split.pdf`
- **Node Classification**: `barplot_prop_*.pdf`, `plausible_node*.pdf`
- **Ancestral Evolution**: `ancestral_reconstruction_by_semantic_meaning.pdf`
- **Trait Influence**: `number_of_traits_influence.pdf`
- **Tree Distance**: `rf_hdi.pdf`, `rf_trait_influence.pdf`

## Software Dependencies

- **R (≥4.0)** with packages:
  - `ape`, `phangorn` - Phylogenetic analysis
  - `TreeSim` - Tree simulation
  - `phytools` - Phylogenetic tools
  - `tidyverse`, `dplyr`, `tidyr` - Data manipulation
  - `broom` - Statistical modeling
  - `parallel` - Parallel computing
  - `here` - Path management
  - `ggplot2`, `patchwork` - Data visualization
  - `Matrix` - Sparse and dense matrix classes
  - `castor` - Phylogenetic comparative analysis
  - `adephylo` - Phylogenetic signal analysis
  - `stringr` - String manipulation
  - `reshape2` - Data reshaping
  - `purrr` - Functional programming tools
  - `gridExtra` - Grid graphics utilities
  - `magrittr` - Pipe operators
  - `xml2` - XML parsing
  - `readr` - Data import
  - `tibble` - Modern data frames
  - `stats` - Statistical functions
  - `dotwhisker` - Coefficient plots
  - `ggeffects` - Marginal effects visualization
  - `beastier` - BEAST interface
  - `treeio` - Tree I/O operations
  - `tracerer` - BEAST trace log analysis
  - `TreeTools` - Tree manipulation utilities
  - `pheatmap` - Heatmap visualization

- **External Software**:
  - **BEAST 2** - Bayesian phylogenetic analysis

## Usage

### Simulation Analysis (`src/simulations/`)

#### Simulation (optional)
To simulate new trees and data, follow these steps:

1. **Generate phylogenetic trees**:

```bash
Rscript src/simulations/01_tree_simulation_1500.R
Rscript src/simulations/01_tree_simulation_6000.R
Rscript src/simulations/01_tree_simulation_12000.R
Rscript src/simulations/01_tree_simulation_main.R
```

2. **Run inferences manually using BEAST**

3. **Compute summary trees** (edit the `simulation_folder` variable in each script first to point at the target dataset):

```bash
Rscript src/simulations/02_compute_consensus.R
Rscript src/simulations/02_compute_mcc.R
```

#### Analysis
To reproduce the simulation analyses, follow these steps:

1. **Analyze results**:

**Important**:

- In `03_compute_number_of_nodes_resumed.R`, `03_marginal_prob_first_split_resumed.R`, `03_resume_to_true_TF.R`, and `03_true_false_uncertain.R`, the input/output dataset is selected by hardcoded path — edit the script to point at the desired `data/simulated-*` folder and the matching `output/results/*_{N_traits}.csv` output name.

- `03_marginal_prob_first_split_ic.R` was run on a cluster and is computationally costly. We recommend testing on a small number of files (~50) or reusing its existing output.

- `03_compute_rf.R` is currently broken (invalid syntax, undefined variable) and needs to be fixed before use.

Run the analysis files in this order:

   ```bash
   Rscript src/simulations/03_compute_number_of_nodes_resumed.R
   Rscript src/simulations/03_marginal_prob_first_split_ic.R
   Rscript src/simulations/03_marginal_prob_first_split_resumed.R
   Rscript src/simulations/03_reconstruction_difficulty.R
   Rscript src/simulations/03_resume_to_true_TF.R
   Rscript src/simulations/03_true_false_uncertain.R
   ```

2. **Consolidation and modeling**:
   ```bash
   Rscript src/simulations/04_simulation_analyses_main.R
   Rscript src/simulations/04_simulation_analyses_trait_influence.R
   ```

### Ancestral Reconstruction (`src/ancestral_reconstruction/`)

To perform ancestral reconstruction, execute the files in the following order:

1. **Compute meaning sets**:
   ```bash
   Rscript src/ancestral_reconstruction/00_compute_meaning_set.R
   ```

2. **Extract tracelogs**:
   ```bash
   Rscript src/ancestral_reconstruction/00_tracelogs.R
   ```

3. **Ancestral reconstruction**:
   ```bash
   Rscript src/ancestral_reconstruction/01_ancestral_state_reconstruction_ie.R
   Rscript src/ancestral_reconstruction/01_ancestral_state_reconstruction_st.R
   ```
   **Note**: Computationally costly, historically run on a cluster.

4. **Recover pruned trees**:
   ```bash
   Rscript src/ancestral_reconstruction/02_tree_pruned_recuperation.R
   ```

5. **Post-processing**:
   ```bash
   Rscript src/ancestral_reconstruction/03_post_process.R
   ```
   **Note**: this script currently reads `meanings_sets_ie.csv` from `data/real/`, but that file is written to `output/results/` by step 1 — verify/fix the path before running.

### Shared Cognate Analysis (`src/shared_cognates/`)

To analyze shared cognates between phylogenetic subgroups, execute:

```bash
Rscript src/shared_cognates/shared_cognate_theorique.R
Rscript src/shared_cognates/simulate_phylogenies.R
```

### Figures

Once the relevant `output/results/` files exist, generate all figures with:

```bash
Rscript src/plots.R
```

## Data Availability

- **Real Data**: Phylogenetic analyses of Indo-European and Sino-Tibetan language families
- **Simulated Data**: Birth-death trees with varying trait counts (1500, 3000, 6000, 12000)
- **Results**: Processed datasets and statistical summaries

## Citation


Sagart, L., Jacques, G., Lai, Y., Ryder, R. J., Thouzeau, V., Greenhill, S. J., & List, J. M. (2019). Dated language phylogenies shed light on the ancestry of Sino-Tibetan. *Proceedings of the National Academy of Sciences*, 116(21), 10317-10322.

Heggarty, P., Anderson, C., Scarborough, M., King, B., Bouckaert, R., Jocz, L., Kümmel, M. J., Jügel, T., Irslinger, B., Pooth, R., & others. (2023). Language trees with sampled ancestors support a hybrid model for the origin of Indo-European languages. *Science*, 381(6656), eabg0818.
