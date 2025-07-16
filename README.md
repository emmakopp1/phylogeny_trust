# Phylogenetic Reconstruction Accuracy in Linguistic Evolution

## Abstract

This repository contains code and data for evaluating phylogenetic reconstruction accuracy in linguistic evolution. The study examines the performance of Maximum Clade Credibility (MCC) and consensus tree methods across varying data quantities and linguistic families, combining simulation studies with real-world phylogenetic analyses.

## Project Structure

```
phylogeny_trust/
├── data/                           # Raw and processed data files
│   ├── real/                       # Real linguistic phylogenetic data
│   │   ├── iecor_ctmc-strict-M1/   # Indo-European (IE) BEAST results
│   │   └── st_ctmc-strict-fbd-uni/ # Sino-Tibetan (ST) BEAST results
│   ├── simulated-2025-05-13/       # Simulation dataset (850 traits)
│   ├── simulated-2025-07-02-6000/  # Simulation dataset (6000 traits)
│   └── simulated-2025-07-08-12000/ # Simulation dataset (12000 traits)
├── src/                            # Source code organized by analysis type
│   ├── simulations/                # Phylogenetic simulation analyses
│   ├── ancestral_reconstruction/   # Ancestral state reconstruction
│   ├── 01_Import.R                 # Data import and processing
│   ├── 02_Tydi.R                   # Data tidying and summarization
├── output/                         # Generated results and figures
│   ├── results/                    # Processed datasets and statistics
│   └── figs/                       # Publication-ready figures
└── archive/                        # Legacy and experimental code
```

## Methodology

### Simulation Study

The simulation study evaluates phylogenetic reconstruction accuracy using birth-death trees with varying evolutionary ages and trait counts:

- **Tree Simulation**: 50 taxa birth-death trees scaled to 8 ages (1-8 time units)
- **Trait Evolution**: Binary trait evolution under CTMC models
- **Data Sizes**: 3000, 6000, and 12000 traits
- **Reconstruction Methods**: MCC and 50% majority-rule consensus trees
- **Evaluation Metrics**: Node accuracy, first split identification, marginal probabilities

### Real Data Analysis

Phylogenetic reconstruction of multiple language families using:
- **Families**: Indo-European (data outputs files taken from CITE ARTICLE), Sino-Tibetan (input data taken from CITE ARTICLE and but we did the inference)
- **Methods**: BEAST with strict clock and birth-death/fossilized birth-death models
- **Ancestral Reconstruction**: Semantic meaning evolution analysis

## Script Documentation

### Core Analysis Pipeline

#### 1. Data Import and Processing (`src/`)

**`01_Import.R`**
- **Purpose**: Imports and processes real phylogenetic data from multiple language families
- **Input**: BEAST tree files, log files, NEXUS alignments from `data/real/`
- **Output**: 
  - `output/results/*/`*`_tipages.csv.bz` - Tip age distributions
  - `output/results/*/`*`_tracelog.csv` - MCMC trace summaries
  - `output/results/ntipschars.csv` - Taxa and character counts

**`02_Tydi.R`**
- **Purpose**: Tidies imported data and computes effective sample sizes
- **Input**: All CSV files from `output/results/*/`
- **Output**:
  - `output/results/tipages_summary.csv` - Summarized tip ages
  - `output/results/tracelog_summary.csv` - MCMC parameter summaries
  - `output/results/ess.csv` - Effective sample size calculations

#### 2. Simulation Analysis (`src/simulations/`)

##### Main Scripts (`src/simulations/`)

**`01_tree_simulation_main.R`, `01_tree_simulation_6000.R`, `01_tree_simulation_12000.R`**
- **Purpose**: Simulation scripts generating phylogenetic trees and BEAST analyses
- **Input**: `data/beast-data-sim.xml`, `data/ctmc-strict-bd-template.xml`
- **Output**: `data/simulated-{date}/` containing:
  - `tree-sim-{sim}-{age}.tree` - True phylogenetic trees
  - `beast-*.xml` - BEAST configuration files
  - `ctmc-strict-bd-*.trees` - Posterior tree distributions

**`03_compute_number_of_nodes_resumed.R`**
- **Purpose**: Counts internal nodes in MCC and consensus trees
- **Input**: `data/simulated-*/consensus-*.tree`, `data/simulated-*/mcc-*.tree`
- **Output**: `output/results/number_nodes_mcc_cs_{N_traits}.csv`

**`03_marginal_prob_first_split_resumed.R`**
- **Purpose**: Extracts posterior support for first splits in summary trees
- **Input**: `data/simulated-*/mcc-*.tree`, `data/simulated-*/consensus-*.tree`
- **Output**: `output/results/marginal_prob_first_split_mcc_consensus_{N_traits}.csv`

**`03_reconstruction_difficulty.R`**
- **Purpose**: Analyzes reconstruction difficulty across different conditions
- **Input**: Simulation results from various analyses
- **Output**: Reconstruction difficulty metrics

**`04_simulation_analyses_main.R`**
- **Purpose**: Consolidates simulation results and performs statistical modeling
- **Input**: Multiple CSV files from `output/results/`
- **Output**:
  - `output/results/marginal_probability_first_split_ic.csv` - First split probabilities
  - `output/results/count_true_to_cs.rds` - Node accuracy counts
  - `output/results/regression_*.csv` - Logistic regression results
  - `output/results/prop_*.csv` - Reconstruction accuracy proportions

**`04_simulation_analyses_trait_influence.R`**
- **Purpose**: Analyzes the influence of trait number on reconstruction accuracy
- **Input**: Simulation results across different trait counts
- **Output**: Trait influence analysis results

**`05_visualization.R`**
- **Purpose**: Generates publication-ready figures
- **Input**: All processed results from `output/results/`
- **Output**: Multiple PDF figures in `output/figs/`
  - `marginal_probability_first_split_ic.pdf`
  - `barplot_prop_resume_to_true.pdf`
  - `plausible_node.pdf`

##### Additional Scripts (`src/simulations/`)

**`02_compute_mcc.R`**
- **Purpose**: Computes Maximum Clade Credibility trees from posterior distributions
- **Input**: `*.trees` files from simulation directories
- **Output**: `mcc-{age}.tree` files in same directories

**`02_compute_consensus.R`**
- **Purpose**: Computes 50% majority-rule consensus trees
- **Input**: `*.trees` files from simulation directories
- **Output**: `consensus-{age}.tree` files in same directories

**`03_marginal_prob_first_split_ic.R`**
- **Purpose**: Calculates marginal probabilities for first split identification with confidence intervals
- **Input**: `tree-sim-*.tree` and `*.trees` from simulation directories
- **Output**: `marginal_prob_first_split_ic_{start}_{end}.csv`

**`03_resume_to_true_TF.R`**
- **Purpose**: Determines if nodes in summary trees exist in true trees
- **Input**: True trees, consensus trees, MCC trees
- **Output**: `resume_to_true_TF_{N_traits}.csv`

**`03_true_false_uncertain.R`**
- **Purpose**: Classifies true tree nodes as true, false, or uncertain in consensus trees
- **Input**: True trees, consensus trees
- **Output**: `true_false_uncertain_nodes_{N_traits}.csv`

#### 3. Ancestral State Reconstruction (`src/ancestral_reconstruction/`)

**`00_compute_meaning_set.R`**
- **Purpose**: Extracts semantic meaning boundaries from NEXUS files
- **Input**: `data/real/*/` NEXUS files
- **Output**:
  - `output/results/meanings_sets_ie.csv`
  - `output/results/meanings_sets_st.csv`

**`01_ancestral_state_reconstruction.R`**
- **Purpose**: Performs ancestral state reconstruction using Markov models
- **Input**: BEAST trees, NEXUS alignments, meaning sets
- **Output**:
  - `output/results/ancestral_reconstruction_ie.csv`
  - `output/results/ancestral_reconstruction_st.csv`

**`02_post_process.R`**
- **Purpose**: Post-processes ancestral reconstruction results
- **Input**: Ancestral reconstruction results
- **Output**:
  - `output/results/ancestral_reconstruction_summary_ie.csv`
  - `output/results/ancestral_reconstruction_summary_st.csv`

**`03_visualization.R`**
- **Purpose**: Visualizes ancestral reconstruction results
- **Input**: Ancestral reconstruction summaries
- **Output**: `output/figs/ancestral_reconstruction_by_semantic_meaning.pdf`

## Key Output Files

### Statistical Results (`output/results/`)

- **Simulation Accuracy**: `marginal_probability_first_split_*.csv`
- **Node Classification**: `true_false_uncertain_nodes_*.csv`
- **Tree Comparison**: `resume_to_true_TF_*.csv`
- **Regression Analysis**: `regression_*.csv`
- **Ancestral States**: `ancestral_reconstruction_*.csv`
- **Theoretical Bounds**: `bounds_real_tb*.csv`

### Figures (`output/figs/`)

- **Reconstruction Accuracy**: `marginal_probability_first_split_*.pdf`
- **Node Classification**: `barplot_*.pdf`
- **Ancestral Evolution**: `ancestral_reconstruction_by_semantic_meaning.pdf`
- **Trait Influence**: `number_of_traits_influence.png`

## Software Dependencies

- **R (≥4.0)** with packages:
  - `ape`, `phangorn` - Phylogenetic analysis
  - `TreeSim` - Tree simulation
  - `phytools` - Phylogenetic tools
  - `tidyverse` - Data manipulation
  - `broom` - Statistical modeling
  - `parallel` - Parallel computing
  - `here` - Path management

- **External Software**:
  - **BEAST 2** - Bayesian phylogenetic analysis

## Usage

### Simulation Analysis (`src/simulations/`)

#### Simulation (optionnal)
To reproduce simulate new trees and data, follow these steps:

1. **Generate phylogenetic trees**:

```bash
Rscript src/simulations/01_tree_simulation_6000.R
Rscript src/simulations/01_tree_simulation_12000.R
Rscript src/simulations/01_tree_simulation_main.R
```

2. **Run manually inferences using BEAST**

3. **Compute summary trees**:

```bash
Rscript src/simulations/02_compute_consensus.R
Rscript src/simulations/02_compute_mcc.R
```

#### Analysis
To reproduce the simulation analyses, follow these steps:

1. **Analyze results**:

In the files `03_compute_number_of_nodes_resumed` and `03_marginal_prob_first_split_resumed.R` select the corresponding output file:

   - `data/simulated-2025-05-13` → `output_path <- here("output/results/number_nodes_mcc_cs.csv")`
   - `data/simulated-2025-07-02-6000` → `output_path <- here("output/results/number_nodes_mcc_cs_6000.csv")`
   - `data/simulated-2025-07-08-12000` → `output_path <- here("output/results/number_nodes_mcc_cs_12000.csv")`

In the file `03_compute_number_of_nodes_resumed` you should set the `age_init_sim` variable:
   - Set to `1` for `data/simulated-2025-05-13`
   - Set to `8` for the other two studies

Run the analysis files in this order:

   ```bash
   Rscript src/simulations/03_compute_number_of_nodes_resumed.R
   Rscript src/simulations/03_marginal_prob_first_split_ic.R
   Rscript src/simulations/03_marginal_prob_first_split_resumed.R
   Rscript src/simulations/03_reconstruction_difficulty.R
   Rscript src/simulations/03_resume_to_true_TF.R
   Rscript src/simulations/03_true_false_uncertain.R
   ```

   **Important**: The file `src/simulations/03_marginal_prob_first_split_ic.R` was run on a cluster and is computationally costly. 
   We recommend making tests by setting the variable `phylo_length_test` to a small number (~50). 

2. **Consolidation and modeling**:
   ```bash
   Rscript src/simulations/04_simulation_analyses_main.R
   Rscript src/simulations/04_simulation_analyses_trait_influence.R
   ```

3. **Generate figures**:
   ```bash
   Rscript src/simulations/05_visualization.R
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
   Rscript src/ancestral_reconstruction/01_ancestral_state_reconstruction.R
   ```
   **Note**: This is a computationally costly script, runned on a separate cluster. For testing purposes, set `length_phylo <- 2` in the script. To reproduce the full results, use `length_phylo <- 200`.

4. **Post-processing**:
   ```bash
   Rscript src/ancestral_reconstruction/02_post_process.R
   ```

5. **Visualization**:
   ```bash
   Rscript src/ancestral_reconstruction/03_visualization.R
   ```


## Data Availability

- **Real Data**: Phylogenetic analyses of Indo-European and Sino-Tibetan language families
- **Simulated Data**: Birth-death trees with varying trait counts (3000, 6000, 12000)
- **Results**: Processed datasets and statistical summaries

## Citation


Sagart, L., Jacques, G., Lai, Y., Ryder, R. J., Thouzeau, V., Greenhill, S. J., & List, J. M. (2019). Dated language phylogenies shed light on the ancestry of Sino-Tibetan. *Proceedings of the National Academy of Sciences*, 116(21), 10317-10322.

Heggarty, P., Anderson, C., Scarborough, M., King, B., Bouckaert, R., Jocz, L., Kümmel, M. J., Jügel, T., Irslinger, B., Pooth, R., & others. (2023). Language trees with sampled ancestors support a hybrid model for the origin of Indo-European languages. *Science*, 381(6656), eabg0818.