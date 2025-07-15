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
│   │   ├── cluster/                # High-performance cluster scripts
│   │   └── local/                  # Local analysis scripts
│   ├── ancestral_reconstruction/   # Ancestral state reconstruction
│   ├── 01_Import.R                 # Data import and processing
│   ├── 02_Tydi.R                   # Data tidying and summarization
│   └── 03_Transform.R              # Theoretical bounds computation
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
- **Data Sizes**: 850, 6000, and 12000 traits
- **Reconstruction Methods**: MCC and 50% majority-rule consensus trees
- **Evaluation Metrics**: Node accuracy, first split identification, marginal probabilities

### Real Data Analysis

Phylogenetic reconstruction of multiple language families using:
- **Families**: Indo-European, Sino-Tibetan, Bantu, and others
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

**`03_Transform.R`**
- **Purpose**: Computes theoretical bounds for phylogenetic reconstruction
- **Input**: `output/results/tipages_summary.csv`, `output/results/tracelog_summary.csv`
- **Output**:
  - `output/results/bounds_real_tb.csv` - Theoretical reconstruction bounds
  - `output/results/bounds_real_tb_by_sens.csv` - Bounds by semantic category

#### 2. Simulation Analysis (`src/simulations/`)

##### Local Scripts (`src/simulations/local/`)

**`01_tree_simulation_main.R`**
- **Purpose**: Master simulation script generating phylogenetic trees and BEAST analyses
- **Input**: `data/beast-data-sim.xml`, `data/ctmc-strict-bd-template.xml`
- **Output**: `data/simulated-{date}/` containing:
  - `tree-sim-{sim}-{age}.tree` - True phylogenetic trees
  - `beast-*.xml` - BEAST configuration files
  - `ctmc-strict-bd-*.trees` - Posterior tree distributions

**`02_compute_number_of_nodes_resumed.R`**
- **Purpose**: Counts internal nodes in MCC and consensus trees
- **Input**: `data/simulated-*/consensus-*.tree`, `data/simulated-*/mcc-*.tree`
- **Output**: `output/results/number_nodes_mcc_cs_{N_traits}.csv`

**`02_marginal_prob_first_split_resumed.R`**
- **Purpose**: Extracts posterior support for first splits in summary trees
- **Input**: `data/simulated-*/mcc-*.tree`, `data/simulated-*/consensus-*.tree`
- **Output**: `output/results/marginal_prob_first_split_mcc_consensus_{N_traits}.csv`

**`03_simulation_analyses_main.R`**
- **Purpose**: Consolidates simulation results and performs statistical modeling
- **Input**: Multiple CSV files from `output/results/`
- **Output**:
  - `output/results/marginal_probability_first_split_ic.csv` - First split probabilities
  - `output/results/count_true_to_cs.rds` - Node accuracy counts
  - `output/results/regression_*.csv` - Logistic regression results
  - `output/results/prop_*.csv` - Reconstruction accuracy proportions

**`04_visualization.R`**
- **Purpose**: Generates publication-ready figures
- **Input**: All processed results from `output/results/`
- **Output**: Multiple PDF figures in `output/figs/`
  - `marginal_probability_first_split_ic.pdf`
  - `barplot_prop_resume_to_true.pdf`
  - `plausible_node.pdf`

##### Cluster Scripts (`src/simulations/cluster/`)

**`compute_mcc.R`**
- **Purpose**: Computes Maximum Clade Credibility trees from posterior distributions
- **Input**: `*.trees` files from simulation directories
- **Output**: `mcc-{age}.tree` files in same directories

**`compute_consensus.R`**
- **Purpose**: Computes 50% majority-rule consensus trees
- **Input**: `*.trees` files from simulation directories
- **Output**: `consensus-{age}.tree` files in same directories

**`marginal_prob_first_split_ic.R`**
- **Purpose**: Calculates marginal probabilities for first split identification with confidence intervals
- **Input**: `tree-sim-*.tree` and `*.trees` from simulation directories
- **Output**: `marginal_prob_first_split_ic_{start}_{end}.csv`

**`resume_to_true_TF.R`**
- **Purpose**: Determines if nodes in summary trees exist in true trees
- **Input**: True trees, consensus trees, MCC trees
- **Output**: `resume_to_true_TF_{N_traits}.csv`

**`true_false_uncertain.R`**
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
  - **TreeAnnotator** - Tree summarization

## Usage

1. **Data Processing**: Run `src/01_Import.R`, `src/02_Tydi.R`, `src/03_Transform.R`
2. **Simulation**: Execute `src/simulations/local/01_tree_simulation_main.R`
3. **Cluster Analysis**: Submit cluster scripts in `src/simulations/cluster/`
4. **Statistical Analysis**: Run `src/simulations/local/03_simulation_analyses_main.R`
5. **Visualization**: Execute `src/simulations/local/04_visualization.R`

## Data Availability

- **Real Data**: Phylogenetic analyses of Indo-European and Sino-Tibetan language families
- **Simulated Data**: Birth-death trees with varying trait counts (850, 6000, 12000)
- **Results**: Processed datasets and statistical summaries

## Citation

[Citation information to be added upon publication]

## License

[License information to be added]