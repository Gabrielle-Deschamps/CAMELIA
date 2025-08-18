This repository contains all scripts and data processing steps used in the paper:

CAMELIA – enhancing species assemblage predictions by integrating community indices (submitted to Methods in Ecology and Evolution).

CAMELIA is a neural network framework that refines stacked species distribution model (S-SDM) predictions by integrating community-level information (species richness and functional trait indices). 
This repository includes:
- Data preprocessing
- Model training and prediction
- Evaluation scripts
---

## Repository structure

```
.
├── CAMELIA.ipynb                # Main notebook to run the CAMELIA workflow
├── data/                        # Input data (species probabilities, traits, community indices, etc.)
├── results/                     # Output folder for trained models and evaluation results
├── utilities/
│   ├── FONCTIONS.R              # Functions used by Evaluation.R and Evaluation_degraded_SSDM.R
│   ├── data_utilities.py        # Functions for loading and preprocessing data
│   ├── evaluation.py            # Functions to compute evaluation metrics
│   ├── neural_model.py          # Neural network architecture and training routines
│   ├── viz_utilities.py         # Helper functions for visualisation
├── Evaluation.R                 # R script to evaluate predictions obtained from different methods and generate the figures presented in the article.
├── Evaluation_degraded_SSDM.R   # R script to evaluate and plot results with degraded SDM inputs
└── README.md                    
```

---

## Requirements

- **Python** ≥ 3.8  
  - Required packages: `numpy`, `pandas`, `scikit-learn`, `tensorflow` `keras`, `matplotlib`, `seaborn`
- **R** ≥ 4.0  
  - Required packages: `tidyverse`, `ggplot2`, `cowplot`, `ggpubr`, `data.table`

You can install Python dependencies with:
```bash
pip install numpy pandas scikit-learn tensorflow keras matplotlib seaborn
```
And R dependencies within R:
```r
install.packages(c("tidyverse", "ggplot2", "cowplot", "ggpubr", "data.table"))
```

---

## Workflow to reproduce the results

1. **Prepare data**  
   Place the required input files in the `data/` folder:  
   - Species probability matrices from S-SDMs (here ssdm_predictions.csv) 
   - Predicted community indices in the `community_indices/` folder (richness (SR), trait means (CM), trait standard deviations(CSTD))  
   - Species trait table (here traits.csv)

2. **Run the CAMELIA model**  
   Open and execute `CAMELIA.ipynb`:
   - Loads input data via `utilities/data_utilities.py`
   - Trains the CAMELIA neural network (`utilities/neural_model.py`)
   - Predicts refined species probabilities
   - Saves results in the `results/` folder

3. **Evaluate predictions**  
   - Run `Evaluation.R` to compute metrics (AUC, TSS, Sørensen index, R², RMSE) and generate main figures from the manuscript.
   - Run `Evaluation_degraded_SSDM.R` to reproduce the analysis with degraded SDM inputs.

---