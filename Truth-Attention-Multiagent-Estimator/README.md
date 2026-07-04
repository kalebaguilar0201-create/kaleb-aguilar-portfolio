# FRC Performance Estimation using Neural Networks and Optimization

Neural network–based model for predicting alliance performance in FRC matches, combined with optimization-based decomposition to estimate individual robot contributions.

> **Status:** Active research project using TBA 2026 season data.

---

## Overview

This project models FRC alliance performance using deep learning and optimization techniques.

It combines:

- A neural network that predicts alliance score distributions
- A structured representation of robots as feature vectors
- An optimization-based estimator to infer individual robot contributions

The goal is to bridge team-level prediction with robot-level performance attribution.

---

## Data Source

Data is collected from The Blue Alliance (TBA) API for the 2026 season.

- ~30,000 qualification matches
- Multiple regional events + championship
- Preprocessed into a unified dataset (`df_concatenado.csv`)

Each record contains:

- Match metadata (event, match number, time)
- Alliance composition (3 robots)
- Scoring breakdown (auto, teleop, endgame)
- Final alliance score

---

## Pipeline Overview

### 1. Data Collection

- Downloaded via TBA API
- Event-level and match-level queries
- Cached per event to avoid redundant calls
- Concatenated into a single dataset

---

### 2. Feature Engineering

Each robot is represented by a high-dimensional vector including:

- One-hot encoding of robot identity
- Historical win-rate statistics per event
- Global performance descriptors (normalized)

These vectors are constructed using only pre-championship data to avoid leakage.

---

### 3. Alliance Modeling

Each match consists of 3 robots per alliance.

A neural network (`AllianceScoreNet`) processes:

- Individual robot embeddings
- Attention mechanism over alliance composition
- Output: discretized score distribution (64 bins)

Training:

- Soft-label targets (Gaussian-smoothed bins)
- Cross-entropy loss + entropy regularization
- SGD optimizer with learning rate scheduling

---

### 4. Individual Performance Estimation (EPA Module)

The model predicts alliance-level performance. To infer individual robot contribution, an inverse optimization problem is solved:

- Linear decomposition of alliance score
- Regularized least squares formulation
- Bayesian-style prior over robot contributions
- Non-negativity constraints

Two solvers are used:

- L-BFGS-B (bounded quasi-Newton optimization)
- Projected Gradient Descent (PGD)

Both are evaluated under different priors:

- Neural network–based prior
- Baseline average-performance prior

---

## Evaluation

Evaluation is performed on 997 championship matches.

Metrics include:

- Match outcome accuracy
- RMSE of predicted alliance scores
- Error in margin prediction
- Calibration of uncertainty

### Results

| Method            | Accuracy | RMSE   |
|------------------|----------|--------|
| PGD + Attention  | 0.802    | 93.32  |
| L-BFGS + Attn    | 0.800    | 93.81  |
| Baseline (avg/3) | 0.768    | 178.99 |

---

## Live Dashboard

Explore per-robot performance estimates, match predictions, and EPA breakdowns interactively:

**[Open Dashboard →](https://datastudio.google.com/reporting/d705f0af-2b10-4587-8c7f-2174319729b9)**

---

## Repository Structure

```
.
├── DataDownloader.ipynb
├── df_concatenado.csv
├── ExpLBFGS.ipynb
├── ExpPGD.ipynb
├── requirements.txt
└── README.md
```

---

## Dependencies

- torch
- numpy
- pandas
- scipy
- scikit-learn
- matplotlib
- requests

---

## Execution Pipeline

1. Run `DataDownloader.ipynb` to generate dataset  
2. Run either experiment notebook (PGD or L-BFGS)  
3. Load trained models and evaluate EPA module  
4. Compute final match-level performance metrics  

---

## Key Contributions

- Neural modeling of FRC alliance performance
- Attention-based aggregation of robot embeddings
- Soft-label training for score distributions
- Optimization-based inverse performance estimation
- Comparison between gradient-based and quasi-Newton solvers

---

## Notes

- Data leakage is prevented via temporal split (pre/post cutoff)
- Models are trained separately for autonomous and teleoperated phases
- EPA module introduces interpretability at robot level

---

## Author

**Kaleb Aguilar**

Physics Engineer  
M.Sc. Student in Computer Science at CIMAT

LinkedIn:  
https://linkedin.com/in/kaleb-aguilar-70971129a
