# Top-k Quantile Regression Trees for Probabilistic Forecasting

Random Forest-based probabilistic forecasting applied to bike-sharing data from EcoBici Guadalajara. This project evaluates the Top-k simplification of Quantile Regression Forests for improving predictive performance and interpretability.

> **Status:** Active research project. Based on the paper *“Simplifying Random Forests' Probabilistic Forecasts”*.

---

## Overview

This repository contains the implementation developed as part of an applied machine learning study on probabilistic forecasting.

The project focuses on predicting bike trip duration using Random Forest models interpreted as conditional distributions, and introduces a Top-k sparsification strategy over the learned weight matrix.

The main objective is to evaluate whether reducing the full distribution to a small number of representative training samples improves prediction quality and stability.

---

## Repository Structure

```
.
├── notebooks/
│   └── topk_ecobici.ipynb
│
├── data/
│   ├── raw/
│   └── processed/
│
├── results/
│   ├── rf_full/
│   ├── topk_k3/
│   ├── topk_k5/
│   └── topk_k10/
│
├── figures/
│   ├── predictions_vs_true.png
│   └── residuals.png
│
├── docs/
│   └── report.pdf
│
├── requirements.txt
├── LICENSE
└── README.md
```

---

## Features

- Random Forest implemented from scratch (1000 trees)
- Bootstrap aggregation with feature subsampling
- Quantile Regression Forest-style weight computation
- Co-occurrence-based leaf similarity matrix
- Top-k sparsification of predictive distributions
- Probabilistic forecasting via weighted empirical distribution
- Evaluation using MAE, RMSE, and R² metrics
- Empirical comparison of k ∈ {3, 5, 10}

---

## Methodology

The model is based on interpreting Random Forest predictions as:

- A weighted combination of training targets
- A conditional empirical distribution over observed outcomes

A weight matrix is constructed by counting co-occurrences of test and training samples in terminal leaves across all trees.

The Top-k method retains only the k largest weights per test sample and renormalizes the resulting distribution.

---

## Technologies

- Python
- NumPy
- Pandas
- Scikit-learn (Decision Trees)
- Matplotlib
- Jupyter Notebook

---

## Results

| Model   | MAE (min) | RMSE (min) | R²    |
|---------|-----------|------------|-------|
| RF full | 5.073     | 6.212      | 0.110 |
| Top-3   | 2.756     | 4.498      | 0.533 |
| Top-5   | 2.774     | 4.390      | 0.555 |
| Top-10  | 2.994     | 4.505      | 0.532 |

The best performance is obtained at k = 5, indicating that the predictive mass concentrates in a small subset of similar training samples.

---

## Research Context

This project is part of an applied study in machine learning for probabilistic forecasting.

It explores:

- Random Forests as kernel-like estimators
- Empirical conditional distributions
- Sparsification of prediction weights
- Trade-off between bias and variance in probabilistic predictions

---

## Future Work

- Formal comparison against classical Quantile Regression Forests
- Extension to temporal models (sequence-aware forecasting)
- Scalability improvements for full dataset weight computation
- Integration with uncertainty calibration metrics
- Application to other mobility datasets

---

## Author

**Kaleb Aguilar**

Physics Engineer  
M.Sc. Student in Computer Science

LinkedIn:  
https://linkedin.com/in/kaleb-aguilar-70971129a