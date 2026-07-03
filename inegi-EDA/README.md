# Exploratory Data Analysis: Birth Date Patterns in Mexico

Analysis of birth date distributions across Mexican states using INEGI survey data. The project uncovers seasonal and regional patterns in birth frequency, explores the 9-month conception lag hypothesis, and includes a Monte Carlo simulation of the arc-sine law for random walks.

> **Language:** R (ggplot2 + tidyverse)

---

## Overview

This project consists of two self-contained analyses:

- **Birth date EDA** — Day-of-year and monthly birth frequency analysis from INEGI microdata, broken down by state. Births are normalized by days per month to remove calendar bias, and a 9-month offset trend line is overlaid to visually test whether birth peaks align with conception windows (holidays, festivals).

- **Arc-sine law simulation** — Monte Carlo estimation of the fraction of time a symmetric random walk stays positive, illustrating the arc-sine distribution from probability theory.

---

## Data Sources

- **`inegi.csv.CSV`** — INEGI survey microdata containing individual-level records with birth day (`fn_dia`), birth month (`fn_mes`), and state code (`estado`).
- **`catalogo_entidades.csv`** — Lookup table mapping INEGI state codes to state names.

Both files are required in the working directory to run the notebook.

---

## Analysis Pipeline

### Part 1 — Birth Date Patterns

1. **Data loading and cleaning**
   - Read INEGI microdata with `readr`.
   - Filter out sentinel values (`fn_dia == 99`, `fn_mes == 99`) and drop remaining `NA` rows.
   - Compute day-of-year from month and day fields.
   - Merge with the state catalog to obtain human-readable state names.

2. **National birth frequency by day of year**
   - 365-bin histogram colored by count using the `viridis` (turbo) palette.
   - Month delimiters overlaid as dashed vertical lines.

3. **State-level breakdowns**
   - Per-state daily histograms for selected states: México, Guanajuato, Nuevo León, Oaxaca, Quintana Roo.
   - Horizontal bar chart of total surveyed individuals per state.

4. **Monthly normalized analysis**
   - Birth counts divided by the number of days in each month to remove calendar length bias.
   - Months reordered starting in April (academic/fiscal alignment).
   - Two trend lines overlaid: raw monthly trend and the same trend shifted by −9 months, testing whether birth peaks trace back to specific conception periods.
   - Repeated for individual states (San Luis Potosí, Nuevo León, Quintana Roo) to check regional variation.

### Part 2 — Arc-Sine Law (Random Walk)

- 10,000 experiments of a 100-step symmetric random walk.
- For each walk, the fraction of time the cumulative sum stays positive is recorded.
- The resulting distribution is plotted as a histogram, illustrating the U-shaped arc-sine law: the walker spends most of its time either mostly winning or mostly losing, rarely near 50%.

---

## Repository Structure

```
.
├── birth_date_eda.ipynb      # Main R notebook
├── inegi.csv.CSV             # INEGI microdata (not included)
├── catalogo_entidades.csv    # State code lookup (not included)
└── README.md
```

---

## Key Visualizations

- Daily birth frequency heatmap-histogram (national and per-state)
- Surveyed population by state (ranked horizontal bar chart)
- Monthly birth rate with 9-month conception offset overlay
- Arc-sine distribution from Monte Carlo random walks

---

## Dependencies

- R ≥ 4.0
- ggplot2, dplyr, tidyr, readr
- viridis

---

## Notes

- The 9-month offset analysis is exploratory and visual — no formal hypothesis test is performed, but the alignment between conception windows and cultural events (holidays, long weekends) is apparent in several states.
- The random walk simulation uses `set.seed(6170)` for reproducibility.
- Data files are not included in the repository due to INEGI redistribution terms; instructions for downloading them can be found at [inegi.org.mx](https://www.inegi.org.mx/).

---

## Author

**Kaleb Aguilar**

Physics Engineer
M.Sc. Student in Computer Science at CIMAT

LinkedIn:
https://linkedin.com/in/kaleb-aguilar-70971129a
