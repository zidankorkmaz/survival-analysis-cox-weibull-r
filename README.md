# Survival Analysis in R: Cox PH, Model Selection, and Weibull AFT Models

## Overview
This project applies survival analysis methods to two clinical datasets: relapse times in HIV-positive patients under two treatment regimens, and time-to-blindness in a diabetic retinopathy trial. The analysis covers Kaplan-Meier estimation, Cox proportional hazards modelling with model selection, comparison of tied-event handling methods, and accelerated failure time (AFT) modelling.

## Data

**Dataset 1: HIV relapse times** (n=34, 17 per group)
Time (in days) until CD4 count relapse in HIV-positive patients under a two-drug regimen (AZT + ddC) versus a triple-drug regimen (AZT + ddC + saquinavir), with right-censored observations.

**Dataset 2: Diabetic Retinopathy Study** (`DRS.csv`, n=197)
A clinical trial of laser treatment versus placebo for diabetic retinopathy, where the event of interest is blindness (visual acuity falling below 5/100 at two consecutive tests). Covariates include treatment, age, and risk group (6-12).

## Methods

### Kaplan-Meier estimation and confidence intervals
Kaplan-Meier survival curves were estimated for each HIV treatment group, with 95% confidence intervals computed for the probability of relapse-free survival at 120 days.

### Cox proportional hazards modelling
A Cox PH model was fitted to the HIV data with treatment as a covariate. For the DRS data, a full model (treatment, age, risk group) was reduced via backward AIC selection to identify the preferred model. The proportional hazards assumption was checked using a complementary log-log plot.

### Tied event handling
Since the DRS dataset contains tied survival times, the preferred Cox model was fitted using both the exact partial likelihood and the Efron approximation, and the resulting estimates were compared.

### Weibull AFT model
An accelerated failure time model with a Weibull baseline distribution was fitted to the DRS data using the same covariates as the preferred Cox model, with model fit assessed via a log-log survival plot.

## Results

### HIV relapse analysis
- Kaplan-Meier curves showed a higher observed relapse-free probability in the triple-drug group throughout follow-up.
- 120-day relapse-free survival: 6.9% (two-drug) vs. 43.1% (triple-drug); confidence intervals overlapped.
- Cox PH model: hazard ratio for triple-drug vs. two-drug = 0.564 (95% CI: 0.255, 1.248; p = 0.157) — a 43.6% lower hazard that was not statistically significant.
- The complementary log-log plot showed approximately parallel curves, supporting the proportional hazards assumption.

### Diabetic retinopathy analysis
Backward AIC selection removed age (p = 0.756) from the full model, retaining treatment and risk group:

| Variable | Hazard Ratio | 95% CI | p-value |
|---|---|---|---|
| Treatment (placebo vs laser) | 1.688 | (1.102, 2.585) | 0.016 |
| Risk group | 1.125 | (0.974, 1.300) | 0.109 |

Placebo patients had a 69% higher hazard of blindness compared to laser treatment.

**Tied events:** 22 distinct tied survival times among 77 distinct events. Exact and Efron methods produced nearly identical estimates (e.g. treatment coefficient: 0.534 vs 0.524), indicating no meaningful bias from the tie-handling method chosen.

**Weibull AFT model:**

| Covariate | Estimate | 95% CI |
|---|---|---|
| Treatment (placebo) | -0.645 | (-1.165, -0.124) |
| Risk group | -0.144 | (-0.317, 0.030) |

The negative treatment coefficient indicates placebo is associated with shorter survival times (factor of ≈0.525) compared to laser, consistent with the Cox model findings. The log-log survival plot showed an approximately linear relationship, supporting the Weibull distribution as a reasonable choice for the baseline hazard.

## Interpretation
Across both datasets, treatment effects were directionally consistent with clinical expectation: more intensive HIV therapy and laser treatment for retinopathy were both associated with better outcomes. The DRS analysis provided statistically significant evidence for a treatment effect, while the smaller HIV dataset (n=34) showed a similar direction of effect but did not reach statistical significance, illustrating how sample size affects the ability to detect real treatment differences.

Comparing exact and Efron methods for tied events confirmed that, for this dataset, the choice of tie-handling approach does not materially affect conclusions — a useful validation step in confirming model robustness.

## Files
- `MATH3085_R code.R` — full R code for all analyses
- `MATH3085_Report.pdf` — full written report with figures
- `DRS.csv` — diabetic retinopathy dataset
