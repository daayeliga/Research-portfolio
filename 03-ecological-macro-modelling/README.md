# Ecological Macroeconomic Modelling: Carbon Tax & Degrowth Scenarios

**Module:** MSc Development Economics, Advanced Macroeconomics — SOAS University of London
**Author:** David Ayuuya Ayeliga

## Overview

This project discusses simulation results from a **stock-flow-fund ecological macroeconomic model** (DEFINE-SIMPLE, extended with climate damage feedbacks), comparing a baseline trajectory against a higher-carbon-tax scenario and a degrowth scenario across macroeconomic, financial, and environmental variables through to 2070.

## Contents

| File | Description |
|---|---|
| `PartB_Ecological_Macro_Model_Discussion.docx` | Full discussion of model scenarios, results, and the underlying R code used to build and run the model |

## Key Findings

- A **carbon tax alone** modestly reduces emissions and warming but slightly depresses output and profits, as firms absorb higher costs without an offsetting demand stimulus.
- Combined with **higher government spending**, the carbon tax scenario stabilises government debt-to-GDP while preserving its environmental gains — supporting a fiscal-environmental policy mix consistent with post-Keynesian ecological macroeconomics (Fontana & Sawyer, 2016).
- A **degrowth scenario** achieves the strongest environmental outcomes of all three (lowest temperature, emissions, and capital depreciation) but at a severe macro-financial cost: the rate of profit collapses toward zero and firm leverage rises sharply, echoing the financial instability channel identified by Dafermos et al. (2017).
- Government debt-to-GDP falls furthest under degrowth, but this reflects a shrinking GDP denominator against strong carbon tax revenues rather than genuine fiscal health — a reminder that headline debt ratios can mask real economic contraction.

## Methods & Frameworks

A custom-coded **stock-flow-fund ecological macroeconomic model in R**, incorporating climate damage functions (Dietz & Stern, 2015; Weitzman, 2012) and a temperature-emissions link calibrated to Matthews et al. (2021), following the DEFINE modelling tradition (Dafermos, Nikolaidi & Galanis, 2017).
