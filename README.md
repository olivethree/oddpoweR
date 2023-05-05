# oddpoweR: Power analysis for generalized linear mixed-effects models

## Overview

You can use this shiny app to compute simple generalized linear mixed-effects models with binomial or count data as the response variable. Under the hood this app relies on the simr R package (Green & MacLeod, 2016)

## Features

-   You can specify the effect size for each predictor term as a Odds Ratio.

-   A reference table helps you set your Odds ratios based on their equivalence to the traditional Cohen's d effect sizes (small, medium, large). The effect size conversions follow the recommendation by Chen, Ciohen, and Chen (2010).

-   You can specify the overall model intercept and the variability for the random effect terms (subject and trials).

-   You get a CSV file with a time stamped table of the results, so you can share as open science materials (for example, as supporting materials for a preregistration).

-   You get a plot of the power curve along with the recommended sample size (subjects) based on a reliable criterion for the estimated power (i.e. lower bound of the 95% CI for power estimated at a given sample size).

## Work in progress

-   Currently the app only calculates power for two-way or three-way interactions, both in a fully crossed design (i.e., fully within-subjects). The option to calculate for a given main effect (predictor A or B) may be implemented in a future version.

-   The random effects are also only random intercepts only by subject (or any unique identifier for a case) and by trial (or any given class). This app does not aim to cover all possible cases of GLMEMs, and the current plan is to keep this limitation in for the sake of simplicity.

-   Many simulation cases can take many hours (if not days) to run. In such cases it is not recommended to run the simulation using a web app such as this one. I'll work on a solution for such cases.

## **Shiny app**

<https://olivethree.shinyapps.io/oddpoweR/>

## **Citation**

-   Oliveira, M. (2023). oddpoweR: Power analysis for generalized linear mixed-effects models. R Shiny application (Version 0.1), <https://olivethree.shinyapps.io/oddpoweR>

## References

-   Chen, H., Cohen, P., & Chen, S. (2010). How big is a big odds ratio? Interpreting the magnitudes of odds ratios in epidemiological studies. *Communications in Statistics --- Simulation and Computation*, *39*(4), 860-864.

-   Green, P., & MacLeod, C. J. (2016). SIMR: An R package for power analysis of generalized linear mixed models by simulation. *Methods in Ecology and Evolution*, *7*(4), 493-498.

## About me

[Manuel Oliveira](https://manueloliveira.nl/) (a.k.a. Manuel Barbosa de Oliveira)

Currently doing research and teaching at Utrecht University

Statistics \| Data Science \| Social Cognition \| Face Perception \| Psychology of AI

[Google Scholar Page](https://scholar.google.pt/citations?user=8BGdMv8AAAAJ&hl=en)
