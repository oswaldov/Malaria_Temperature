# Malaria_Temperature
Data files and R code for manuscript entitled "Temperature impacts the environmental suitability for malaria transmission by Anopheles gambiae and Anopheles stephensi"

Cite the code: [![DOI](https://zenodo.org/badge/442938145.svg)](https://zenodo.org/badge/latestdoi/442938145)

# Authors
1. Oswaldo C Villena
2. Sadie J. Ryan
3. Courtney C. Murdock
4. Leah R. Johnson

# Modeling the temperature dependence of suitability 

First, fit a thermal response to each individual trait (e.g., mosquito development rate, MDR) to independent data using a Bayesian approach where Markov Chain Monte Carlo sampling is used to obtain the posterior distribution for each trait (step 1).
Second, incorporate the posterior distributions of each trait into the suitability metric, S(T), model (step 2).
Next, perform a consistency analysis to determine wheter the model accurately represents the behaviour of the study system (step 3).
All code was executed in R-3.6.3.

# Data files
1. traits.csv for step 1
2. Angamb_Pfalc_samps.Rsave for step 2
3. Prevalence_Pfalciparum.csv, STPosterior_AnGamb_Pfalc.Rsave, STPosterior_AnGamb_Pvivax.Rsave, STPosterior_AnSteph_Pfalc.Rsave, STPosterior_AnSteph_Pvivax.Rsave for step 3.
