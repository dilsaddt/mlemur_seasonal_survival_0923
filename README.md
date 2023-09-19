# mlemur_seasonal_survival_0923
Seasonal survival analysis of the gray mouse lemur population in Kirindy Forest, Madagascar.
Code for the multistate, Jolly-Seber capture-recapture analysis.
States are juveniles in dry season (JD), juveniles in wet season (JW), adults in dry season (AD), and adults in wet seasin (AW).
Males and females are analysed seperately to include sex effect on survival.
Environmental factors used in the model to check survival patterns under these factors are rainfall, maximum temperature, and population density.


## Files

data/mlemur_example_data_092023 = data for running models and predictions

*** Data here is representative just to make code running, not the actual data used in the analyses. ***

model_run/mlemur_additive_github_0923.R = additive model code
model_run/mlemur_interaction1_github_0923.R = interaction model 1 code
model_run/mlemur_interaction2_github_0923.R = interaction model 2 code

model_diagnostics/mlemur_diagnostics_github_0923.R = model diagnostics for convergence and parameter distributions
here only diagnostics of additive model is shown as an example, but the process is the same for interaction models

prediction/mlemur_prediction_additive_github_0923.R = predictions and visualization from the additive model
prediction/mlemur_prediction_interaction1_github_0923.R = predictions and visualization from the interaction model 1
prediction/mlemur_prediction_interaction1_github_0923.R = predictions and visualization from the interaction model 2

## Software
R version  4.2.1
NIMBLE version 0.13.2 (through R-package 'nimble')
