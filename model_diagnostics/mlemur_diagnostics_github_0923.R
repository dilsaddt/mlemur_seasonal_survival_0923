############################################################################
# Multi-state model with capture-recapture analysis
# with Jolly-Seber model framework
# Code adapted from Kery and Schaub, BPA Book, Chapter 10
# Part 10.3.2

# Checking for model diagnostics
# Same for all three models
# Here shown for additive models

# Date: September 2023
# Author: Dilsad Dagtekin
###########################################################################

## 1. House keeping ----
##########################

rm(list = ls())

setwd("/path_to_files/")

## 1.1. Load libraries ----
###########################

load.libraries <- function(){
  library(nimble)
  library(parallel)
  library(MCMCvis)
  library(ggplot2)
}

load.libraries()

## 1.2. Load individual data ----
#################################


load("/path_to_files/model_output.RData")

str(run_mlemur_additive_0923$samples)

## 2. WAIC ----
###############################################

additive_WAIC <- run_mlemur_additive_0923$WAIC$WAIC
additive_WAIC

## 3. Convergence and distribution check ----
###############################################

# let's give a burn-in of 20000 (decided on checking all samples and seeing where the convergence began)
ni2 <- 100000-20000

run_mlemur_additive_0923_list_cutoff <- list(run_mlemur_additive_0923$samples[[1]][20000:100000,], run_mlemur_additive_0923$samples[[2]][20000:100000,])

## 3.1. Traceplots ----
###############################################

# Traceplots for each parameter in the model
# For checking convergence visually, also checking Rhat values

MCMCtrace(run_mlemur_additive_0923_list_cutoff, params = c("alpha_phi_JD", "alpha_phi_JW", "alpha_phi_AD", "alpha_phi_AW",
                                                           "alpha_p_JD", "alpha_p_JW", "alpha_p_AD", "alpha_p_AW"), pdf=F, ind = T, Rhat = T, n.eff = T, iter=ni2)

MCMCtrace(run_mlemur_additive_0923_list_cutoff, params =  c("beta1_JD", "beta2_JD", "beta3_JD"), pdf=F, ind = T, Rhat = T, n.eff = T, iter=ni2)

MCMCtrace(run_mlemur_additive_0923_list_cutoff, params =  c("beta1_JW", "beta2_JW","beta3_JW"), pdf=F, ind = T, Rhat = T, n.eff = T, iter=ni2)

MCMCtrace(run_mlemur_additive_0923_list_cutoff, params =  c("beta1_AD", "beta2_AD", "beta3_AD"), pdf=F, ind = T, Rhat = T, n.eff = T, iter=ni2)

MCMCtrace(run_mlemur_additive_0923_list_cutoff, params =  c("beta1_AW", "beta2_AW", "beta3_AW"), pdf=F, ind = T, Rhat = T, n.eff = T, iter=ni2)

MCMCtrace(run_mlemur_additive_0923_list_cutoff, params = c("gamma"), pdf=F, ind = T, Rhat = T, n.eff = T, iter=ni2)

MCMCtrace(run_mlemur_additive_0923_list_cutoff, params = c("eps_p_JD"), pdf=F, ind = T, Rhat = T, n.eff = T, iter=ni2)
MCMCtrace(run_mlemur_additive_0923_list_cutoff, params = c("eps_p_JW"), pdf=F, ind = T, Rhat = T, n.eff = T, iter=ni2)
MCMCtrace(run_mlemur_additive_0923_list_cutoff, params = c("eps_p_AD"), pdf=F, ind = T, Rhat = T, n.eff = T, iter=ni2)
MCMCtrace(run_mlemur_additive_0923_list_cutoff, params = c("eps_p_AW"), pdf=F, ind = T, Rhat = T, n.eff = T, iter=ni2)

# summary

(run_mlemur_additive_0923_sum_all <- MCMCsummary(run_mlemur_additive_0923_list_cutoff, 
                                                 params = c("alpha_phi_JD", "alpha_phi_JW", "alpha_phi_AD", "alpha_phi_AW",
                                                            "alpha_p_JD", "alpha_p_JW", "alpha_p_AD", "alpha_p_AW", "gamma",
                                                            "beta1_JD", "beta2_JD", "beta3_JD",
                                                            "beta1_JW", "beta2_JW", "beta3_JW",
                                                            "beta1_AD", "beta2_AD", "beta3_AD",
                                                            "beta1_AW", "beta2_AW", "beta3_AW",
                                                            "eps_p_JD", "eps_p_JW", "eps_p_AD", "eps_p_AW"), round=3))


## 3.2. Rhats and sample size ----
###############################################

# Checking how what is the percentage of Rhats that are over 1.1
(length(which(run_mlemur_additive_0923_sum_all$Rhat>1.1))/length(run_mlemur_additive_0923_sum_all$mean))*100

# Plotting Rhat distribution and efficient sample numbers
par(mfrow = c(1,1))
hist(run_mlemur_additive_0923_sum_all$Rhat)
plot(run_mlemur_additive_0923_sum_all$n.eff)
abline(h=100, col="red")

## 3.3. Parameter distributions ----
###############################################
# Posterior distribution plots to see which parameters overlap with zero (informative or not informative)

# plot settings

panel_background <- panel_bg(fill = 'white')

color_scheme_set("red")

## 3.3.1. Juveniles in Dry Season ----
######################################

JD_add_pardist <- mcmc_areas(run_js_model_samper_F_14_list_cutoff,
                                  pars = c("alpha_phi_JD","beta1_JD","beta2_JD", "beta3_JD"),
                                  area_method = "equal height",
                                  prob = 0.95) + panel_background +
  scale_y_discrete(labels = rev(c("Intercept", "Rainfall", "Maximum Temperature", "Population Density"))) + 
  theme_classic() +
  theme(axis.text = element_text(size=26),
        axis.title = element_text(size=28))

JD_add_pardist


## 3.3.2. Juveniles in Wet Season ----
######################################

JW_add_pardist <- mcmc_areas(run_js_model_samper_F_14_list_cutoff,
                                 pars = c("alpha_phi_JW","beta1_JW","beta2_JW", "beta3_JW"),
                                 area_method = "equal height",
                                 prob = 0.95) + panel_background +
  scale_y_discrete(labels = rev(c("Intercept", "Rainfall", "Maximum Temperature", "Population Density"))) + 
  
  theme_classic() +
  theme(axis.text = element_text(size=26),
        axis.title = element_text(size=28))

JW_add_pardist

## 3.3.3. Adults in Dry Season ----
######################################

AD_add_pardist <- mcmc_areas(run_js_model_samper_F_14_list_cutoff,
                                  pars = c("alpha_phi_AD","beta1_AD","beta2_AD", "beta3_AD"),
                                  area_method = "equal height",
                                  prob = 0.95) + panel_background +
  scale_y_discrete(labels = rev(c("Intercept", "Rainfall", "Maximum Temperature", "Population Density"))) + 
  
  theme_classic() +
  theme(axis.text = element_text(size=26),
        axis.title = element_text(size=28))

AD_add_pardist

## 3.3.4. Adults in Wet Season ----
######################################

AW_add_pardist <- mcmc_areas(run_js_model_samper_F_14_list_cutoff,
                                 pars = c("alpha_phi_AW","beta1_AW","beta2_AW", "beta3_AW"),
                                 area_method = "equal height",
                                 prob = 0.95) + panel_background +
  scale_y_discrete(labels = rev(c("Intercept", "Rainfall", "Maximum Temperature", "Population Density"))) + 
  
  theme_classic() +
  theme(axis.text = element_text(size=26),
        axis.title = element_text(size=28))

AW_add_pardist

