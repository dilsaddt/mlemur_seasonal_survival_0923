############################################################################
# Multi-state model with capture-recapture analysis
# with Jolly-Seber model framework
# Code adapted from Kery and Schaub, BPA Book, Chapter 10
# Part 10.3.2

# Predictions for additive model

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

## 1.2. Load covariate data ----
#################################


load("path_to_data/mlemur_example_data_092023.RData") # the data file is called:

str() # covariate data


# 1.3. Get the model output ----
##################################

load("/path_to_files/model_output.RData")

str(run_mlemur_additive_0923$samples)

run_mlemur_additive_0923_list_cutoff <- list(run_mlemur_additive_0923$samples[[1]][20001:100000,], run_mlemur_additive_0923$samples[[2]][20001:100000,])

(run_mlemur_additive_0923_sum_all <- MCMCsummary(run_mlemur_additive_0923_list_cutoff, 
                                                 params = c("alpha_phi_JD", "alpha_phi_JW", "alpha_phi_AD", "alpha_phi_AW",
                                                            "alpha_p_JD", "alpha_p_JW", "alpha_p_AD", "alpha_p_AW", "gamma",
                                                            "beta1_JD", "beta2_JD", "beta3_JD",
                                                            "beta1_JW", "beta2_JW", "beta3_JW",
                                                            "beta1_AD", "beta2_AD", "beta3_AD",
                                                            "beta1_AW", "beta2_AW", "beta3_AW",
                                                            "eps_p_JD", "eps_p_JW", "eps_p_AD", "eps_p_AW",
                                                            "N_JD", "N_JW", "N_AD", "N_AW"), round=3))



###########################################################################
# 2. Predictions
###########################################################################

## 2.1. SURVIVAL ----
###############################################

# logit(phi_JD) <- alpha + rain_prev_wet + tmax_dry + popsize
# logit(phi_JW)  <- alpha + rain_wet + tmax_prev_dry + popsize
# logit(phi_AD) <- alpha + rain_prev_wet + tmax_dry + popsize
# logit(phi_AW)  <- alpha + rain_wet + tmax_prev_dry + popsize

# JD
(alpha_phi_JD <- c(run_mlemur_additive_0923$samples$chain1[20001:100000,'alpha_phi_JD'],run_mlemur_additive_0923$samples$chain2[20001:100000,'alpha_phi_JD']))
(beta1_JD <- c(run_mlemur_additive_0923$samples$chain1[20001:100000,'beta1_JD'],run_mlemur_additive_0923$samples$chain2[20001:100000,'beta1_JD']))
(beta2_JD <- c(run_mlemur_additive_0923$samples$chain1[20001:100000,'beta2_JD'],run_mlemur_additive_0923$samples$chain2[20001:100000,'beta2_JD']))
(beta3_JD <- c(run_mlemur_additive_0923$samples$chain1[20001:100000,'beta3_JD'],run_mlemur_additive_0923$samples$chain2[20001:100000,'beta3_JD']))

# JW
(alpha_phi_JW <- c(run_mlemur_additive_0923$samples$chain1[20001:100000,'alpha_phi_JW'],run_mlemur_additive_0923$samples$chain2[20001:100000,'alpha_phi_JW']))
(beta1_JW <- c(run_mlemur_additive_0923$samples$chain1[20001:100000,'beta1_JW'],run_mlemur_additive_0923$samples$chain2[20001:100000,'beta1_JW']))
(beta2_JW <- c(run_mlemur_additive_0923$samples$chain1[20001:100000,'beta2_JW'],run_mlemur_additive_0923$samples$chain2[20001:100000,'beta2_JW']))
(beta3_JW <- c(run_mlemur_additive_0923$samples$chain1[20001:100000,'beta3_JW'],run_mlemur_additive_0923$samples$chain2[20001:100000,'beta3_JW']))

# AD
(alpha_phi_AD <- c(run_mlemur_additive_0923$samples$chain1[20001:100000,'alpha_phi_AD'],run_mlemur_additive_0923$samples$chain2[20001:100000,'alpha_phi_AD']))
(beta1_AD <- c(run_mlemur_additive_0923$samples$chain1[20001:100000,'beta1_AD'],run_mlemur_additive_0923$samples$chain2[20001:100000,'beta1_AD']))
(beta2_AD <- c(run_mlemur_additive_0923$samples$chain1[20001:100000,'beta2_AD'],run_mlemur_additive_0923$samples$chain2[20001:100000,'beta2_AD']))
(beta3_AD <- c(run_mlemur_additive_0923$samples$chain1[20001:100000,'beta3_AD'],run_mlemur_additive_0923$samples$chain2[20001:100000,'beta3_AD']))

# AW
(alpha_phi_AW <- c(run_mlemur_additive_0923$samples$chain1[20001:100000,'alpha_phi_AW'],run_mlemur_additive_0923$samples$chain2[20001:100000,'alpha_phi_AW']))
(beta1_AW <- c(run_mlemur_additive_0923$samples$chain1[20001:100000,'beta1_AW'],run_mlemur_additive_0923$samples$chain2[20001:100000,'beta1_AW']))
(beta2_AW <- c(run_mlemur_additive_0923$samples$chain1[20001:100000,'beta2_AW'],run_mlemur_additive_0923$samples$chain2[20001:100000,'beta2_AW']))
(beta3_AW <- c(run_mlemur_additive_0923$samples$chain1[20001:100000,'beta3_AW'],run_mlemur_additive_0923$samples$chain2[20001:100000,'beta3_AW']))

rain_prev_wet_post_d <- rain_prev_wet_post

rain_wet_d <- rain_wet

tmax_dry_d <- tmax_dry

tmax_prev_dry_pre_d <- tmax_prev_dry_pre

popsize_d <- popsizeModel

## 2.1.1. JD ----
###############################################
# logit(phi_JD) <- alpha + rain_prev_wet + tmax_dry + popsize

###############################################
## RAINFALL OF THE PREVIOUS WET SEASON ----
###############################################

# length.out rain_prev_wet X mcmc samples
surv_JD_rain_prev_wet <- array(NA, dim = c(length(rain_prev_wet_post_d$prec_std), (nrow(run_mlemur_additive_0923$samples$chain1[20001:100000,])*2)))

for (i in 1:length(rain_prev_wet_post_d$prec_std)){
  
  surv_JD_rain_prev_wet[i,] <- plogis(alpha_phi_JD 
                                         + beta1_JD*rain_prev_wet_post_d$prec_std[i]
                                         + beta2_JD*0 # mean tmax value = 0
                                         + beta3_JD*0) # mean popsize value = 0
  
  
}


str(surv_JD_rain_prev_wet)

# then we take the mean of the mcmc list
pm.surv_JD_rain_prev_wet <- apply(surv_JD_rain_prev_wet, 1, mean)
str(pm.surv_JD_rain_prev_wet)

# then calculate the credible intervals
CRI.surv_JD_rain_prev_wet <- apply(surv_JD_rain_prev_wet, 1, function(x) quantile(x, c(0.025, 0.975)))
str(CRI.surv_JD_rain_prev_wet)

surv_prob_JD_rain_prev_wet <- data.frame(state = "JD", 
                                            rain_prev_wet_std = rain_prev_wet_post_d$prec_std,
                                            rain_prev_wet = rain_prev_wet_post_d$prec,
                                            param = "survival",
                                            pred = pm.surv_JD_rain_prev_wet,
                                            lower = CRI.surv_JD_rain_prev_wet[1,],
                                            upper = CRI.surv_JD_rain_prev_wet[2,],
                                            sex = "Females")


surv_prob_JD_rain_prev_wet

surv_prob_JD_rain_prev_wet_add <- surv_prob_JD_rain_prev_wet

JD_rain_prev_wet_plot <- ggplot(surv_prob_JD_rain_prev_wet)+
  geom_ribbon(aes(x= rain_prev_wet, ymin = lower, ymax = upper), fill="blue",alpha=0.3) +
  geom_line(aes(rain_prev_wet, pred), col="blue", lwd=3, linetype=1) +
  theme_classic()+
  ylim(0,1) +
  ylab(expression(paste("Survival Probabilitiy", " ", (phi)))) +
  xlab('Rainfall of the previous wet season (mm)') +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        strip.text = element_text(size=16, face = "bold"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.key.width = unit(1,"cm"),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10))

JD_rain_prev_wet_plot

###############################################
## MAXIMUM TEMPERATURE OF THE DRY SEASON ----
###############################################

# length.out tmax_dry X mcmc samples
surv_JD_tmax_dry <- array(NA, dim = c(length(tmax_dry_d$tmax_std), (nrow(run_mlemur_additive_0923$samples$chain1[20001:100000,])*2)))

for (i in 1:length(tmax_dry_d$tmax_std)){
  
  surv_JD_tmax_dry[i,] <- plogis(alpha_phi_JD 
                                    + beta1_JD*0  # mean prev_rain value = 0
                                    + beta2_JD*tmax_dry_d$tmax_std[i]
                                    + beta3_JD*0) # mean popsize value = 0
}


str(surv_JD_tmax_dry)

# then we take the mean of the mcmc list
pm.surv_JD_tmax_dry <- apply(surv_JD_tmax_dry, 1, mean)
str(pm.surv_JD_tmax_dry)

# then calculate the credible intervals
CRI.surv_JD_tmax_dry <- apply(surv_JD_tmax_dry, 1, function(x) quantile(x, c(0.025, 0.975)))
str(CRI.surv_JD_tmax_dry)

surv_prob_JD_tmax_dry <- data.frame(state = "JD", 
                                       tmax_dry_std = tmax_dry_d$tmax_std,
                                       tmax_dry = tmax_dry_d$tmax,
                                       param = "survival",
                                       pred = pm.surv_JD_tmax_dry,
                                       lower = CRI.surv_JD_tmax_dry[1,],
                                       upper = CRI.surv_JD_tmax_dry[2,],
                                       sex = "Females")


surv_prob_JD_tmax_dry

surv_prob_JD_tmax_dry_add <- surv_prob_JD_tmax_dry

JD_tmax_dry_plot <- ggplot(surv_prob_JD_tmax_dry)+
  geom_ribbon(aes(x= tmax_dry, ymin = lower, ymax = upper), fill = "red",alpha=0.3) +
  geom_line(aes(tmax_dry, pred), col="red", lwd=3, linetype=1) +
  theme_classic()+
  ylim(0,1) +
  ylab(expression(paste("Survival Probabilitiy", " ", (phi)))) +
  xlab('Maximum temperature of the dry season (°C)') +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        strip.text = element_text(size=16, face = "bold"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.key.width = unit(1,"cm"),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10))

JD_tmax_dry_plot

###############################################
## POPULATION SIZE ----
###############################################

# length.out tmax_dry X mcmc samples
surv_JD_popsize <- array(NA, dim = c(length(popsize_d$pop_std), (nrow(run_mlemur_additive_0923$samples$chain1[20001:100000,])*2)))

for (i in 1:length(popsize_d$pop_std)){
  
  surv_JD_popsize[i,] <- plogis(alpha_phi_JD 
                                   + beta1_JD*0  # mean prev_rain value = 0
                                   + beta2_JD*0  # mean tmax value = 0
                                   + beta3_JD*popsize_d$pop_std[i])
}


str(surv_JD_popsize)

# then we take the mean of the mcmc list
pm.surv_JD_popsize <- apply(surv_JD_popsize, 1, mean)
str(pm.surv_JD_popsize)

# then calculate the credible intervals
CRI.surv_JD_popsize <- apply(surv_JD_popsize, 1, function(x) quantile(x, c(0.025, 0.975)))
str(CRI.surv_JD_popsize)

surv_prob_JD_popsize <- data.frame(state = "JD", 
                                      popsize_std = popsize_d$pop_std,
                                      popsize= popsize_d$pop,
                                      param = "survival",
                                      pred = pm.surv_JD_popsize,
                                      lower = CRI.surv_JD_popsize[1,],
                                      upper = CRI.surv_JD_popsize[2,],
                                      sex = "Females")


surv_prob_JD_popsize

surv_prob_JD_popsize_add <- surv_prob_JD_popsize

JD_popsize_plot <- ggplot(surv_prob_JD_popsize_add)+
  geom_ribbon(aes(x= popsize, ymin = lower, ymax = upper), fill = "red",alpha=0.3) +
  geom_line(aes(popsize, pred), col="red", lwd=3, linetype=1) +
  theme_classic()+
  ylim(0,1) +
  ylab(expression(paste("Survival Probabilitiy", " ", (phi)))) +
  xlab('Population Size') +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        strip.text = element_text(size=16, face = "bold"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.key.width = unit(1,"cm"),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10))

JD_popsize_plot


## 2.1.2. JW ----
###############################################
# logit(phi_JW)  <- alpha + rain_wet + tmax_prev_dry + popsize

###############################################
## RAINFALL OF THE WET SEASON ----
###############################################

# length.out rain_prev_wet X mcmc samples
surv_JW_rain_wet <- array(NA, dim = c(length(rain_wet_d$prec_std), (nrow(run_mlemur_additive_0923$samples$chain1[20001:100000,])*2)))

for (i in 1:length(rain_wet_d$prec_std)){
  
  surv_JW_rain_wet[i,] <- plogis(alpha_phi_JW 
                                   + beta1_JW*rain_wet_d$prec_std[i]
                                   + beta2_JW*0 # mean tmax_prev value = 0
                                   + beta3_JW*0) # mean popsize value = 0
}


str(surv_JW_rain_wet)

# then we take the mean of the mcmc list
pm.surv_JW_rain_wet <- apply(surv_JW_rain_wet, 1, mean)
str(pm.surv_JW_rain_wet)

# then calculate the credible intervals
CRI.surv_JW_rain_wet <- apply(surv_JW_rain_wet, 1, function(x) quantile(x, c(0.025, 0.975)))
str(CRI.surv_JW_rain_wet)

surv_prob_JW_rain_wet <- data.frame(state = "JW", 
                                      rain_wet_std = rain_wet_d$prec_std,
                                      rain_wet = rain_wet_d$prec,
                                      param = "survival",
                                      pred = pm.surv_JW_rain_wet,
                                      lower = CRI.surv_JW_rain_wet[1,],
                                      upper = CRI.surv_JW_rain_wet[2,],
                                      sex = "Females")


surv_prob_JW_rain_wet

surv_prob_JW_rain_wet_add <- surv_prob_JW_rain_wet

JW_rain_wet_plot <- ggplot(surv_prob_JW_rain_wet)+
  geom_ribbon(aes(x= rain_wet, ymin = lower, ymax = upper), fill="darkgreen",alpha=0.3) +
  geom_line(aes(rain_wet, pred), col="darkgreen", lwd=3, linetype=1) +
  theme_classic()+
  ylim(0,1) +
  ylab(expression(paste("Survival Probabilitiy", " ", (phi)))) +
  xlab('Rainfall of the wet season (mm)') +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        strip.text = element_text(size=16, face = "bold"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.key.width = unit(1,"cm"),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10))

JW_rain_wet_plot

#######################################################
## MAXIMUM TEMPERATURE OF THE PREVIOUS DRY SEASON ----
#######################################################

# length.out tmax_prev_dry_pre X mcmc samples
surv_JW_tmax_prev_dry <- array(NA, dim = c(length(tmax_prev_dry_pre_d$tmax_std), (nrow(run_mlemur_additive_0923$samples$chain1[20001:100000,])*2)))

for (i in 1:length(tmax_prev_dry_pre_d$tmax_std)){
  
  surv_JW_tmax_prev_dry[i,] <-plogis(alpha_phi_JW 
                                       + beta1_JW*0 # mean rain_wet value = 0
                                       + beta2_JW*tmax_prev_dry_pre_d$tmax_std[i]
                                       + beta3_JW*0) # mean popsize value = 0
}


str(surv_JW_tmax_prev_dry)

# then we take the mean of the mcmc list
pm.surv_JW_tmax_prev_dry <- apply(surv_JW_tmax_prev_dry, 1, mean)
str(pm.surv_JW_tmax_prev_dry)

# then calculate the credible intervals
CRI.surv_JW_tmax_prev_dry <- apply(surv_JW_tmax_prev_dry, 1, function(x) quantile(x, c(0.025, 0.975)))
str(CRI.surv_JW_tmax_prev_dry)

surv_prob_JW_tmax_prev_dry <- data.frame(state = "JW", 
                                           tmax_prev_dry_std = tmax_prev_dry_pre_d$tmax_std,
                                           tmax_prev_dry = tmax_prev_dry_pre_d$tmax,
                                           param = "survival",
                                           pred = pm.surv_JW_tmax_prev_dry,
                                           lower = CRI.surv_JW_tmax_prev_dry[1,],
                                           upper = CRI.surv_JW_tmax_prev_dry[2,],
                                           sex = "Females")


surv_prob_JW_tmax_prev_dry

surv_prob_JW_tmax_prev_dry_add <- surv_prob_JW_tmax_prev_dry

JW_tmax_prev_dry_plot <- ggplot(surv_prob_JW_tmax_prev_dry)+
  geom_ribbon(aes(x= tmax_prev_dry, ymin = lower, ymax = upper), fill="blue",alpha=0.3) +
  geom_line(aes(tmax_prev_dry, pred), col="blue", lwd=3, linetype=1) +
  theme_classic()+   ylim(0,1) +
  ylab(expression(paste("Survival Probabilitiy", " ", (phi)))) +
  xlab('Maximum temperature of the previous dry season (°C)') +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        strip.text = element_text(size=16, face = "bold"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.key.width = unit(1,"cm"),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10))

JW_tmax_prev_dry_plot

#######################################################
## POPULATION SIZE ----
#######################################################

# length.out tmax_prev_dry_pre X mcmc samples
surv_JW_popsize <- array(NA, dim = c(length(popsize_d$pop_std), (nrow(run_mlemur_additive_0923$samples$chain1[20001:100000,])*2)))

for (i in 1:length(popsize_d$pop_std)){
  
  surv_JW_popsize[i,] <-plogis(alpha_phi_JW 
                                 + beta1_JW*0 # mean rain_wet value = 0
                                 + beta2_JW*0 # mean tmax_prev value = 0
                                 + beta3_JW*popsize_d$pop_std[i]) 
}


str(surv_JW_popsize)

# then we take the mean of the mcmc list
pm.surv_JW_popsize <- apply(surv_JW_popsize, 1, mean)
str(pm.surv_JW_popsize)

# then calculate the credible intervals
CRI.surv_JW_popsize <- apply(surv_JW_popsize, 1, function(x) quantile(x, c(0.025, 0.975)))
str(CRI.surv_JW_popsize)

surv_prob_JW_popsize <- data.frame(state = "JW", 
                                     popsize_std = popsize_d$pop_std,
                                     popsize = popsize_d$pop,
                                     param = "survival",
                                     pred = pm.surv_JW_popsize,
                                     lower = CRI.surv_JW_popsize[1,],
                                     upper = CRI.surv_JW_popsize[2,],
                                     sex = "Females")


surv_prob_JW_popsize

surv_prob_JW_popsize_add <- surv_prob_JW_popsize

JW_popsize_plot <- ggplot(surv_prob_JW_popsize)+
  geom_ribbon(aes(x= popsize, ymin = lower, ymax = upper), fill="blue",alpha=0.3) +
  geom_line(aes(popsize, pred), col="blue", lwd=3, linetype=1) +
  theme_classic()+
  ylim(0,1) +
  ylab(expression(paste("Survival Probabilitiy", " ", (phi)))) +
  xlab('Population Size') +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        strip.text = element_text(size=16, face = "bold"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.key.width = unit(1,"cm"),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10))

JW_popsize_plot

## 2.1.3. AD ----
###############################################
# logit(phi_AD) <- alpha + rain_prev_wet + tmax_dry + popsize

###############################################
## RAINFALL OF THE PREVIOUS WET SEASON ----
###############################################

# length.out rain_prev_wet X mcmc samples
surv_AD_rain_prev_wet <- array(NA, dim = c(length(rain_prev_wet_post_d$prec_std), (nrow(run_mlemur_additive_0923$samples$chain1[20001:100000,])*2)))

for (i in 1:length(rain_prev_wet_post_d$prec_std)){
  
  surv_AD_rain_prev_wet[i,] <- plogis(alpha_phi_AD 
                                         + beta1_AD*rain_prev_wet_post_d$prec_std[i]
                                         + beta2_AD*0 # mean tmax value = 0
                                         + beta3_AD*0) # mean popsize value = 0
  
  
}


str(surv_AD_rain_prev_wet)

# then we take the mean of the mcmc list
pm.surv_AD_rain_prev_wet <- apply(surv_AD_rain_prev_wet, 1, mean)
str(pm.surv_AD_rain_prev_wet)

# then calculate the credible intervals
CRI.surv_AD_rain_prev_wet <- apply(surv_AD_rain_prev_wet, 1, function(x) quantile(x, c(0.025, 0.975)))
str(CRI.surv_AD_rain_prev_wet)

surv_prob_AD_rain_prev_wet <- data.frame(state = "AD", 
                                            rain_prev_wet_std = rain_prev_wet_post_d$prec_std,
                                            rain_prev_wet = rain_prev_wet_post_d$prec,
                                            param = "survival",
                                            pred = pm.surv_AD_rain_prev_wet,
                                            lower = CRI.surv_AD_rain_prev_wet[1,],
                                            upper = CRI.surv_AD_rain_prev_wet[2,],
                                            sex = "Females")


surv_prob_AD_rain_prev_wet

surv_prob_AD_rain_prev_wet_add <- surv_prob_AD_rain_prev_wet

AD_rain_prev_wet_plot <- ggplot(surv_prob_AD_rain_prev_wet)+
  geom_ribbon(aes(x= rain_prev_wet, ymin = lower, ymax = upper), fill="blue",alpha=0.3) +
  geom_line(aes(rain_prev_wet, pred), col="blue", lwd=3, linetype=1) +
  theme_classic()+
  ylim(0,1) +
  ylab(expression(paste("Survival Probabilitiy", " ", (phi)))) +
  xlab('Rainfall of the previous wet season (mm)') +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        strip.text = element_text(size=16, face = "bold"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.key.width = unit(1,"cm"),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10))

AD_rain_prev_wet_plot

###############################################
## MAXIMUM TEMPERATURE OF THE DRY SEASON ----
###############################################

# length.out tmax_dry X mcmc samples
surv_AD_tmax_dry <- array(NA, dim = c(length(tmax_dry_d$tmax_std), (nrow(run_mlemur_additive_0923$samples$chain1[20001:100000,])*2)))

for (i in 1:length(tmax_dry_d$tmax_std)){
  
  surv_AD_tmax_dry[i,] <- plogis(alpha_phi_AD 
                                    + beta1_AD*0  # mean prev_rain value = 0
                                    + beta2_AD*tmax_dry_d$tmax_std[i]
                                    + beta3_AD*0) # mean popsize value = 0
}


str(surv_AD_tmax_dry)

# then we take the mean of the mcmc list
pm.surv_AD_tmax_dry <- apply(surv_AD_tmax_dry, 1, mean)
str(pm.surv_AD_tmax_dry)

# then calculate the credible intervals
CRI.surv_AD_tmax_dry <- apply(surv_AD_tmax_dry, 1, function(x) quantile(x, c(0.025, 0.975)))
str(CRI.surv_AD_tmax_dry)

surv_prob_AD_tmax_dry <- data.frame(state = "AD", 
                                       tmax_dry_std = tmax_dry_d$tmax_std,
                                       tmax_dry = tmax_dry_d$tmax,
                                       param = "survival",
                                       pred = pm.surv_AD_tmax_dry,
                                       lower = CRI.surv_AD_tmax_dry[1,],
                                       upper = CRI.surv_AD_tmax_dry[2,],
                                       sex = "Females")


surv_prob_AD_tmax_dry

surv_prob_AD_tmax_dry_add <- surv_prob_AD_tmax_dry

AD_tmax_dry_plot <- ggplot(surv_prob_AD_tmax_dry)+
  geom_ribbon(aes(x= tmax_dry, ymin = lower, ymax = upper), fill = "red",alpha=0.3) +
  geom_line(aes(tmax_dry, pred), col="red", lwd=3, linetype=1) +
  theme_classic()+
  ylim(0,1) +
  ylab(expression(paste("Survival Probabilitiy", " ", (phi)))) +
  xlab('Maximum temperature of the dry season (°C)') +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        strip.text = element_text(size=16, face = "bold"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.key.width = unit(1,"cm"),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10))

AD_tmax_dry_plot

###############################################
## POPULATION SIZE ----
###############################################

# length.out tmax_dry X mcmc samples
surv_AD_popsize <- array(NA, dim = c(length(popsize_d$pop_std), (nrow(run_mlemur_additive_0923$samples$chain1[20001:100000,])*2)))

for (i in 1:length(popsize_d$pop_std)){
  
  surv_AD_popsize[i,] <- plogis(alpha_phi_AD 
                                   + beta1_AD*0  # mean prev_rain value = 0
                                   + beta2_AD*0  # mean tmax value = 0
                                   + beta3_AD*popsize_d$pop_std[i])
}


str(surv_AD_popsize)

# then we take the mean of the mcmc list
pm.surv_AD_popsize <- apply(surv_AD_popsize, 1, mean)
str(pm.surv_AD_popsize)

# then calculate the credible intervals
CRI.surv_AD_popsize <- apply(surv_AD_popsize, 1, function(x) quantile(x, c(0.025, 0.975)))
str(CRI.surv_AD_popsize)

surv_prob_AD_popsize <- data.frame(state = "AD", 
                                      popsize_std = popsize_d$pop_std,
                                      popsize= popsize_d$pop,
                                      param = "survival",
                                      pred = pm.surv_AD_popsize,
                                      lower = CRI.surv_AD_popsize[1,],
                                      upper = CRI.surv_AD_popsize[2,],
                                      sex = "Females")


surv_prob_AD_popsize

surv_prob_AD_popsize_add <- surv_prob_AD_popsize

AD_popsize_plot <- ggplot(surv_prob_AD_popsize_add)+
  geom_ribbon(aes(x= popsize, ymin = lower, ymax = upper), fill = "red",alpha=0.3) +
  geom_line(aes(popsize, pred), col="red", lwd=3, linetype=1) +
  theme_classic()+
  ylim(0,1) +
  ylab(expression(paste("Survival Probabilitiy", " ", (phi)))) +
  xlab('Population Size') +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        strip.text = element_text(size=16, face = "bold"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.key.width = unit(1,"cm"),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10))

AD_popsize_plot

## 2.1.4. AW ----
###############################################
# logit(phi_AW)  <- alpha + rain_wet + tmax_prev_dry + popsize

###############################################
## RAINFALL OF THE WET SEASON ----
###############################################

# length.out rain_prev_wet X mcmc samples
surv_AW_rain_wet <- array(NA, dim = c(length(rain_wet_d$prec_std), (nrow(run_mlemur_additive_0923$samples$chain1[20001:100000,])*2)))

for (i in 1:length(rain_wet_d$prec_std)){
  
  surv_AW_rain_wet[i,] <- plogis(alpha_phi_AW 
                                   + beta1_AW*rain_wet_d$prec_std[i]
                                   + beta2_AW*0 # mean tmax_prev value = 0
                                   + beta3_AW*0) # mean popsize value = 0
}


str(surv_AW_rain_wet)

# then we take the mean of the mcmc list
pm.surv_AW_rain_wet <- apply(surv_AW_rain_wet, 1, mean)
str(pm.surv_AW_rain_wet)

# then calculate the credible intervals
CRI.surv_AW_rain_wet <- apply(surv_AW_rain_wet, 1, function(x) quantile(x, c(0.025, 0.975)))
str(CRI.surv_AW_rain_wet)

surv_prob_AW_rain_wet <- data.frame(state = "AW", 
                                      rain_wet_std = rain_wet_d$prec_std,
                                      rain_wet = rain_wet_d$prec,
                                      param = "survival",
                                      pred = pm.surv_AW_rain_wet,
                                      lower = CRI.surv_AW_rain_wet[1,],
                                      upper = CRI.surv_AW_rain_wet[2,],
                                      sex = "Females")


surv_prob_AW_rain_wet

surv_prob_AW_rain_wet_add <- surv_prob_AW_rain_wet

AW_rain_wet_plot <- ggplot(surv_prob_AW_rain_wet)+
  geom_ribbon(aes(x= rain_wet, ymin = lower, ymax = upper), fill="darkgreen",alpha=0.3) +
  geom_line(aes(rain_wet, pred), col="darkgreen", lwd=3, linetype=1) +
  theme_classic()+
  ylim(0,1) +
  ylab(expression(paste("Survival Probabilitiy", " ", (phi)))) +
  xlab('Rainfall of the wet season (mm)') +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        strip.text = element_text(size=16, face = "bold"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.key.width = unit(1,"cm"),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10))

AW_rain_wet_plot

#######################################################
## MAXIMUM TEMPERATURE OF THE PREVIOUS DRY SEASON ----
#######################################################

# length.out tmax_prev_dry_pre X mcmc samples
surv_AW_tmax_prev_dry <- array(NA, dim = c(length(tmax_prev_dry_pre_d$tmax_std), (nrow(run_mlemur_additive_0923$samples$chain1[20001:100000,])*2)))

for (i in 1:length(tmax_prev_dry_pre_d$tmax_std)){
  
  surv_AW_tmax_prev_dry[i,] <-plogis(alpha_phi_AW 
                                       + beta1_AW*0 # mean rain_wet value = 0
                                       + beta2_AW*tmax_prev_dry_pre_d$tmax_std[i]
                                       + beta3_AW*0) # mean popsize value = 0
}


str(surv_AW_tmax_prev_dry)

# then we take the mean of the mcmc list
pm.surv_AW_tmax_prev_dry <- apply(surv_AW_tmax_prev_dry, 1, mean)
str(pm.surv_AW_tmax_prev_dry)

# then calculate the credible intervals
CRI.surv_AW_tmax_prev_dry <- apply(surv_AW_tmax_prev_dry, 1, function(x) quantile(x, c(0.025, 0.975)))
str(CRI.surv_AW_tmax_prev_dry)

surv_prob_AW_tmax_prev_dry <- data.frame(state = "AW", 
                                           tmax_prev_dry_std = tmax_prev_dry_pre_d$tmax_std,
                                           tmax_prev_dry = tmax_prev_dry_pre_d$tmax,
                                           param = "survival",
                                           pred = pm.surv_AW_tmax_prev_dry,
                                           lower = CRI.surv_AW_tmax_prev_dry[1,],
                                           upper = CRI.surv_AW_tmax_prev_dry[2,],
                                           sex = "Females")


surv_prob_AW_tmax_prev_dry

surv_prob_AW_tmax_prev_dry_add <- surv_prob_AW_tmax_prev_dry

AW_tmax_prev_dry_plot <- ggplot(surv_prob_AW_tmax_prev_dry)+
  geom_ribbon(aes(x= tmax_prev_dry, ymin = lower, ymax = upper), fill="blue",alpha=0.3) +
  geom_line(aes(tmax_prev_dry, pred), col="blue", lwd=3, linetype=1) +
  theme_classic()+   ylim(0,1) +
  ylab(expression(paste("Survival Probabilitiy", " ", (phi)))) +
  xlab('Maximum temperature of the previous dry season (°C)') +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        strip.text = element_text(size=16, face = "bold"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.key.width = unit(1,"cm"),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10))

AW_tmax_prev_dry_plot

#######################################################
## POPULATION SIZE ----
#######################################################

# length.out tmax_prev_dry_pre X mcmc samples
surv_AW_popsize <- array(NA, dim = c(length(popsize_d$pop_std), (nrow(run_mlemur_additive_0923$samples$chain1[20001:100000,])*2)))

for (i in 1:length(popsize_d$pop_std)){
  
  surv_AW_popsize[i,] <-plogis(alpha_phi_AW 
                                 + beta1_AW*0 # mean rain_wet value = 0
                                 + beta2_AW*0 # mean tmax_prev value = 0
                                 + beta3_AW*popsize_d$pop_std[i]) 
}


str(surv_AW_popsize)

# then we take the mean of the mcmc list
pm.surv_AW_popsize <- apply(surv_AW_popsize, 1, mean)
str(pm.surv_AW_popsize)

# then calculate the credible intervals
CRI.surv_AW_popsize <- apply(surv_AW_popsize, 1, function(x) quantile(x, c(0.025, 0.975)))
str(CRI.surv_AW_popsize)

surv_prob_AW_popsize <- data.frame(state = "AW", 
                                     popsize_std = popsize_d$pop_std,
                                     popsize = popsize_d$pop,
                                     param = "survival",
                                     pred = pm.surv_AW_popsize,
                                     lower = CRI.surv_AW_popsize[1,],
                                     upper = CRI.surv_AW_popsize[2,],
                                     sex = "Females")


surv_prob_AW_popsize

surv_prob_AW_popsize_add <- surv_prob_AW_popsize

AW_popsize_plot <- ggplot(surv_prob_AW_popsize)+
  geom_ribbon(aes(x= popsize, ymin = lower, ymax = upper), fill="blue",alpha=0.3) +
  geom_line(aes(popsize, pred), col="blue", lwd=3, linetype=1) +
  theme_classic()+
  ylim(0,1) +
  ylab(expression(paste("Survival Probabilitiy", " ", (phi)))) +
  xlab('Population Size') +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        strip.text = element_text(size=16, face = "bold"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.key.width = unit(1,"cm"),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10))

AW_popsize_plot

## 2.2. RECAPTURE ----
###############################################

# logit(p) <- alpha + eps(year)

(alpha_p_JD <- c(run_mlemur_additive_0923$samples$chain1[20001:100000,'alpha_p_JD'],run_mlemur_additive_0923$samples$chain2[20001:100000,'alpha_p_JD']))
(alpha_p_JW <- c(run_mlemur_additive_0923$samples$chain1[20001:100000,'alpha_p_JW'],run_mlemur_additive_0923$samples$chain2[20001:100000,'alpha_p_JW']))
(alpha_p_AD <- c(run_mlemur_additive_0923$samples$chain1[20001:100000,'alpha_p_AD'],run_mlemur_additive_0923$samples$chain2[20001:100000,'alpha_p_AD']))
(alpha_p_AW <- c(run_mlemur_additive_0923$samples$chain1[20001:100000,'alpha_p_AW'],run_mlemur_additive_0923$samples$chain2[20001:100000,'alpha_p_AW']))

eps_p_JD <- eps_p_JW <- eps_p_AD <- eps_p_AW <- list(NA)

for (i in 1:length(unique(yearModel$yearCat))){
  eps_p_JD[[i]] <- c(run_mlemur_additive_0923$samples$chain1[20001:100000,paste0('eps_p_JD[',i,']')],run_mlemur_additive_0923$samples$chain2[20001:100000,paste0('eps_p_JD[',i,']')])
  eps_p_JW[[i]] <- c(run_mlemur_additive_0923$samples$chain1[20001:100000,paste0('eps_p_JW[',i,']')],run_mlemur_additive_0923$samples$chain2[20001:100000,paste0('eps_p_JW[',i,']')])
  eps_p_AD[[i]] <- c(run_mlemur_additive_0923$samples$chain1[20001:100000,paste0('eps_p_AD[',i,']')],run_mlemur_additive_0923$samples$chain2[20001:100000,paste0('eps_p_AD[',i,']')])
  eps_p_AW[[i]] <- c(run_mlemur_additive_0923$samples$chain1[20001:100000,paste0('eps_p_AW[',i,']')],run_mlemur_additive_0923$samples$chain2[20001:100000,paste0('eps_p_AW[',i,']')])
  
}
eps_p_JD
eps_p_JW
eps_p_AD
eps_p_AW

recap <- array(NA, dim = c(4, (nrow(run_mlemur_additive_0923$samples$chain1[20001:100000,])*2)))

# taking the random year effect in account for calculation, but not plotting it
for(i in 1:length(unique(yearModel$yearCat))){
  
  recap[1,] <- plogis(alpha_p_JD + eps_p_JD[[i]])
  recap[2,] <- plogis(alpha_p_JW + eps_p_JW[[i]])
  recap[3,] <- plogis(alpha_p_AD + eps_p_AD[[i]])
  recap[4,] <- plogis(alpha_p_AW + eps_p_AW[[i]])
  
}

str(recap)

# then we take the mean of the mcmc list
pm.recap <- apply(recap, 1, mean)
str(pm.recap)

# then calculate the credible intervals
CRI.recap <- apply(recap, 1, function(x) quantile(x, c(0.025, 0.975)))
str(CRI.recap)


recap_prob <- data.frame(state = c("JD", "JW", "AD", "AW"),
                         param = "recapture")
recap_prob$pred <- pm.recap
recap_prob$lower <- CRI.recap[1,]
recap_prob$upper <- CRI.recap[2,]
recap_prob$sex <- "Females"

recap_prob

recap_prob_add <- recap_prob


ggplot(recap_prob)+
  geom_errorbar(aes(x= state, ymin = lower, ymax = upper, col=state), lwd = 0.5, width = 0.3) +
  geom_point(aes(state, pred, col=state), size=3) +
  theme_classic()+
  ylim(0,1) +
  ylab(expression(paste("Recapture Probabilitiy", " ", (p)))) +
  xlab('State - Age x Season') +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        strip.text = element_text(size=16, face = "bold"),
        legend.position = "none",
        legend.direction = "horizontal",
        legend.key.width = unit(1,"cm"),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10))

## 2.3. DERIVED ABUNDANCES ----
###############################################

rownames(run_mlemur_additive_0923_sum_all)[189:412]

derived_abundance_mean <- as.data.frame(run_mlemur_additive_0923_sum_all[189:412,c(1,3,5)])
derived_abundance_mean$year <- rep(yearModel$year)
derived_abundance_mean$yearCat <- rep(yearModel$yearCat)
derived_abundance_mean$season <- rep(c("Postwet", "Prewet"))
derived_abundance_mean$season2 <- rep(c("Dry", "Wet"))
derived_abundance_mean$season3 <- rep(c("Dry", "Wet"))
derived_abundance_mean$season3 <- ifelse(derived_abundance_mean$mean == 0, "out_season", derived_abundance_mean$season2)
derived_abundance_mean$season3 <- factor(derived_abundance_mean$season3, levels = c("Dry", "Wet", "out_season"))

derived_abundance_mean$age <- rep(c("Juvenile", "Adult"), each = 112)
derived_abundance_mean$state <- rep(c("Juvenile Dry", "Juvenile Wet", "Adult Dry", "Adult Wet"), each = 56)
derived_abundance_mean$state <- factor(derived_abundance_mean$state, levels = c("Juvenile Dry", "Juvenile Wet", "Adult Dry", "Adult Wet"))

derived_abundance_mean$sex <- "Females"

colnames(derived_abundance_mean) <- c("mean", "lower", "upper", "year", "yearCat", "season", "season2", "season3", "age", "state", "sex")


ggplot(derived_abundance_mean)+
  #geom_errorbar(aes(x= year, ymin = lower, ymax = upper, col=season2), lwd = 0.5, width = 0.3) +
  geom_ribbon(aes(x= year, ymin = lower, ymax = upper, group = state, fill = state), alpha=0.3) +
  #F8766D - red
  #619CFF - blue
  #B5A9A8 - grey
  geom_line(aes(year, mean, col=state), lwd=3) +
  # scale_fill_manual(values = c("#F8766D", "#619CFF", "#B5A9A8")) +
  # scale_color_manual(values = c("#F8766D", "#619CFF", "#B5A9A8")) +
  theme_classic()+
  # stat_smooth(method=lm)+
  facet_wrap(~season) +
  # ylim(0,100) +
  ylab(expression(paste("Recapture Probabilitiy", " ", (p)))) +
  xlab('Derived population size') +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        strip.text = element_text(size=16, face = "bold"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.key.width = unit(1,"cm"),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10))
