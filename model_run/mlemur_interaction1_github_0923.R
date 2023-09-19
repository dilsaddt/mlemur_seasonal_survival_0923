############################################################################
# Multi-state model with capture-recapture analysis
# with Jolly-Seber model framework
# Code adapted from Kery and Schaub, BPA Book, Chapter 10
# Part 10.3.2

# The aim of this script is to run the multistate JS models
# with 5 states
# Juvenile-Dry = JD (1) 
# Juvenile-Wet = JW (2)
# Adult-Dry = AD (3)
# Adult-Wet = AW (4)
# not seen = NS (5)

# Interaction model I: rainfall*population size

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

## 1.2. Load data ----
#################################

load("path_to_data/mlemur_example_data_092023.RData") # the data file is called:

str(mlemur_CH) # capture history
str(mlemur_inits) # intial values for the latent state
# covariates are:
str(rain_prev_wet_post)
str(rain_wet)
str(tmax_prev_dry_pre)
str(tmax_dry)
str(popsize)
str(year)

###########################################################################
# 2. Analysis of the model
###########################################################################

## NIMBLE FUNCTION ----
##########################

# nimble function for the gamma values
nimF.gamma <- nimbleFunction(
  run = function(n.occ = double(0) # total number of occasions
                 , postwet.occasion = double(1) # vector indicating which occasions are postwet (1) or prewet (0)
                 , gamma.prior = double(1)
  ){
    
    for(t in 1:n.occ){
      if(postwet.occasion[t]==0){
        gamma.prior[t] <- 0
      }
    }
    
    returnType(double(1)) # we return an numeric vector
    return(gamma.prior)
    
  }
)

assign("nimF.gamma", nimF.gamma, envir = .GlobalEnv) # important for parallelization

# gamma.prior <- runif((nimble_constants$n.occasions-1), 0, 1)
# nimF.gamma(n.occ = (nimble_constants$n.occasions-1), postwet.occasion = postwet_season, gamma.prior = gamma.prior)

## NIMBLE MODEL ----
##########################

mlemur_interactionI_0923 <- nimbleCode( {
  
  #--------------------------------------
  # Parameters:
  # phi: survival probability
  # gamma: removal entry probability
  # p: capture probability
  #--------------------------------------
  # States (S): 6 states
  # 1 not yet entered NE
  # 2 Juvenile Dry JD
  # 3 Juvenile Wet JW
  # 4 Adult Dry AD
  # 5 Adult Wet AW
  # 6 dead
  #
  # Observations (O): 5 states
  # 1 seen as JD
  # 2 seen as JW
  # 3 seen as AD
  # 4 seen as AW
  # 5 not seen
  #
  #-------------------------------------------------
  # Covariates:
  # rainfall of previous wet season on dry season = rain_prev_wet_post
  # maximum temperature of previous dry season on wet season = tmax_prev_dry_pre
  # maximum temperature of dry season = tmax_dry
  # rainfall of wet season = rain_wet
  # population size = popsize
  #-------------------------------------------------
  
  ### PRIORS AND CONSTRAINTS ###
  
  gamma[1:(n.occasions-1)] <- nimF.gamma(n.occ = (n.occasions-1), 
                                         postwet.occasion =  postwet_season[1:(n.occasions-1)],
                                         gamma.prior =  gamma.prior[1:(n.occasions-1)])
  
  # priors for intercepts
  alpha_phi_JD ~ dnorm(0, 1)    # Prior for mean survival of JD
  alpha_phi_JW ~ dnorm(0, 1)    # Prior for mean survival of JW
  alpha_phi_AD ~ dnorm(0, 1)    # Prior for mean survival of AD
  alpha_phi_AW ~ dnorm(0, 1)    # Prior for mean survival of AW
  
  alpha_p_JD ~ dnorm(0, 1)    # Prior for mean capture of JD
  alpha_p_JW ~ dnorm(0, 1)   # Prior for mean capture of JW
  alpha_p_AD ~ dnorm(0, 1)    # Prior for mean capture of AD
  alpha_p_AW ~ dnorm(0, 1)    # Prior for mean capture of AW
  
  # priors for slopes
  
  # JD
  beta1_JD ~ dnorm(0, 1)
  beta2_JD ~ dnorm(0, 1)
  beta3_JD ~ dnorm(0, 1)
  
  # JW
  beta1_JW ~ dnorm(0, 1)
  beta2_JW ~ dnorm(0, 1)
  beta3_JW ~ dnorm(0, 1)
  
  # AD
  beta1_AD ~ dnorm(0, 1)
  beta2_AD ~ dnorm(0, 1)
  beta3_AD ~ dnorm(0, 1)
  
  # AW
  beta1_AW ~ dnorm(0, 1)
  beta2_AW ~ dnorm(0, 1)
  beta3_AW ~ dnorm(0, 1)
  
  # priors for random effect #
  
  sigma_p_JD ~ dunif(0, 10)    # Prior for standard deviation of p
  sigma_p_JW ~ dunif(0, 10)
  sigma_p_AD ~ dunif(0, 10)
  sigma_p_AW ~ dunif(0, 10)
  
  for (t in 1:(n.years)){
    
    eps_p_JD[t] ~ dnorm(0, sd = sigma_p_JD)
    eps_p_JW[t] ~ dnorm(0, sd = sigma_p_JW)
    eps_p_AD[t] ~ dnorm(0, sd = sigma_p_AD)
    eps_p_AW[t] ~ dnorm(0, sd = sigma_p_AW)
    
  }
  
  for (t in 1:(n.occasions-1)){ # occasions are seasons
    
    gamma.prior[t] ~ dunif(0, 1) # Prior for entry probabilities
    
    ### LIKELIHOOD ###
    
    logit(phi_JD[t]) <- alpha_phi_JD + beta1_JD*rain_prev_wet_post[t] + beta2_JD*tmax_dry[t] + beta3_JD*popsize[t] + 
      beta4_JD*rain_prev_wet_post[t]*popsize[t]
    
    logit(phi_JW[t]) <- alpha_phi_JW   + beta1_JW*rain_wet[t] + beta2_JW*tmax_prev_dry_pre[t] + beta3_JW*popsize[t] + 
      beta4_JW*rain_wet[t]*popsize[t]
    
    logit(phi_AD[t]) <- alpha_phi_AD + beta1_AD*rain_prev_wet_post[t] + beta2_AD*tmax_dry[t] + beta3_AD*popsize[t] + 
      beta4_AD*rain_prev_wet_post[t]*popsize[t]
    
    logit(phi_AW[t]) <- alpha_phi_AW   + beta1_AW*rain_wet[t] + beta2_AW*tmax_prev_dry_pre[t] + beta3_AW*popsize[t] + 
      beta4_AW*rain_wet[t]*popsize[t]
    
    logit(p_JD[t]) <- alpha_p_JD + eps_p_JD[year[t]]
    logit(p_JW[t]) <- alpha_p_JW + eps_p_JW[year[t]]
    logit(p_AD[t]) <- alpha_p_AD + eps_p_AD[year[t]]
    logit(p_AW[t]) <- alpha_p_AW + eps_p_AW[year[t]]
    
  }
  
  ### Define state-transition and observation matrices ###
  
  for (i in 1:M){  
    
    ## Define probabilities of state S(t+1) given S(t) ##
    
    for (t in 1:(n.occasions-1)){
      
      # from NE state
      ps[1,i,t,1] <- 1-gamma[t]  # NE to NE
      ps[1,i,t,2] <- gamma[t]    # NE to JD 
      # (in closed populations, individuals can only enter by birth, in our case here the juveniles are born in January then get capture in postwet)
      ps[1,i,t,3] <- 0           # NE to JW
      ps[1,i,t,4] <- 0           # NE to AD
      ps[1,i,t,5] <- 0           # NE to AW
      ps[1,i,t,6] <- 0           # NE to dead
      
      # from JD state
      ps[2,i,t,1] <- 0              # JD to NE
      ps[2,i,t,2] <- 0              # JD to JD
      ps[2,i,t,3] <- phi_JD[t]      # JD to JW
      ps[2,i,t,4] <- 0              # JD to AD
      ps[2,i,t,5] <- 0              # JD to AW
      ps[2,i,t,6] <- 1-phi_JD[t]    # JD to dead
      
      # from JW state
      ps[3,i,t,1] <- 0              # JW to NE
      ps[3,i,t,2] <- 0              # JW to JD
      ps[3,i,t,3] <- 0              # JW to JW
      ps[3,i,t,4] <- phi_JW[t]      # JW to AD
      ps[3,i,t,5] <- 0              # JW to AW
      ps[3,i,t,6] <- 1-phi_JW[t]    # JW to dead
      
      # from AD state
      ps[4,i,t,1] <- 0              # AD to NE
      ps[4,i,t,2] <- 0              # AD to JD
      ps[4,i,t,3] <- 0              # AD to JW
      ps[4,i,t,4] <- 0              # AD to AD
      ps[4,i,t,5] <- phi_AD[t]      # AD to AW
      ps[4,i,t,6] <- 1-phi_AD[t]    # AD to dead
      
      # from AW state
      ps[5,i,t,1] <- 0              # AW to NE
      ps[5,i,t,2] <- 0              # AW to JD
      ps[5,i,t,3] <- 0              # AW to JW
      ps[5,i,t,4] <- phi_AW[t]      # AW to AD
      ps[5,i,t,5] <- 0              # AW to AW
      ps[5,i,t,6] <- 1-phi_AW[t]    # AW to dead
      
      # from dead state
      ps[6,i,t,1] <- 0              # dead to NE
      ps[6,i,t,2] <- 0              # dead to JD
      ps[6,i,t,3] <- 0              # dead to JW
      ps[6,i,t,4] <- 0              # dead to AD
      ps[6,i,t,5] <- 0              # dead to AW
      ps[6,i,t,6] <- 1              # dead to dead
      
      ## Define probabilities of O(t) given S(t) ##
      
      
      # NE state
      po[1,i,t,1] <- 0              # seen as JD
      po[1,i,t,2] <- 0              # seen as JW
      po[1,i,t,3] <- 0              # seen as AD
      po[1,i,t,4] <- 0              # seen as AW
      po[1,i,t,5] <- 1              # not seen
      
      # JD state
      po[2,i,t,1] <- p_JD[t]        # seen as JD
      po[2,i,t,2] <- 0              # seen as JW
      po[2,i,t,3] <- 0              # seen as AD
      po[2,i,t,4] <- 0              # seen as AW
      po[2,i,t,5] <- 1-p_JD[t]      # not seen
      
      # JW state
      po[3,i,t,1] <- 0              # seen as JD
      po[3,i,t,2] <- p_JW[t]        # seen as JW
      po[3,i,t,3] <- 0              # seen as AD
      po[3,i,t,4] <- 0              # seen as AW
      po[3,i,t,5] <- 1-p_JW[t]      # not seen
      
      # AD state
      po[4,i,t,1] <- 0              # seen as JD
      po[4,i,t,2] <- 0              # seen as JW
      po[4,i,t,3] <- p_AD[t]        # seen as AD
      po[4,i,t,4] <- 0              # seen as AW
      po[4,i,t,5] <- 1-p_AD[t]      # not seen
      
      # AW state
      po[5,i,t,1] <- 0              # seen as JD
      po[5,i,t,2] <- 0              # seen as JW
      po[5,i,t,3] <- 0              # seen as AD
      po[5,i,t,4] <- p_AW[t]        # seen as AW
      po[5,i,t,5] <- 1-p_AW[t]      # not seen
      
      
      # dead state
      po[6,i,t,1] <- 0              # seen as JD
      po[6,i,t,2] <- 0              # seen as JW
      po[6,i,t,3] <- 0              # seen as AD
      po[6,i,t,4] <- 0              # seen as AW
      po[6,i,t,5] <- 1              # not seen
      
    } # t
    #} # i
    
    
    ### LIKELIHOOD ##
    
    #for (i in 1:M){
    
    ## Define latent state at first occasion ##
    z[i,1] <- 1   # Make sure that all M individuals are in state 1 at t=1
    
    for (t in 2:n.occasions){
      
      ## State process: draw S(t) given S(t-1) ##
      z[i,t] ~ dcat(ps[z[i,t-1], i, t-1, 1:6])
      
      ## Observation process: draw O(t) given S(t) ##
      y[i,t] ~ dcat(po[z[i,t], i, t-1, 1:5])
      
      
      ### DERIVED PARAMETERS ##
      
      # total JD
      z_JD[i, t-1] <- equals(z[i,t], 2) 
      
      # total JW
      z_JW[i, t-1] <- equals(z[i,t], 3) 
      
      # total AD
      z_AD[i, t-1] <- equals(z[i,t], 4) 
      
      # total AW
      z_AW[i, t-1] <- equals(z[i,t], 5) 
      
    } # t
    
    
  } # i
  
  for (t in 2:n.occasions){
    N_JD[t-1] <- sum(z_JD[1:M,t-1])
    N_JW[t-1] <- sum(z_JW[1:M,t-1])
    N_AD[t-1] <- sum(z_AD[1:M,t-1])
    N_AW[t-1] <- sum(z_AW[1:M,t-1])
  }
})

# prepare initial values for z (true state):
head(mlemur_inits)

# ADDING DUMMY OCCASION -----
#############################

# Add dummy occasion to initial values: 1992_2 is the dummy occasion
# the dummy has to be 1 = not yet entered
my.z.init <- cbind(rep(1, dim(mlemur_inits)[1]), mlemur_inits)

nz <- 200 # 200 additional individuals for 100 individuals that were sampled 
my.z.init.ms <- rbind(my.z.init, matrix(0, ncol = dim(my.z.init)[2], nrow = nz))
my.z.init.ms[my.z.init.ms==0] <- 1 # fake new ones are not entered
str(my.z.init.ms)


# Add dummy occasion to the data: a column of "not seen" (5) before the first occasion
CH.du <- cbind(rep(5, dim(mlemur_ms_CH)[1]), mlemur_ms_CH)
str(CH.du)

# Augment data with fake individuals, add 5 for "not seen" again
nz <- 200
CH.ms <- rbind(CH.du, matrix(5, ncol = dim(CH.du)[2], nrow = nz))
str(CH.ms)

# for the nimble function for the gamma values
postwet_season <- (1:(dim(CH.ms)[2]-1))%%2 # we give 1 to postwet season
length(postwet_season)

# Bundle data
str(nimble_constants <- list(n.occasions = dim(CH.ms)[2], M = dim(CH.ms)[1], postwet_season = postwet_season, 
                             n.years = length(unique(year$yearCat)), year = year$yearCat,
                             rain_wet = rain_wet$prec_std, 
                             rain_prev_wet_post = rain_prev_wet_post$prec_std,
                             tmax_dry = tmax_dry$tmax_std, 
                             tmax_prev_dry_pre= tmax_prev_dry_pre$tmax_std,
                             popsize = popsize$pop_std))

str(nimble_data <- list(y = CH.ms))


set.seed(123) # for reproducibility

# create initial values for gamma
postwet_index = seq(1,(dim(CH.ms)[2]-1),2)
n.postwetindex = length(postwet_index)

gamma_init <- rep(0, (dim(CH.ms)[2]-1))
for (t in postwet_index[1:n.postwetindex]){
  # they can only enter during wet seasons
  gamma_init[t] <- runif(1, 0, 1)
  
}

inits <- function(){list(alpha_phi_JD = runif(1, 0, 1), alpha_phi_JW = runif(1, 0, 1), 
                         alpha_phi_AD = runif(1, 0, 1), alpha_phi_AW = runif(1, 0, 1), 
                         alpha_p_JD = runif(1, 0, 1), alpha_p_JW = runif(1, 0, 1), 
                         alpha_p_AD = runif(1, 0, 1), alpha_p_AW = runif(1, 0, 1),
                         
                         gamma.prior = gamma_init,
                         
                         beta1_JD = 0, beta2_JD = 0, beta3_JD = 0,
                         
                         beta1_JW = 0, beta2_JW = 0, beta3_JW = 0,
                         
                         beta1_AD = 0, beta2_AD = 0, beta3_AD = 0,
                         
                         beta1_AW = 0, beta2_AW = 0, beta3_AW = 0,
                         
                         sigma_p_JD = runif(1, 0, 0.5), sigma_p_JW = runif(1, 0, 0.5), 
                         sigma_p_AD = runif(1, 0, 0.5), sigma_p_AW = runif(1, 0, 0.5), 
                         
                         eps_p_JD = rep(0,length(unique(year$yearCat))), eps_p_JW = rep(0,length(unique(year$yearCat))),
                         eps_p_AD = rep(0,length(unique(year$yearCat))), eps_p_AW = rep(0,length(unique(year$yearCat))),
                         
                         z = my.z.init.ms)}

params <- c("alpha_phi_JD", "alpha_phi_JW", "alpha_phi_AD", "alpha_phi_AW",
            "alpha_p_JD", "alpha_p_JW", "alpha_p_AD", "alpha_p_AW", "gamma",
            "beta1_JD", "beta2_JD", "beta3_JD",
            "beta1_JW", "beta2_JW", "beta3_JW",
            "beta1_AD", "beta2_AD", "beta3_AD",
            "beta1_AW", "beta2_AW", "beta3_AW",
            "eps_p_JD", "eps_p_JW", "eps_p_AD", "eps_p_AW",
            "N_JD", "N_JW", "N_AD", "N_AW")

start1 <- Sys.time()

# Compile model
model_mlemur_interactionI_0923 <- nimbleModel(code = mlemur_interactionI_0923, constants = nimble_constants, data = nimble_data,inits = inits())

end1 <- Sys.time()

end1-start1

start2 <- Sys.time()

# Run model
run_mlemur_interactionI_0923 <- nimbleMCMC(model = model_mlemur_interactionI_0923,
                                           monitors = params,
                                           thin = 1, niter = 100000, nburnin = 0, nchains = 2, samplesAsCodaMCMC = TRUE, WAIC = TRUE)

end2 <- Sys.time()

end2-start2

# not putting any burn-in to see where it starts to converge
# not putting any thinning to get all samples
# trial MCMC settings are:  thin = 1, niter = 10, nburnin = 0, nchains = 2
# running MCMC settings are:  thin = 1, niter = 100000, nburnin = 0, nchains = 2


str(run_mlemur_interactionI_0923$samples)

