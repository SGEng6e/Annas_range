
# Nimble model code for:
# Supplemental feeding as a driver of population expansion and morphological change in Anna’s Hummingbirds

# NM Alexandre, FG Romero, SG English, E Grames, F Garzón-Agudelo, K Epperly, T Barnes, DR Powers, 
## AE Smith, Z Migicovsky, L Stein, S Akalu, H Sridhar, G Montross, E Collins, A Rico-Guevara

##### LIBRARIES #############################################################################################

library(nimble)
library(nimbleHMC)

##### LOAD DATA #############################################################################################

load("./nimble_data.rda")
load("./nimble_constants.rda")

##### MODEL SCRIPT ##########################################################################################

code <- nimbleCode({
    
  for (i in 1:n.count){
    ## LINEAR PROCESS MODEL ##
    log(lam[i]) <- 
      beta_in +
      beta_fa * FA[i] +
      beta_ea * EA[i] +
      beta_hp * HP[i] +
      beta_lt * LT[i] * YR[i] +
      beta_yr * YR[i] +
      beta_ltln[route[i]] +
      beta_ef * EF[i] +
      eps[i]
    
    ### overdispersion
    eps[i] ~ dnorm(0, sd = sd_ovds)
    
    ## LIKELIHOOD ##
    count[i] ~ dpois(lam[i])
    sim_count[i] ~ dpois(lam[i])
  } # i
  
  ##################################### PRIORS ##############################################################
  
  beta_in ~ dnorm(0, sd = 2)
  beta_fa ~ dnorm(0, sd = 2)
  beta_ea ~ dnorm(0, sd = 2)
  beta_hp ~ dnorm(0, sd = 2)
  beta_yr ~ dnorm(0, sd = 2)
  beta_lt ~ dnorm(0, sd = 2)
  beta_ef ~ dnorm(0, sd = 2)
  sd_ovds ~ dunif(0, 10)
  
  ##################################### SPATIAL GAM EFFECTS #################################################
  
  sd_ltln ~ dgamma(1,1)
  
  for (k in 1:n.knots) {
    BX_ltln[k] ~ dnorm(0, sd = sd_ltln)
  }
  
  ##################################### DERIVED / MARGINAL EFFECTS ##########################################
  for (i in 1:n.route) {
    beta_ltln[i] <- inprod(BX_ltln[1:n.knots], latlon[i,1:n.knots])
  }
  
  baseline[1:n.route] <- exp(beta_in + beta_ltln[1:n.route])
  marginal_fa[1:n.count] <- exp(beta_in + beta_fa * FA[1:n.count])
  marginal_ea[1:n.count] <- exp(beta_in + beta_ea * EA[1:n.count])
  marginal_yr[1:n.count] <- exp(beta_in + beta_yr * YR[1:n.count] + beta_lt * LT[1:n.count] * YR[1:n.count])
  marginal_hp[1:n.count] <- exp(beta_in + beta_hp * HP[1:n.count])
  
  ### Freeman–Tukey residuals
  Rn1i[1:n.count] <- pow(sqrt(count[1:n.count]) - sqrt(lam[1:n.count]), 2)
  Rn2i[1:n.count] <- pow(sqrt(sim_count[1:n.count]) - sqrt(lam[1:n.count]), 2)
  
  Rn1 <- sum(Rn1i[1:n.count])
  Rn2 <- sum(Rn2i[1:n.count])
  # Bayesian p-value
  bayes.p <- Rn2 > Rn1
})

##### MODEL PARAMS ##########################################################################################

params <- c("beta_in","beta_fa","beta_ea","beta_yr","beta_lt","beta_hp","beta_ef","BX_ltln","sd_ltln","sd_ovds",
            "baseline","marginal_fa","marginal_ea","marginal_yr","marginal_hp","bayes.p","Rn2","Rn1")

inits <- list(
  beta_in = 4,
  beta_fa = 0,
  beta_ea = 0,
  beta_hp = 0,
  beta_lt = 0,
  beta_yr = 0,
  beta_ef = 0.5,
  BX_ltln = rep(0, nimble.constants$n.knots),
  eps = rep(0, nimble.constants$n.count),
  sd_ovds = 1,
  sd_ltln = 1,
  sim_count = round(nimble.data$count * 1.3)
)

##### CONFIGURE AND RUN MODEL ###############################################################################

Rmodel <- nimbleModel(code, nimble.constants, nimble.data, inits, buildDerivs = T)
config <- configureHMC(Rmodel, monitors = params)
Rmcmc <- buildMCMC(config)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
samples <- runMCMC(
  Cmcmc,
  nchains = 1,
  niter = 20000,
  nburnin= 5000,
  thin = 2,
  samplesAsCodaMCMC = T
)

idx <- length(list.files("./", pattern = "anna_glmnuts_chain"))+1
assign(paste0("anna_glmnuts_chain", idx), value = samples)

save(list = paste0("anna_glmnuts_chain", idx), file = paste0("./anna_glmnuts_chain", idx, ".rda"))
