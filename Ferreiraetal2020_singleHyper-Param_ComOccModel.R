###################################################################################################

# Multi-species occupancy models to investigate effect of habitat protection on Cerrado mammals #
# Ferreira et al 2020, Biological Conservation https://doi.org/10.1016/j.biocon.2020.108762

# model with single hyper-parameter; used for spp-level inference and species richness

# Based largely on Zipkin's et al 2010 code with inputs from Kery & Royle 2016
# Earlier versions influenced by Kery & Royle 2016 (Chap 11 in AHME)

###################################################################################################

##### Loading and preparing data #####

# Step 1 - load data and construct 3D array (sites X occasions X species)
# loading 1/0 matrix for 501 sites (rows) stacked by 34 spp (27 recorded + 7 augmented)
# columns are occasions (12)
# Dasypus novemncinctus, Dasypus septemcinctus and Dasypus sp joined under the latter

tmp2 <- read.csv("Augmented_34spp_SVP-2012-2017_6dMatrix_NAiflessThan6d.csv") 
tmp2 <- tmp2[,4:15] #getting only det/nondet columns


# dividing original matrix into a list of several matrices for each species
# 501 is the number of sites by which the mtx should be divided
Mtxlist <- lapply(split(seq_len(nrow(tmp2)),(seq_len(nrow(tmp2))-1) %/%501 +1), 
                  function(i) tmp2[i,])

# using abind (package that joins matrices into arrays) to create a 3D array based on the list of matrices above: "Mtxlst"
# dimensions are: sites, occasions(rep), species
library(abind)
X <- abind(Mtxlist, along=3)
str(X)
X[1:20,,1] # just checking the array

n <- dim(X)[3]-7 # number of observed spp (total n of species minus the number of augmented species)
nzeroes <- dim(X)[3]-27 # number of augmented (all 0) species in the data set
J <- dim(X)[1] # number of CT sites
K <- read.csv("validOccasions_VECTOR_NAiflessThan6d.csv") # number of occasions in each CT site

# Step 2 - load covariates data
covs <- read.csv("Scaled_SVP_sitecovs_501sites.csv", header=T) # continuous variables are scaled
head(covs)

pa_type <- covs$pa_type # one PA type must be 0 to eliminate its term in the occurence equation
NDVI <- covs$NDVImean_500m
trail <- covs$trail
Dist_road <- covs$Dist_road
Dist_water <- covs$Dist_water
CTcode2 <- as.numeric(covs$CTcode2) #as.numeric works here because none of the CTmodels must be 0

# Step 3 - Bundle and summarize data

str(sp.data <- list(n = n, nzeroes = nzeroes, J = J, K = K, X = X,
                    pa_type = pa_type, 
                    Dist_road = Dist_road, Dist_water = Dist_water, NDVI = NDVI,
                    trail = trail, CTcode2 = CTcode2) )


# Step 4 - Defining initial values for the MCMC

wst <- rep(1, n+nzeroes)                   # Simply set everybody at occurring
zst <- array(1, dim = c(J, n+nzeroes)) # ditto
# wst and zst are as suggested by Kery & Royle 2015; this is diff from Zipkin et al 2010

sp.inits <- function() {
  omegaGuess = runif(1, n/(n+nzeroes), 1)
  psi.meanGuess = runif(1, .25,1)
  list(omega=omegaGuess, Z = zst, w = wst, 
       psiMulti = rnorm(n = n+nzeroes), psiStrict = rnorm(n = n+nzeroes), 
       alphaRoad = rnorm(n = n+nzeroes), 
       alphaWat = rnorm(n = n+nzeroes), alphaNDVI = rnorm(n = n+nzeroes), 
       betaNtrail = rnorm(n = n+nzeroes), betalYtrail = rnorm(n = n+nzeroes), 
       betalCT2 = rnorm(n = n+nzeroes), betalCT3 = rnorm(n = n+nzeroes))
}
#### Multi-spp model, 1 Hyper-parameter only ####
sink("modelB_PostViva.txt")
cat("
    model {
    
    # Prior distribution for community-level params - hyperpriors
    omega ~ dunif(0,1)   # 'hyper-psi'
    
    mu.Multi  ~ dnorm(0, 0.001)
    mu.Strict ~ dnorm(0, 0.001)
    mu.alphaRoad ~ dnorm(0, 0.001)
    mu.alphaWat ~ dnorm(0, 0.001)
    mu.alphaNDVI ~ dnorm(0, 0.001)
    mu.betaNtrail ~ dnorm(0, 0.001)
    mu.betaYtrail ~ dnorm(0, 0.001)
    mu.betaCT2 ~ dnorm(0, 0.001)
    mu.betaCT3 ~ dnorm(0, 0.001)
    
    # tau and sd based on Kery & Royle; Zipkin specify tau directly, without sd
    
    tau.Multi <- pow(sd.Multi,-2) 
    tau.Strict <- pow(sd.Strict,-2) 
    tau.alphaRoad <- pow(sd.alphaRoad,-2) 
    tau.alphaWat <- pow(sd.alphaWat,-2) 
    tau.alphaNDVI <- pow(sd.alphaNDVI,-2) 
    tau.betaNtrail <- pow(sd.betaNtrail,-2) 
    tau.betaYtrail <- pow(sd.betaYtrail,-2) 
    tau.betaCT2 <- pow(sd.betaCT2,-2) 
    tau.betaCT3 <- pow(sd.betaCT3,-2) 
    
    sd.Multi ~ dunif(0,2) 
    sd.Strict ~ dunif(0,2) 
    sd.alphaRoad ~ dunif(0,2) 
    sd.alphaWat ~ dunif(0,2) 
    sd.alphaNDVI ~ dunif(0,2) 
    sd.betaNtrail ~ dunif(0,2) 
    sd.betaYtrail ~ dunif(0,2) 
    sd.betaCT2 ~ dunif(0,2) 
    sd.betaCT3 ~ dunif(0,2) 
    
    
    # specify species-level priors for species i (out of 34) influenced by comm-level priors
     for (i in 1:(n+nzeroes)) {
    
    w[i] ~ dbern(omega)   # Superpopulation process: Ntotal species sampled out of all spp available
    psiMulti[i] ~ dnorm(mu.Multi, tau.Multi)
    psiStrict[i] ~ dnorm(mu.Strict, tau.Strict)
    alphaRoad[i] ~ dnorm(mu.alphaRoad, tau.alphaRoad)
    alphaWat[i] ~ dnorm(mu.alphaWat, tau.alphaWat)
    alphaNDVI[i] ~ dnorm(mu.alphaNDVI, tau.alphaNDVI)
    
    betaNtrail[i] ~ dnorm(mu.betaNtrail, tau.betaNtrail)
    betaYtrail[i] ~ dnorm(mu.betaYtrail, tau.betaYtrail)
    betaCT2[i] ~ dnorm(mu.betaCT2, tau.betaCT2)
    betaCT3[i] ~ dnorm(mu.betaCT3, tau.betaCT3)
    
    
    
    # Ecological model for true occurrence (process model) of sp i at site j
    # loop to define Z-matrix ('true' matrix of 1-0), obs: spp loop above still open
    
    for (j in 1:J) {
    logit(psi[j,i]) <- psiMulti[i]*(1 - pa_type[j]) + psiStrict[i]*pa_type[j] +
    alphaRoad[i]*Dist_road[j] + alphaWat[i]*Dist_water[j] + alphaNDVI[i]*NDVI[j]
    
    mu.psi[j,i] <- psi[j,i] * w[i]
    Z[j,i] ~ dbern(mu.psi[j,i])
    
    
    
    # Observation model for replicated detection/nondetection observations; observed 1-0 matrix
    # detetection of species i at site j for survey occasion k; obs spp and site loops still open
    
    for (k in 1:K[j]) {
    logit(p[j,k,i]) <- betaNtrail[i] + betaYtrail[i] * trail[j] + 
    betaCT2[i] * equals(CTcode2[j],2) + betaCT3[i] * equals(CTcode2[j],3)
    
    mu.p[j,k,i] <- p[j,k,i] * Z[j,i]
    X[j,k,i] ~ dbern(mu.p[j,k,i])
    Xnew[j,k,i] ~ dbern(mu.p[j,k,i])  # replicate data
    
    # Observed dataset
    chi2.actual[j,k,i] <- pow(X[j,k,i] - mu.p[j,k,i], 2)/ (mu.p[j,k,i] + 0.0001)  # Add small value to denominator to prevent division by zero
    
    # Expected dataset  
    chi2.sim[j,k,i] <- pow(Xnew[j,k,i] - mu.p[j,k,i], 2)/ (mu.p[j,k,i] + 0.0001)  # Add small value to denominator to prevent division by zero
    
    } #  temporal replicate loop
    
    chi2.actual.sum[j,i] <- sum(chi2.actual[j,1:K[j],i])
    chi2.sim.sum[j,i] <- sum(chi2.sim[j,1:K[j],i])  
    
    } #  site loop
    } # close species loop
    
    
    # Calculate the discrepancy measure (Pearsons chi-squared) for each spp, which is then defined as the mean(p.fit > p.fitnew)
    # note that this loop does not include the all-zero spp
    for(g in 1:n){
    fit.sp.actual[g] <- sum(chi2.actual.sum[,g])      # Species-specific fit statistic for actual dataset
    fit.sp.sim[g] <- sum(chi2.sim.sum[,g])            # Species-specific fit statistic for simulated dataset
    
    c.hat.sp[g] <- fit.sp.actual[g]/fit.sp.sim[g]
    bpv.sp[g] <- step(fit.sp.sim[g] - fit.sp.actual[g])
    }
    
    
    # Calculate overall Pearson's chi-squared discrepency measure, defined post-hoc as: mean(p.fit>p.fit.new)
    
    fit.actual <- sum(chi2.actual.sum[1:J, 1:n])
    fit.sim <- sum(chi2.sim.sum[1:J, 1:n])
    
    c.hat <- fit.actual/fit.sim
    
    bpv <- step(fit.sim - fit.actual)    
    
    ###### Derived quantities ########
    
    for(o in 1:n){
    spp.Pvalue[o] <- mean(fit.sp.actual[o] > fit.sp.sim[o])
    }
    overall.Pvalue <- mean(fit.actual > fit.sim)
    
    # calculating effect of protection: psi in Strict PA - psi in multi PA
    for(u in 1:n){
    PAeffect[u] <- psiStrict[u] - psiMulti[u]
    }
    
    #### estimates of spp richness ###
    #
    for(q in 1:(n+nzeroes)){
    Nocc.fs[q] <- sum(Z[,q])       # Number of occupied sites among the 501
    Nocc.fs.PAst[q] <- sum(Z[c(1:38,108:166,224:331,454:501),q])  #number occupied sitesin strict PA
    Nocc.fs.PAmu[q] <- sum(Z[c(39:107,167:223,332:453),q])  #number of occupied sites in multi-use PA
    }
    
    # Site species richness overall, large (>15kg) and threatened (according to MMA 2014)     
    for (s in 1:J){
    Nsite[s] <- sum(Z[s,])          # Number of occurring species at each site 
    Nsite.large[s] <- sum(Z[s,c(1,11,16,17,19,21,22,25,26,29,30,32,33)])      # Number of >15kg spp at each site
    Nsite.small[s] <- sum(Z[s,c(2:10,12:15,18,20,23,24,27,28,31,34)]) # number of spp < 15kg
    Nsite.threat[s] <- sum(Z[s,c(1,12,14,15,17,19,21,23,24,25,26,28,29,30,32,34)]) # Numb of threat spp (MMA2014) at each site
    Nsite.Glob.threat[s] <- sum(Z[s,c(14,17,25,26,28,29,32)]) # Numb of globally spp (IUCN) at each site
    Nsite.non.threat[s] <- sum(Z[s,c(2:11,13,16,18,20,22,27,31,33)]) #number of spp not present in either national or global redlist
    }
    
    ## Mean site spp richness for each type of PA
    # All species
    Nsite.strict <- Nsite[c(1:38,108:166,224:331,454:501)]    # sprich per site for strict PA; no need to monitor
    Nsite.multi <-Nsite[c(39:107,167:223,332:453)]          # sprich per site for multi-use PA; no need to monitor 
    
    mean.Nsite.strict <- mean(Nsite.strict) # mean SPrich per site for strict PAs, param to be monitored
    mean.Nsite.multi <- mean(Nsite.multi)   # mean SPrich per site for multi-use PAs, param to be monitored
    
    
    # Large species
    Nlarge.strict <- Nsite.large[c(1:38,108:166,224:331,454:501)]    # large sprich per site for strict PA; no need to monitor
    Nlarge.multi <-Nsite.large[c(39:107,167:223,332:453)]          # large sprich per site for multi-use PA; no need to monitor 
    
    mean.Nlarge.strict <- mean(Nlarge.strict) # mean large SPrich per site for strict PAs, param to be monitored
    mean.Nlarge.multi <- mean(Nlarge.multi)   # mean large SPrich per site for multi-use PAs, param to be monitored
    
    # Small species
    Nsmall.strict <- Nsite.small[c(1:38,108:166,224:331,454:501)]    # small sprich per site for strict PA; no need to monitor
    Nsmall.multi <-Nsite.small[c(39:107,167:223,332:453)]          # small sprich per site for multi-use PA; no need to monitor 
    
    mean.Nsmall.strict <- mean(Nsmall.strict) # mean small SPrich per site for strict PAs, param to be monitored
    mean.Nsmall.multi <- mean(Nsmall.multi)   # mean small SPrich per site for multi-use PAs, param to be monitored
    
    # Threatened species (according to MMA 2014)
    Nthreat.strict <- Nsite.threat[c(1:38,108:166,224:331,454:501)]    # threat sprich per site for strict PA; no need to monitor
    Nthreat.multi <-Nsite.threat[c(39:107,167:223,332:453)]          # threat sprich per site for multi-use PA; no need to monitor 
    
    mean.Nthreat.strict <- mean(Nthreat.strict) # mean threat SPrich per site for strict PAs, param to be monitored
    mean.Nthreat.multi <- mean(Nthreat.multi)   # mean threat SPrich per site for multi-use PAs, param to be monitored

    # Globally threatened species (according to IUCN)
    NGlobthreat.strict <- Nsite.Glob.threat[c(1:38,108:166,224:331,454:501)]    # threat sprich per site for strict PA; no need to monitor
    NGlobthreat.multi <-Nsite.Glob.threat[c(39:107,167:223,332:453)]          # threat sprich per site for multi-use PA; no need to monitor 
    
    mean.NGlobthreat.strict <- mean(NGlobthreat.strict) # mean threat SPrich per site for strict PAs, param to be monitored
    mean.NGlobthreat.multi <- mean(NGlobthreat.multi)   # mean threat SPrich per site for multi-use PAs, param to be monitored

    # Non-threatened species
    Nnon.threat.strict <- Nsite.non.threat[c(1:38,108:166,224:331,454:501)]    # non-threat sprich per site for strict PA; no need to monitor
    Nnon.threat.multi <-Nsite.non.threat[c(39:107,167:223,332:453)]          # non-threat sprich per site for multi-use PA; no need to monitor 
    
    mean.Nnon.threat.strict <- mean(Nnon.threat.strict) # mean non-threat SPrich per site for strict PAs, param to be monitored
    mean.Nnon.threat.multi <- mean(Nnon.threat.multi)   # mean non-threat SPrich per site for multi-use PAs, param to be monitored
    
    
    Ntotal <- sum(w[])                 # Total metacommunity size
    
    }
    
    ",fill = TRUE)
sink()

# parameters to monitor

params1 <- c("omega", "mu.Multi", "mu.Strict", 
             "mu.alphaRoad", "mu.alphaWat", "mu.alphaNDVI", 
             "mu.betaNtrail", "mu.betaYtrail", "mu.betaCT2", "mu.betaCT3",
             "psiMulti", "psiStrict",
             "alphaRoad", "alphaWat", "alphaNDVI",
             "betaNtrail", "betaYtrail", "betaCT2", "betaCT3",
             "PAeffect",
             "Nocc.fs", "Nocc.fs.PAst", "Nocc.fs.PAmu", 
             "Nsite", "Nsite.large", "Nsite.small", "Nsite.threat", "Nsite.Glob.threat",  "Nsite.non.threat", 
             "mean.Nsite.strict", "mean.Nsite.multi", 
             "mean.Nlarge.strict", "mean.Nlarge.multi",
             "mean.Nsmall.strict", "mean.Nsmall.multi",
             "mean.Nthreat.strict", "mean.Nthreat.multi",
             "mean.NGlobthreat.strict", "mean.NGlobthreat.multi",
             "mean.Nnon.threat.strict", "mean.Nnon.threat.multi",
             "Ntotal", "Omega",
             "c.hat.sp", "c.hat", "bpv.sp", "bpv", 
             "spp.Pvalue", "overall.Pvalue")



# MCMC settings
ni <- 120000   ;   nt <- 10   ;   nb <- 30000   ;   nc <- 3

# Run JAGS, check convergence and summarize posteriors
library(jagsUI)

outB_PostViva <- jags(sp.data, sp.inits, params1, "modelB_PostViva.txt", 
                      n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)



