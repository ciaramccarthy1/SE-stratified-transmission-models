pars <- list()
pars <- within(pars, {
  
  Disease <- "RSV-illness"
  Vaccination <- "No"
  Incidence   <- pset$Incidence
  
  #Clinical responses
  #Age susceptibility = SUSCEPTIBLE-weighted mean of exposure-group susceptibility:
  #   u_a = sum_k d_k S_k^a / sum_k S_k^a   (d_k from David's fitted d1,d2,d3;
  #   S_k^a from his post-burn-in state). 
  u   <- c(0.6392, 0.2614, 0.2408, 0.2401, 0.2400, 0.2400, 0.2400, 0.2400, 0.2400)
  #Symptomatic (clinical) fraction = David rsvie fitted 1-pA (posterior medians), infection-weighted for 0-4.
  y   <- c(0.8476, 0.4776, 0.2564, 0.2564, 0.2564, 0.2564, 0.2564, 0.2564, 0.2564)
  #Clinical fraction - by age and IMD group
  y45 <- rep(y,5)
  
  #Hospitalisation fraction
  #h    <-
  #Mortality fraction (in hospital)
  #m    <-
  
  age = c(mean(0:4),mean(5:14),mean(15:24),mean(25:34),mean(35:44),mean(45:54),mean(55:64),mean(65:74),mean(75:90))
  
  #Hospitalisation & mortality: computed as a POST-PROCESSING step (infections x
  #age-specific outcome risk), exactly as in David's rsvie - NOT within the dynamic model.
  #So the transmission dynamics carry no H/D removal: h = m = mH = 0. This leaves the
  #symptomatic fit unchanged (I still recovers at rIR) and keeps the population conserved
  #under ageing (no disease deaths). Apply outcome risks to the infection output downstream.
  m  <- rep(0, 9)
  h  <- rep(0, 9)
  mH <- rep(0, 9)
  rH <- 1/6                     #1/hospital stay length (days^-1); unused while h = 0
  
  # Natural waning of post-infection immunity (R -> S) is set once with the
  # other vaccine/waning parameters below (see rW_nat near rW).
  
  #temporal
  dt     <- 0.1             #0.01 #time step (days)
  times  <- 0:180 #365      #days sequence
  nt     <- (max(times)-min(times))/dt + 1       #no. time points, iterations
  nw     <- ceiling((max(times)-min(times))/7)   #weeks length of model run
  nd     <- ceiling((max(times)-min(times)))+1   #days length of model run
  
  #demography
  ages   <- c("0 to 4","5 to 14","15 to 24","25 to 34","35 to 44","45 to 54","55 to 64","65 to 74","75+")
  na     <- 9               #number of age groups 
  nimd   <- 5               #number of SE groups
  # ageons (age proportions) is now computed in modelrun.r from demographics_9age.csv
  
  #natural history
  #see also Reis and Sharma 2016, 2018 (consistent parameters, but simpler model)
  rEI    <- 1/4.98          #latency (single; David si=4.93 ~ Hodgson 4.98)
  rI1I2  <- 2/6.16          #legacy Erlang recovery (kept for summary only)
  rI2R   <- 2/6.16          #legacy Erlang recovery (kept for summary only)
  # Age-varying infectious period: infected-weighted mean over David's exposure
  # groups. Group 0 (first infection) clears in 6.10d; groups 2-3 in 4.29d
  # (ga0=6.10, g1=0.880, g2=0.799 from posteriors). Adults (all group 3) -> 4.29d.
  rIR    <- c(0.1961, 0.2314, 0.2330, 0.2330, 0.2330, 0.2330, 0.2330, 0.2330, 0.2330) #I->R rate by age
  rUR    <- c(0.1985, 0.2314, 0.2330, 0.2330, 0.2330, 0.2330, 0.2330, 0.2330, 0.2330) #U->R rate by age
  #rIH    <-                #hospitalisation
  #rHR    <-                #recovery rate in hospital 
  #rHD    <-                #death rate in hospital
  
  R0     <- 2.8             #Reis and Sharma 2018 (but simpler model, not nec consitent with other parameters)
  R0     <- 4.5             #lit rev within Reis and Sharma 2018
  R0     <- 4.8             #R0() consistent with beta below
  f      <- 0.620749        #relative transmission of U group = David rsvie alpha_i (median)
  beta   <- 0.0972          #probability of infectivity upon contact - Hodgson 2020
  
  #Seasonal forcing (Gaussian; from David's rsvie fit)
  # beta(t) = beta*(1 + b1*(1 + exp(-((doy - phi)^2)/(2 psi^2)))),  doy = day-of-year/365
  b1     <- 1.9984          #amplitude
  phi    <- 0.6146          #phase (fraction of year at peak; ~day 224, ~week 32)
  psi    <- 0.2387          #width  (fraction of year)
  
  #Initial condition (disease-agnostic latent seed). NOW UNUSED: modelrun.r
  # overwrites the initial state with David's post-burn-in equilibrium. Kept as
  # zeros only so summary diagnostics still resolve.
  pE1g0   <- rep(0,na*nimd)
  
  #rate of reporting by age, Hodgson 2020
  # TODO PLACEHOLDER older-adult reporting (55-64,65-74,75+) = 0.000147. Not used in symptomatic fit.
  rrep <- c(0.0023321, 0.0000305, 0.0000305, 0.0000305, 0.0000305, 0.0000305, 0.000147, 0.000147, 0.000147)
  
    #vaccines: leaky V compartment with breakthrough shadow chain (V/Ev/Iv/Uv/Hv/Rv/Dv)
    #  rV   - vaccination rate (per day) for fully eligible
    #  vcov - eligibility/coverage fraction per (age x IMD), length ng
    #  VE_inf  - efficacy against infection (reduces FOI on V)
    #  VE_sym  - efficacy against symptoms (reduces E2_v -> I1_v vs U1_v clinical fork)
    #  VE_hosp - efficacy against hospitalisation (reduces I2_v -> H_v)
    #  VE_mort  - efficacy against mortality given hospitalised (reduces H_v -> D)
    #  rW     - vaccine waning rate V -> S (per day)
    #  rW_nat - natural waning rate R/Rv -> S (per day); 0 disables
    # UK RSV 75+ programme: real uptake by IMD quintile (UKHSA Jan 2026 report,
    # population-weighted from deciles by codes/prepare_model_inputs.R).
    # vcov targets band 9 (75+, the last band) only; all other ages have vcov = 0.
    
    ## As a placeholder, from: https://www.sciencedirect.com/science/article/pii/S2666776226000323#appsec1
    VE_hosp_obs <- 0.74
    # https://www.ecdc.europa.eu/en/news-events/rsv-vaccines-safe-and-effective-cochrane-review-finds
    # VE against RSV-associated LRTI - 0.77
    # VE against RSV-associated acute respiratory disease - 0.67
    VE_sym_obs <- 0.67
    VE_mort_obs <- 0.74 # placeholder
    VE_inf_obs <- 0 # placeholder
    
    .uptake <- read.csv(file.path("data", "rsv_uptake_by_imd_quintile.csv"))
    .uptake <- .uptake[order(.uptake$quintile), ]
    stopifnot(.uptake$quintile == 1:5)
    .vc <- matrix(0, na, nimd)
    .vc[na, ] <- .uptake$uptake_pct / 100     # band 9 = last band (75+); fraction in [0,1]
    vcov    <- as.vector(.vc)                 # IMD-major (length ng); col-by-col flatten
    VE_inf  <- rep(0.0, na*nimd)        
    VE_sym  <- rep(1 -  (1 - VE_sym_obs) / (1 - VE_inf_obs), na*nimd) 
    VE_hosp <- rep(1 - (1 - VE_hosp_obs) / (1 - VE_sym_obs), na*nimd) 
    VE_mort  <- rep(1 - (1 - VE_mort_obs) / (1 - VE_hosp_obs), na*nimd)  
    rV      <- 1/180                     #rate of immunisation (per day)
    rW      <- 1/365                     #vaccine waning rate (1/year); TODO source-paper value
    rW_nat  <- 1/358.058                 #natural waning R->S; David's rsvie immunity duration ~358 days (om)
    
})


