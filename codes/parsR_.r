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

    #Mortality fraction if clinically infected - derived and age-adjusted from IFR/y in Hodgson 2020
    # TODO(9-age David-aligned): 15-24 keeps old 15-19 value; 25-34+ = old 20+ adult value. PLACEHOLDER (not used in symptomatic fit).
    m <- c(0.002609698, 0.001407748, 0.165154365, 0.405802227, 0.405802227, 0.405802227,
    0.405802227, 0.405802227, 0.405802227)

    #Hospitalisation pathway (H compartment): I2 -> H at rate h*rI2R; H -> D at rate mH*rH; H -> R at rate (1-mH)*rH
    # TODO(H scaffold): h (hospitalisation fraction given clinical) and mH (mortality given hospitalised) - placeholder derivation:
    #   h  = pmin(10*m, 0.8)      (assume hospitalisation roughly 10x mortality, capped)
    #   mH = m / h                (so overall mortality given clinical = h*mH = m, preserving the legacy m semantics)
    # Replace h and mH with literature values when available (e.g. Hodgson 2020 hospitalisation rates).
    h  <- pmin(10*m, 0.8)
    mH <- m / h
    rH <- 1/6                     #1/hospital stay length (days^-1); TODO use RSV-specific value
    
    #Natural waning of post-infection immunity (R -> S); 0 = lifelong
    # From David's rsvie fit: immunity duration ~358 days (om), so rate = 1/358
    rW_nat <- 1/358.058
	
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

})


