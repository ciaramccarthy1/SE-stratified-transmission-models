pars <- list()
pars <- within(pars, {
  
    Disease <- "RSV-illness"
    Vaccination <- "No"
    Incidence   <- pset$Incidence

    #Clinical responses
    #Susceptibility - Secondary infection (relative to primary) Hodgson 2020
    # - doesn't model exposures sequentially over years (7 years of historical data)
    #Susceptibility - age-adjusted Henderson 1979, Waterlow 2021
    u   <- c(0.85, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65)
    #Critically infected fraction - age-adjusted Hodgson 2020 - age adjusted
    y   <- c(0.8656, 0.4840, 0.3486, 0.2470, 0.2470, 0.2470, 0.2470, 0.2470, 0.2470, 0.2470)
    #Clinical fraction - by age and IMD group
    y45 <- rep(y,5)

    #Hospitalisation fraction
    #h    <-
    #Mortality fraction (in hospital)
    #m    <-

    age = c(mean(0:4),mean(5:14),mean(15:19),mean(20:29),mean(30:39),mean(40:49),mean(50:59),mean(60:64),mean(65:74),mean(75:90))

    #Mortality fraction if clinically infected - derived and age-adjusted from IFR/y in Hodgson 2020
    # TODO(10-age scaffold): 60-64 inherits old 60-69; 65-74 mean(old 60-69, 70+); 75+ inherits old 70+
    m <- c(0.002609698, 0.001407748, 0.165154365, 0.405802227, 0.405802227, 0.405802227, 0.405802227, 0.405802227,
    0.405802227, 0.405802227)

    #Hospitalisation pathway (H compartment): I2 -> H at rate h*rI2R; H -> D at rate mH*rH; H -> R at rate (1-mH)*rH
    # TODO(H scaffold): h (hospitalisation fraction given clinical) and mH (mortality given hospitalised) - placeholder derivation:
    #   h  = pmin(10*m, 0.8)      (assume hospitalisation roughly 10x mortality, capped)
    #   mH = m / h                (so overall mortality given clinical = h*mH = m, preserving the legacy m semantics)
    # Replace h and mH with literature values when available (e.g. Hodgson 2020 hospitalisation rates).
    h  <- pmin(10*m, 0.8)
    mH <- m / h
    rH <- 1/6                     #1/hospital stay length (days^-1); TODO use RSV-specific value
    
    #Natural waning of post-infection immunity (R -> S); 0 = lifelong (default for single-season)
    rW_nat <- 0
	
    #temporal
    dt     <- 0.1             #0.01 #time step (days)
    times  <- 0:180 #365      #days sequence
    nt     <- (max(times)-min(times))/dt + 1       #no. time points, iterations
    nw     <- ceiling((max(times)-min(times))/7)   #weeks length of model run
    nd     <- ceiling((max(times)-min(times)))+1   #days length of model run
    
    #demography
    ages   <- c("0 to 4","5 to 14","15 to 19","20 to 29","30 to 39","40 to 49","50 to 59","60 to 64","65 to 74","75+")
    na     <- 10              #number of age groups
    nimd   <- 5               #number of SE groups
    # ageons (age proportions) is now computed in modelrun.r from demographics_10age.csv
    
    #natural history
    #see also Reis and Sharma 2016, 2018 (consistent parameters, but simpler model)
    rEI    <- 1/4.98          #latency,  Hodgson 2020: first exposure model
    rI1I2  <- 2/6.16          #recovery, Hodgson 2020
    rI2R   <- 2/6.16          #recovery, Hodgson 2020
    rIR    <- 1/(1/rI1I2 + 1/rI2R)
    rUR    <- 1/6.16 #rIR     #recovery, Hodgson 2020
    #rIH    <-                #hospitalisation
    #rHR    <-                #recovery rate in hospital 
    #rHD    <-                #death rate in hospital

    R0     <- 2.8             #Reis and Sharma 2018 (but simpler model, not nec consitent with other parameters)
    R0     <- 4.5             #lit rev within Reis and Sharma 2018
    R0     <- 4.8             #R0() consistent with beta below
    f      <- 0.634           #relative transmission of U group - Hodgson 2020
    beta   <- 0.0972          #probability of infectivity upon contact - Hodgson 2020
    
    #Initial condition
    pE1g0   <- rep(0,na*nimd)  #initialise proportion latently infected across age x SES groups
    pE1g0[5] = (1/10^5)        #1/100,000 latent infections in age group 5 in SES 1 (age 30 to 39, imd=1)
    
    #rate of reporting by age, Hodgson 2020
    # TODO(10-age scaffold): 60-64 inherits old 60-69; 65-74 mean of old 60-69 & 70+; 75+ inherits old 70+
    rrep <- c(0.0023321, 0.0000305, 0.0000305, 0.0000305, 0.0000305, 0.0000305, 0.0000888, 0.000147, 0.000147, 0.000147)

})


