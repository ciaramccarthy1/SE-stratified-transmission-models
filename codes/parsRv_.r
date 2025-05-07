pars <- list()
pars <- within(pars, {
  
    Disease     <- "RSV-illness"
    Vaccination <- "Yes"
    Incidence   <- pset$Incidence
    
    #Clinical responses
    #Susceptibility - Secondary infection (relative to primary) Hodgson 2020
    # - don't use as not modelling sequential exposures over years (7 years of historical data)
    #u <- rep(1,9) #1ry exposure #rep(0.89,9), rep(0.81,9), rep(0.33,9) #2ry-4th exposure #(relative to primary)
    #Susceptibility - age-adjusted Henderson 1979, Waterlow 2021
    u   <- c(0.85, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65)
    #Critically infected fraction - age-adjusted Hodgson 2020 - SRV_parameter_age-adjustment_20mar25.r
    y   <- c(0.8656, 0.4840, 0.3486, 0.2470, 0.2470, 0.2470, 0.2470, 0.2470, 0.2470)
    #Clinical fraction - by age and IMD group
    y45 <- rep(y,5)
    
    #Hospitalisation fraction
    #h    <- 
    #Mortality fraction (in hospital)
    #m    <- 

    age = c(mean(0:4),mean(5:11),mean(12:17),mean(18:29),mean(30:39),mean(40:49),mean(50:59),mean(60:69),mean(70:90)) #85))

    #Mortality fraction if clinically infected - derived and age-adjusted from IFR/y in Hodgson 2020
    m <- c(0.002609698, 0.001407748, 0.165154365, 0.405802227, 0.405802227, 0.405802227, 0.405802227, 0.405802227,
    0.405802227) #c(0.1,rep(0,8)
    #temporal
    dt     <- 0.1             #0.01 #time step (days) #smaller than Baguelin 2013 (0.25)
    times  <- 0:180 #365      #days sequence
    nt     <- (max(times)-min(times))/dt + 1       #no. time points, iterations
    nw     <- ceiling((max(times)-min(times))/7)   #weeks length of model run
    nd     <- ceiling((max(times)-min(times)))+1   #days length of model run
    
    #demography
    ages   <- c("0 to 4","5 to 11","12 to 17","18 to 29","30 to 39","40 to 49","50 to 59","60 to 69", "70+")
    na     <- 9               #number of age groups
    nimd   <- 5               #number of SE groups
    urban  <- T               #area: urban (T), rural (F)
    ageons <- c(0.0466, 0.0873, 0.0693, 0.14997, 0.1337, 0.1258, 0.1351, 0.1058, 0.1358); ageons=ageons/sum(ageons) #2020 mid
    
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
    #rC     <-                #rate of loss of positivity

    R0     <- 2.8             #Reis and Sharma 2018 (but simpler model, not nec consitent with other parameters)
    R0     <- 4.5             #lit rev within Reis and Sharma 2018
    R0     <- 4.8             #R0() consistent with beta below
    f      <- 0.634           #relative transmission of U group - Hodgson 2020
    beta   <- 0.0972          #probability of infectivity upon contact - Hodgson 2020
    
    #Initial condition
    #imd=1, age 30 to 39", 1/100,000 latent infections
    pE1g0   <- rep(0,na*nimd)  #initialise proportion latently infected across age x SES groups
    pE1g0[5] = (1/10^5)        #1/100,000 latent infections in age group 5 in SES 1
    
    #rate of reporting by age, Hodgson 2020
    rrep <- c(0.0023321, 0.0000305, 0.0000305, 0.0000305, 0.0000305, 0.0000305, 0.0000888, 0.000147, 0.000147)

    #vaccines
    ve   <- rep(0.5, na)        #efficacy 0.5
    veff <- rep(ve,  nimd)
    vc   <- rep(1, na)          #coverage 0.1
    vcov <- rep(vc,  nimd)
    vcln <- 0.15                #reduction in clinical fraction
    rV   <- 1/180               #rate of immunisation
    
    #NB parameters
    #k      <- 1              #dispersion/shape par of NB likelihood - #derived sh fit

})


