pars <- list()
pars <- within(pars, {
  
    Disease <- "COVID-19"
    Vaccination <- "Yes"
    Incidence   <- pset$Incidence
    
    #Clinical responses
    #Susceptibility
    # TODO(10-age scaffold): 60-64 inherits old 60-69; 65-74 mean(old 60-69, 70+); 75+ inherits old 70+
    u    <- c(0.40, 0.39, 0.38, 0.72, 0.86, 0.80, 0.82, 0.88, (0.88+0.74)/2, 0.74) #age-adjusted Davies Nat Med 2020
    #Critically infected fraction (initial values; overridden below by exp formula)
    # TODO(10-age scaffold): same as above
    y    <- c(0.29, 0.27, 0.21, 0.26, 0.33, 0.40, 0.49, 0.63, (0.63+0.69)/2, 0.69) #age-adjusted Davies Nat Med 2020
    #sh estimates
    age = c(mean(0:4),mean(5:14),mean(15:19),mean(20:29),mean(30:39),mean(40:49),mean(50:59),mean(60:64),mean(65:74),mean(75:90))
    yA_0=0.438 #Med: 0.438, CI: [0.305,0.663]
    yr_0=0.012 #Med: 0.012, CI: [0.007,0.019]
    y[3:10]=yA_0*exp((age[3:10]-age[10])*yr_0)
    #Clinical fraction - by age and IMD group
    y45  <- rep(y,5)
    #TODO: Re-weight u y h when have final demography

    #Hospitalisation fraction
    #h    <-
    #Mortality fraction (in hospital)
    #m    <-
    #Mortality fraction if clinically infected (derived from Verity 2020 IFR and sh est)
    # TODO(10-age scaffold): 60-64 inherits old 60-69; 65-74 mean(old 60-69, 70+); 75+ inherits old 70+
    m <- c(0.000056, 0.000116, 0.000348, 0.001210, 0.003327, 0.005627, 0.018450, 0.053066, (0.053066+0.139813)/2, 0.139813)

    #Hospitalisation pathway (H compartment): I2 -> H at rate h*rI2R; H -> D at rate mH*rH; H -> R at rate (1-mH)*rH
    # TODO(H scaffold): placeholder h = pmin(10*m, 0.8); mH = m/h so total mortality given clinical preserved.
    h  <- pmin(10*m, 0.8)
    mH <- m / h
    rH <- 1/6                     #1/hospital stay length (days^-1); TODO use COVID-specific value
    
    #temporal
    dt     <- 0.1             #time step (days)
    times  <- 0:180 #365      #days sequence
    nt     <- (max(times)-min(times))/dt + 1       #no. time points, iterations
    nw     <- ceiling((max(times)-min(times))/7)   #week length of model run
    nd     <- ceiling((max(times)-min(times)))+1   #days length of model run
    
    #demography
    ages   <- c("0 to 4","5 to 14","15 to 19","20 to 29","30 to 39","40 to 49","50 to 59","60 to 64","65 to 74","75+")
    na     <- 10              #number of age groups
    nimd   <- 5               #number of SE groups
    urban  <- T               #area: urban (T), rural (F)
    # ageons (age proportions) is now computed in modelrun.r from demographics_10age.csv
    
    #natural history            
    rEI    <- 1/3             #latency = rEU, Davies 2020 Nat Med
    rI1I2  <- 1/2.1           #recovery, Davies 2020 Nat Med
    rI2R   <- 1/2.9           #recovery, Davies 2020 Nat Med
    rIR    <- 1/(1/rI1I2 + 1/rI2R)  #1/3
    rUR    <- 1/5 #rIR        #recovery, Davies 2020 Nat Med
    #rIH    <-                #hospitalisation
    #rHR    <-                #recovery rate in hospital
    #rHD    <-                #death rate in hospital

    R0     <- 2.5             #Assumed, close to Knock 2021 and Davies 2020
    f      <- 0.5             #relative transmission of U group, Davies 2020 Nat Med
    beta   <- 0.06            #transmission prob upon contact, LG 2023

    #Initial condition
    pE1g0   <- rep(0,na*nimd) #initialise proportion latently infected across age x SES groups
    pE1g0[5] = (1/10^5)       #1/100,000 latent infections in age group 5 in SES 1 (age 30 to 39, imd=1)
    
    #rate of reporting by age #Assumed
    rrep <- rep(0.5, na)

    #vaccines: leaky V compartment with breakthrough shadow chain (V/Ev/Iv/Uv/Hv/Rv/Dv)
    # TODO(V scaffold): all VEs uniform across (age, IMD) at placeholder values; review for COVID programme.
    vc      <- rep(1, na)
    vcov    <- rep(vc, nimd)
    VE_inf  <- rep(0.0, na*nimd)         #placeholder: no infection blocking
    VE_sym  <- rep(0.5, na*nimd)         #placeholder
    VE_hosp <- rep(0.7, na*nimd)         #placeholder
    VE_mort  <- rep(0.5, na*nimd)         #placeholder
    rV      <- 1/180                     #rate of immunisation (per day)
    rW      <- 1/365                     #vaccine waning rate (1/year); TODO source-paper value
    rW_nat  <- 0                         #natural waning rate; default 0
    
})


