pars <- list()
pars <- within(pars, {
  
    Disease <- "COVID-19"
    Vaccination <- "No"
    Incidence   <- pset$Incidence
    
    #Clinical responses
    #Susceptibility
    u    <- c(0.40, 0.39, 0.38, 0.72, 0.86, 0.80, 0.82, 0.88, 0.74) #age-adjusted Davies Nat Med 2020
    #Critically infected fraction
    y    <- c(0.29, 0.27, 0.21, 0.26, 0.33, 0.40, 0.49, 0.63, 0.69) #age-adjusted Davies Nat Med 2020
    #sh estimates
    age = c(mean(0:4),mean(5:11),mean(12:17),mean(18:29),mean(30:39),mean(40:49),mean(50:59),mean(60:69),mean(70:90))
    yA_0=0.438 #Med: 0.438, CI: [0.305,0.663]
    yr_0=0.012 #Med: 0.012, CI: [0.007,0.019]
    y[3:9]=yA_0*exp((age[3:9]-age[9])*yr_0)                         
    #0.2900 0.2700 0.1996 0.2223 0.2537 0.2861 0.3225 0.3637 0.4380 #derived sh est
    #Clinical fraction - by age and IMD group
    y45  <- rep(y,5)
    #TODO: Re-weight u y h when have final demography
    
    #Hospitalisation fraction
    #h    <- 
    #Mortality fraction (in hospital)
    #m    <- 
    #Mortality fraction if clinically infected (derived from Verity 2020 IFR and sh est)
    m <- c(0.000056, 0.000116, 0.000348, 0.001210, 0.003327, 0.005627, 0.018450, 0.053066, 0.139813)
    
    #temporal
    dt     <- 0.1             #time step (days)
    times  <- 0:180 #365      #days sequence
    nt     <- (max(times)-min(times))/dt + 1       #no. time points, iterations
    nw     <- ceiling((max(times)-min(times))/7)   #week length of model run
    nd     <- ceiling((max(times)-min(times)))+1   #days length of model run
    
    #demography
    ages   <- c("0 to 4","5 to 11","12 to 17","18 to 29","30 to 39","40 to 49","50 to 59","60 to 69", "70+")
    na     <- 9               #number of age groups
    nimd   <- 5               #number of SE groups
    urban  <- T               #area: urban (T), rural (F)
    ageons <- c(0.0466, 0.0873, 0.0693, 0.14997, 0.1337, 0.1258, 0.1351, 0.1058, 0.1358); ageons=ageons/sum(ageons) #2020 mid
    
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
    rrep <- rep(0.5,9) 

})


