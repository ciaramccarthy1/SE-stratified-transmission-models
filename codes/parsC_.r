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
    #j1n - sh estimates
    age = c(mean(0:4),mean(5:11),mean(12:17),mean(18:29),mean(30:39),mean(40:49),mean(50:59),mean(60:69),mean(70:90)) #85))
    yA_0=0.438 #Med: 0.438, CI: [0.305,0.663]
    yr_0=0.012 #Med: 0.012, CI: [0.007,0.019]
    y[3:9]=yA_0*exp((age[3:9]-age[9])*yr_0)                         
    #0.2900 0.2700 0.1996 0.2223 0.2537 0.2861 0.3225 0.3637 0.4380 #derived sh est
    #Clinical fraction - by age and IMD group
    y45  <- rep(y,5)
    #TODO: estimate by IMD for ages 1:9 
    #TODO: Re-weight u y h again with final demography
    
    #Hospitalisation fraction
    #h    <- c(0, 0.0023, 0.0080, 0.0080, 0.0100, 0.0190, 0.0540, 0.1510, 0.4755) #derived age-adjusted Davies Nat Med 2020
    #        h[1] <- 0.0017   #age1, log lin regression
    #Mortality fraction (in hospital)
    #m    <- c(0.0000, 0.004, 0.005, 0.006, 0.018, 0.048, 0.094, 0.193, 0.411) #derived sh est
    #        m[1] <- 0.0021   #age1, log lin regression
    #m_0  <- c(0.0000, 0.015, 0.009, 0.005, 0.016, 0.040, 0.084, 0.182, 0.409) #derived sh est

    #Mortality fraction if clinically infected (derived from Verity 2020 IFR and y_sh_est mort-critically-Infected_18mar25.r)
    m <- c(0.000056, 0.000116, 0.000348, 0.001210, 0.003327, 0.005627, 0.018450, 0.053066, 0.139813)
    
    #temporal
    dt     <- 0.1             #0.01 #time step (days)
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
    rEI    <- 1/3 #1/2        #latency = rEU, Davies 2020 Nat Med
    rI1I2  <- 1/2.1           #recovery, Davies 2020 Nat Med
    rI2R   <- 1/2.9           #recovery, Davies 2020 Nat Med
    rIR    <- 1/(1/rI1I2 + 1/rI2R)  #1/3
    rUR    <- 1/5 #rIR        #recovery, Davies 2020 Nat Med
    #rIH    <- 1/8.5          #hospitalisation, Davies Lancet PH
    #rHR    <- 1/12.00        #recovery rate in hospital - #derived sh est
    #rHD    <- 1/13.91        #death rate in hospital	   - #derived sh est
    #rC     <- 1/8.5          #rate of loss of positivity, Russell et al; Davies Lancet ID

    R0     <- 2.5 #1.8        #Assumed, close to Knock 2021 and Davies 2020
    f      <- 0.5             #relative transmission of U group, Davies 2020 Nat Med
    beta   <- 0.06            #transmission rate between two given individuals, LG 2023

    #Initial condition
    #imd=1, age 30 to 39", 1/100,000 latent infections
    pE1g0   <- rep(0,na*nimd)  #initialise proportion latently infected across age x SES groups
    pE1g0[5] = (1/10^5)        #1/100,000 latent infections in age group 5 in SES 1
    
    #rate of reporting by age #Assumed
    rrep <- rep(0.5,9) 
    
    #NB parameters
    #k      <- 1               #dispersion/shape par of NB likelihood - #derived sh fit

})


