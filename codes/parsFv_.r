pars <- list()
pars <- within(pars, {
  
    Disease     <- "Influenza"
    Vaccination <- "Yes"
    Incidence   <- pset$Incidence
    
    #TODO: check y and year chosen from Baguelin
    #Clinical responses
    #Susceptibility - variant & age-adjusted Baguelin 2013, Fig 22, 36, S52-54 2007-08, H3N2 dominant over H1N1, B
    u   <- c(0.63750, 0.63750, 0.50625, 0.37500, 0.37500, 0.37500, 0.37500, 0.37500, 0.37500)
    #Critically infected fraction
    y   <- rep(1,9)    #Treat clinical ans sub-clin infections similarly, with some potentially causing death
    #y   <- rep(0.55,9) #most common value in flu studies in IFRdone, CFRsdone, multipliers
    #Clinical fraction - by age and IMD group
    y45 <- rep(y,5)

    #Hospitalisation fraction
    #h    <- 
    #Mortality fraction (in hospital)
    #m    <- 
    
    age = c(mean(0:4),mean(5:11),mean(12:17),mean(18:29),mean(30:39),mean(40:49),mean(50:59),mean(60:69),mean(70:90)) #85))

    #Mortality fraction if clinically infected 
    #-IFR from LG Global paper - age-adjusted from 4 age groups
    IFR <- c(0.000053, 0.000008, 0.000008, 0.000116, 0.000137, 0.000137, 0.000137, 0.002690, 0.005243)
    m <- IFR/y
    #m <- c(0.000053, 0.000008, 0.000008, 0.000116, 0.000137, 0.000137, 0.000137, 0.002690, 0.005243) #y=1
    #m <- c(0.000096, 0.000015, 0.000015, 0.000211, 0.000249, 0.000249, 0.000249, 0.004891, 0.009533) #y=psym=0.55

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
    rEI    <- 1/0.8           #latency,  Baguelin 2013
    rI1I2  <- 2/1.8           #recovery, Baguelin 2013
    rI2R   <- 2/1.8           #recovery, Baguelin 2013
    rIR    <- 1/(1/rI1I2 + 1/rI2R)
    rUR    <- 1/1.8 #rIR      #recovery, Baguelin 2013
    #rIH    <-                #hospitalisation
    #rHR    <-                #recovery rate in hospital 
    #rHD    <-                #death rate in hospital
    #rC     <-                #rate of loss of positivity

    R0     <- 1.95            #variant adjusted from Baguelin 2013 Fig 22, 36, S52-54 2007-08
    f      <- 0 #0.5          #relative transmission of U group
    beta   <- 0.16            #variant adjusted from Baguelin 2013 Fig 22, 36, S52-54 2007-08
    
    #Initial condition
    #imd=1, age 30 to 39", 1/100,000 latent infections
    pE1g0   <- rep(0,na*nimd)  #initialise proportion latently infected across age x SES groups
    pE1g0[5] = (1/10^5)        #1/100,000 latent infections in age group 5 in SES 1
    
    #rate of reporting - variant & age-adjusted Baguelin 2013, Fig 22, 36, S52-54 2007-08, H3N2 dominant over H1N1, B
    rrep <- c(0.004000, 0.004000, 0.014500, 0.025000, 0.025000, 0.025000, 0.025000, 0.018125, 0.011250)
    
    #vaccines
    ve   <- rep(0.5, na)        #efficacy 0.5
    veff <- rep(ve,  nimd)
    vc   <- rep(1, na)          #coverage 0.1
    vcov <- rep(vc,  nimd)
    vcln <- 0.15                #reduction in clinical fraction
    rV   <- 1/180               #rate of immunisation
  
    #NB parameters
    #k      <- 1               #dispersion/shape par of NB likelihood - #derived sh fit

})


