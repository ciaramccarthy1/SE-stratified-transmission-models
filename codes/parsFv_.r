pars <- list()
pars <- within(pars, {
  
    Disease     <- "Influenza"
    Vaccination <- "Yes"
    Incidence   <- pset$Incidence

    #Clinical responses
    #Susceptibility - variant & age-adjusted Baguelin 2013, Fig 22, 36, S52-54 2007-08, H3N2 dominant over H1N1, B
    # TODO(10-age scaffold): 60-64 inherits old 60-69; 65-74 mean(old 60-69, 70+); 75+ inherits old 70+
    u   <- c(0.63750, 0.63750, 0.50625, 0.37500, 0.37500, 0.37500, 0.37500, 0.37500, 0.37500, 0.37500)
    #Critically infected fraction
    y   <- rep(1,10)   #Treat clinical and sub-clin infections similarly, with some infections potentially causing death
    #y   <- rep(0.55,10) #most common value in flu studies in multipliers
    #Clinical fraction - by age and IMD group
    y45 <- rep(y,5)

    #Hospitalisation fraction
    #h    <-
    #Mortality fraction (in hospital)
    #m    <-

    age = c(mean(0:4),mean(5:14),mean(15:19),mean(20:29),mean(30:39),mean(40:49),mean(50:59),mean(60:64),mean(65:74),mean(75:90))

    #Mortality fraction if clinically infected
    #-IFR from LG Global paper - age-adjusted from 4 age groups
    # TODO(10-age scaffold): 60-64 inherits old 60-69; 65-74 mean(old 60-69, 70+); 75+ inherits old 70+
    IFR <- c(0.000053, 0.000008, 0.000008, 0.000116, 0.000137, 0.000137, 0.000137, 0.002690, (0.002690+0.005243)/2, 0.005243)
    m <- IFR/y

    #Hospitalisation pathway (H compartment): I2 -> H at rate h*rI2R; H -> D at rate mH*rH; H -> R at rate (1-mH)*rH
    # TODO(H scaffold): placeholder h = pmin(10*m, 0.8); mH = m/h so total mortality given clinical preserved.
    h  <- pmin(10*m, 0.8)
    mH <- m / h
    rH <- 1/6                     #1/hospital stay length (days^-1); TODO use flu-specific value
    #m <- c(0.000053, 0.000008, 0.000008, 0.000116, 0.000137, 0.000137, 0.000137, 0.002690, 0.005243) #y=1
    #m <- c(0.000096, 0.000015, 0.000015, 0.000211, 0.000249, 0.000249, 0.000249, 0.004891, 0.009533) #y=psym=0.55

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
    urban  <- T               #area: urban (T), rural (F)
    # TODO(10-age scaffold): old 60-69 split 50/50 -> 60-64 + half of 65-74; old 70+ split 5:15 -> half of 65-74 + 75+. Replace with ONS source.
    ageons <- c(0.0466, 0.0873, 0.0693, 0.14997, 0.1337, 0.1258, 0.1351, 0.1058/2, 0.1058/2 + 0.1358/4, 0.1358*3/4); ageons=ageons/sum(ageons) #2020 mid
    
    #natural history            
    rEI    <- 1/0.8           #latency,  Baguelin 2013
    rI1I2  <- 2/1.8           #recovery, Baguelin 2013
    rI2R   <- 2/1.8           #recovery, Baguelin 2013
    rIR    <- 1/(1/rI1I2 + 1/rI2R)
    rUR    <- 1/1.8 #rIR      #recovery, Baguelin 2013
    #rIH    <-                #hospitalisation
    #rHR    <-                #recovery rate in hospital 
    #rHD    <-                #death rate in hospital

    R0     <- 1.95            #variant adjusted from Baguelin 2013 Fig 22, 36, S52-54 2007-08
    f      <- 0 #0.5          #relative transmission of U group
    beta   <- 0.16            #variant adjusted from Baguelin 2013 Fig 22, 36, S52-54 2007-08
    
    #Initial condition
    pE1g0   <- rep(0,na*nimd)  #initialise proportion latently infected across age x SES groups
    pE1g0[5] = (1/10^5)        #1/100,000 latent infections in age group 5 in SES 1 (age 30 to 39, imd=1)
    
    #rate of reporting - variant & age-adjusted Baguelin 2013, Fig 22, 36, S52-54 2007-08, H3N2 dominant over H1N1, B
    # TODO(10-age scaffold): 60-64 inherits old 60-69; 65-74 mean of old 60-69 & 70+; 75+ inherits old 70+
    rrep <- c(0.004000, 0.004000, 0.014500, 0.025000, 0.025000, 0.025000, 0.025000, 0.018125, (0.018125+0.011250)/2, 0.011250)

    #vaccines: leaky V compartment with breakthrough shadow chain (V/Ev/Iv/Uv/Hv/Rv/Dv)
    # TODO(V scaffold): all VEs uniform across (age, IMD) at placeholder values; review for flu programme.
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


