pset <- list()
pset <- within(pset, {
    TODAY <- format(Sys.Date(), "%d-%m-%Y")
    TIME  <- format(Sys.time(),'%H.%M.%S_%d-%m-%Y')

    Disease        <- "Influenza" #"RSV-illness" #"COVID-19"#
    Vaccination    <- 0 #1 #
    DailyIncidence <- 1 #0
    if(DailyIncidence==1){Incidence="Daily"}   else {Incidence="Weekly"}
    if(Vaccination==1)   {Namevacc="vaccine_"} else {Namevacc=""}
    
    DIAGNOSTIC   <- 1 #0
    ncomparisons <- 1

    FIGURES      <- 1 #0 #1
    SUMMARY      <- 1 #0 #1
    
	  COMPILE      <- 1
	platform       <- "repo" # "pc"
	
	if(platform=="repo") TODAY=""
})