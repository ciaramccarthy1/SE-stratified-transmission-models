pset <- list()
pset <- within(pset, {
    TODAY <- format(Sys.Date(), "%d-%m-%Y")
    TIME  <- format(Sys.time(),'%H.%M.%S_%d-%m-%Y')

    Disease        <- "RSV-illness"
    Vaccination    <- 0 #1 #
    DailyIncidence <- 1 #0
    if(DailyIncidence==1){Incidence="Daily"}   else {Incidence="Weekly"}
    if(Vaccination==1)   {Namevacc="vaccine_"} else {Namevacc=""}
    
    DIAGNOSTIC   <- 1 #0
    ncomparisons <- 1

    FIGURES      <- 1 #0 #1
    SUMMARY      <- 1 #0 #1
    
	  COMPILE      <- 1
	# TEMP (testing): use David's rsvie contact matrix by default (Mas45_david.csv).
	# Set FALSE to revert to the Reconnect matrix (Mas45.csv).
	DavidContacts  <- TRUE
	# Reproduce David's results: use his (stationary) populationAgeGroup instead of
	# real ONS demographics. TRUE (with DavidContacts=TRUE) closely reproduces his
	# age-stratified output; FALSE (default) uses real ONS demographics_9age.csv.
	DavidDemog     <- FALSE
	# Continuous demographic ageing (constant per-band rates: births into the
	# youngest band, ageing up, deaths out of the top). FALSE = frozen population (fine for
	# single-season reproduction). Set TRUE for multi-year runs where cohorts must age.
	Ageing         <- TRUE
	platform       <- "repo" # "pc"
	
	if(platform=="repo") TODAY=""
})