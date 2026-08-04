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

	##########################################################################
	# MASTER SWITCH ---------------------------------------------------------
	#   MatchDavid = TRUE  -> reproduce David/rsvie EXACTLY (validation benchmark)
	#   MatchDavid = FALSE -> our extended model (ONS demog, Reconnect contacts,
	#                         realistic age-specific deaths, hospitalisation &
	#                         mortality in the dynamics, recalibrated beta)
	# Flip this one flag to switch every David-vs-ours fork together. Each fork
	# below defaults to MatchDavid; set a literal TRUE/FALSE to DECOUPLE one fork
	# for testing (e.g. ONS demog with David contacts).
	#
	# Forks driven by MatchDavid (David <-> ours):
	#   DavidContacts : David transmission matrix   <-> Reconnect contacts
	#   DavidDemog    : David rectangular population <-> ONS demographics
	#   [TODO death rates] top-band-only <-> ONS age-specific   (modelrun.r)
	#   [TODO h,m,mH]      0 / post-processed <-> in the dynamics (parsR_.r/parsRv_.r)
	#   [TODO beta]        David's fitted qp <-> recalibrated on Reconnect
	##########################################################################
	MatchDavid     <- FALSE

	DavidContacts  <- MatchDavid   # FALSE -> Reconnect matrix (Mas45.csv)
	DavidDemog     <- MatchDavid   # FALSE -> real ONS demographics_9age.csv
	# Continuous demographic ageing: TRUE in BOTH versions (David ages too - his
	# eta = 1/(365*width) + births/deaths). Only the DEATH structure differs
	# between versions (top-band-only vs ONS age-specific), wired via MatchDavid
	# in modelrun.r once the age-specific death rates are added.
	Ageing         <- TRUE
	platform       <- "repo" # "pc"
	
	if(platform=="repo") TODAY=""
})