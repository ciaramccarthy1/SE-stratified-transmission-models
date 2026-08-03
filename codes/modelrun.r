#pipeline SEIRD Rcpp
require(bench)
require(magrittr)
require(ggplot2)
require(ggtext)
require(gridExtra)
require(Rcpp)
require(tidyverse)

#TODO: vaccination refinement, Risk groups

### folders
input_dir0 <- here()
input_dir  <- here("data")
if (pset$platform=="repo"){
source_dir <- here("codes")
output_dir <- here("output")
TODAY      <- pset$TODAY
} else {
source_dir <- here()
output_dir <- here()
TODAY      <- format(Sys.Date(), "%d-%m-%Y")
}

#Disease choice
#source(paste0(source_dir,"/setup.r"))

## Contact matrix is read below, after pars.

## Parameters 
if(pset$Vaccination==0){
   if(pset$Disease=="RSV-illness") source(paste0(source_dir,"/parsR_.r"))
}else{
   if(pset$Disease=="RSV-illness") source(paste0(source_dir,"/parsRv_.r"))
}
# Scenario hook: codes/scenarios.R sets `scenario_overrides` as a named list
# (e.g. list(vcov = ...)) to override per-run pars entries before the cpp call.
if (exists("scenario_overrides", inherits = TRUE)) {
  for (.nm in names(scenario_overrides)) pars[[.nm]] <- scenario_overrides[[.nm]]
}
print(paste0("Disease:     ", pars$Disease))
print(paste0("Vaccination: ", pars$Vaccination))
print(paste0("Incidence  : ", pars$Incidence))
##pars$rV=0 #testing vs non-vaccine code

# number of age groups
na   = pars$na
# number of SES
nimd = pars$nimd
# number of groups
ng   = na*nimd

## Contact matrix (ng x ng): David's transmission matrix (cnt_matrix_p, built by
## codes/build_david_contacts.R -> Mas45_david.csv) when pset$DavidContacts=TRUE,
## else the Reconnect matrix (Mas45.csv).
# Ensemble hook: pars$cm_override (a ng x ng matrix) bypasses the file read
# (used by codes/uncertainty_ensemble.R to pass a per-posterior-draw contact matrix).
if (!is.null(pars$cm_override)) {
  cm45 <- as.matrix(pars$cm_override); dimnames(cm45) <- NULL
} else {
  cm_file <- if (isTRUE(pset$DavidContacts)) "/Mas45_david.csv" else "/Mas45.csv"
  cm45 <- as.matrix(read.csv(paste0(input_dir, cm_file),header=F))
}
cm45dim1 = dim(cm45)[1]

## Demography
## Demographics: real ONS by default, or David's stationary populationAgeGroup when
## pset$DavidDemog=TRUE (to reproduce David's results; built by build_david_contacts.R).
demog_file <- if (isTRUE(pset$DavidDemog)) "/demographics_9age_david.csv" else "/demographics_9age.csv"
demog2021 <- read.csv(paste0(input_dir, demog_file),header=T)
# proportion by age (summed across IMD)
pa<-vector(); for (i in 1:na){pa[i]=sum(demog2021$Population[which(demog2021$Age==pars$ages[i])])/sum(demog2021$Population)}
# Assign back into pars so anything reading pars$ageons gets the CSV-derived value
pars$ageons <- pa

## Initial state: S, E, I, U, R, D  (single-stage, no Erlang)
oNg  <- vector();   # 1/Population
Sg0  <- vector();   # Susceptible - Initial population, unless there's acquired immunity
Eg0  <- rep(0,ng);  # Exposed     - seed of epidemic
Ug0  <- rep(0,ng);  # Asymptomatic cases
Ig0  <- rep(0,ng);  # Clinical cases
Rg0  <- rep(0,ng);  # Recovered
Dg0  <- rep(0,ng);  # Dead
# One row per (IMD, age) in demog2021 - look up by filter
for (is in 1:nimd) {
  for (ia in 1:na) {
    pop <- demog2021$Population[demog2021$IMD == is & demog2021$Age == pars$ages[ia]]
    Sg0[(is-1)*na + ia] <- pop
    oNg[(is-1)*na + ia] <- 1/pop
  }
}
## Population by age group (over SES), by SES (over age), overall
Na<-rep(0,na)
Ns<-rep(0,nimd)
for (ia in 1:na) { for (is in 1:nimd) {
    Na[ia] = Na[ia] + 1/oNg[(is-1)*na + ia] 
    Ns[is] = Ns[is] + 1/oNg[(is-1)*na + ia] }}
Npop = sum(1/oNg);

## Demographic ageing rates (continuous; constant per-band rates). OFF unless
## pset$Ageing is TRUE. When OFF, eta/births are zero and the ageing terms in the cpp vanish.
##   eta[a] = 1/(365*width_a) for a<na  -> flow from band a to a+1
##   eta[na]                            -> death rate from the top band, set so total deaths
##                                         = total births (keeps David's population stationary)
##   births[s] = daily births into IMD s's youngest band (= its 0-4 outflow)
if (isTRUE(pset$Ageing)) {
  ag_lo <- c(0,5,15,25,35,45,55,65,75); ag_hi <- c(5,15,25,35,45,55,65,75,90)  # 9-band edges
  stopifnot(na == length(ag_lo))
  eta      <- 1/(365*(ag_hi-ag_lo))                    # ageing-out rate per band
  flow     <- eta[1]*Na[1]                             # constant demographic flow = daily births
  eta[na]  <- flow / Na[na]                            # top-band death rate balances births
  births   <- sapply(1:nimd, function(s) eta[1] * (1/oNg[(s-1)*na + 1]))  # per-IMD births into 0-4
} else {
  eta    <- rep(0, na)
  births <- rep(0, nimd)
}


## Season-start initial conditions
if (pars$Disease == "RSV-illness") {
  ## From David's rsvie post-burn-in state (all exposure groups). Collapse his 25
  ## age bands -> our 9 (weighted by his band pop Ntot over age overlap), then
  ## apply the fractions to our demographics_9age.csv populations. Our bands are
  ## aligned to David's so each of his 25 bands nests cleanly in one of ours.
  ## Pairs with the age-susceptibility u = sigma_a in parsR_.r (same exposure mix).
  suppressPackageStartupMessages(require(data.table))
  # Ensemble hook: pars$ic_states (25-band data.frame with age_group, Ntot,
  # frac_R/frac_E/frac_A/frac_I) bypasses the file read (per-posterior-draw IC).
  david <- if (!is.null(pars$ic_states)) data.table::as.data.table(pars$ic_states) else
           data.table::fread(paste0(input_dir, "/init_conditions_allgroups.csv"))
  data.table::setorder(david, age_group); stopifnot(nrow(david) == 25)
  # David's 25 band edges [lo,hi) in years (uk_data$ageGroupBoundary; last hi=90)
  d_lo <- c((0:11)/12, 1,2,3,4, 5,10, 15,25,35,45,55,65,75)
  d_hi <- c((1:12)/12, 2,3,4,5, 10,15, 25,35,45,55,65,75,90)
  # Our 9 David-aligned bands (pars$ages order)
  o_lo <- c(0,5,15,25,35,45,55,65,75)
  o_hi <- c(5,15,25,35,45,55,65,75,90)
  wmean <- function(x, w) if (sum(w) > 0) sum(x*w)/sum(w) else 0
  fcols <- c("frac_R","frac_E","frac_A","frac_I")
  ic9 <- sapply(fcols, function(cc) vapply(1:na, function(b) {
    ov <- pmax(0, pmin(o_hi[b], d_hi) - pmax(o_lo[b], d_lo))   # overlap years
    wmean(david[[cc]], david$Ntot * ov/(d_hi - d_lo)) }, numeric(1)))
  stopifnot(all(ic9 >= 0), all(rowSums(ic9) <= 1))
  for (is in 1:nimd) for (ia in 1:na) {
    g   <- (is-1)*na + ia
    pop <- 1/oNg[g]
    Rg0[g] <- pop * ic9[ia, "frac_R"]
    Eg0[g] <- pop * ic9[ia, "frac_E"]
    Ug0[g] <- pop * ic9[ia, "frac_A"]
    Ig0[g] <- pop * ic9[ia, "frac_I"]
    Sg0[g] <- pop - Rg0[g] - Eg0[g] - Ug0[g] - Ig0[g]   # = pop*frac_S, exact conservation
  }
} else {
  # COVID-19 / Influenza: single-cell latent seed (age 30-39, imd=1)
  Eg0 = (1/oNg)*pars$pE1g0
  Sg0 = Sg0 - Eg0
}


## R0 and average contacts
source(paste0(source_dir,"/R0_.r"))      #outputs av contact rate
betanew = R0(pars,as.numeric(pars$R0),0) #default 2.5
# Fitting hook: set pars$beta_override (e.g. via scenario_overrides) to use a
# directly-chosen beta instead of the R0-derived value (for eyeball fitting).
if (!is.null(pars$beta_override)) betanew <- as.numeric(pars$beta_override)
print(paste0("Assuming R0 = ", pars$R0 ,"... beta is ", round(betanew,4)) )

## Parameters
# Note: pars already contains h, mH, rH (and for vacc: vcov, VE_inf/sym/hosp/sev, rV, rW, rW_nat) - they flow through via within().
parscpp45 = within(parscpp45 <- pars, {
                 cm=as.vector(cm45); cmdim1=cm45dim1; beta=betanew;
                 Sg0=Sg0; Eg0=Eg0; Ig0=Ig0; Ug0=Ug0;
                 Rg0=Rg0; Dg0=Dg0; oNg=oNg;
                 # recovery rates age-varying (length na); rep_len tolerates scalars
                 rIR=rep_len(as.numeric(rIR), na); rUR=rep_len(as.numeric(rUR), na);
                 eta=eta; births=births })   # demographic ageing (zeros unless pset$Ageing)
#  for output
parsum = parscpp45
#  remove what's not needed for Rcpp:
parscpp45 <- parscpp45 %>% magrittr::inset(c('age', 'ages', 'ageons', 'm'), NULL)  #parscpp45[['age']] <- NULL; etc

## Model output (for the proposed parameters)
if (pset$COMPILE==1) {
#pset$Vaccination==0 or pset$Vaccination==1
                if(pset$DailyIncidence==0) sourceCpp(file = paste0(source_dir,"/","SEIRDas_.cpp"))
                if(pset$DailyIncidence==1) sourceCpp(file = paste0(source_dir,"/","SEIRDasday_.cpp"))
                mas <- model(parscpp45)
## scale of plots without vaccination
                                           IUw_novacc    = max(10^5*c(mas$byw$Iw, mas$byw$Uw)/Npop)
                                           # mas$byw$Iw_s: matrix (nd|nw) x nimd; mas$byaw$Iw_a: matrix x na
                                           Iw_imd_novacc = max(10^5 * sweep(mas$byw$Iw_s,  2, Ns, "/"))
                                           Iw_age_novacc = max(10^5 * sweep(mas$byaw$Iw_a, 2, Na, "/"))
				
if (pset$Vaccination==1) {
                if(pset$DailyIncidence==0) sourceCpp(file = paste0(source_dir,"/","SEIRDasvacc_.cpp"))
                if(pset$DailyIncidence==1) sourceCpp(file = paste0(source_dir,"/","SEIRDasvaccday_.cpp"))
				mas <- model(parscpp45)
} }
Iwpeakval = max(mas$byw$Iw)*10^(-6)
Iwpeakloc = mas$byw$time[which(mas$byw$Iw==max(mas$byw$Iw))]
print(paste0("Peak:  ", round(Iwpeakval,3) ," million at ", Iwpeakloc, " days"))
cat("\n")


## Figures
if (pset$FIGURES==1){
  
ar=1 #aspect ratio

if(pars$Incidence=="Daily"){daily="daily_"} else {daily=""}
filename=paste0(parsum$Disease,"_SEIRD_epidemic_",daily,pset$Namevacc,TODAY)
pdf(file=paste0(output_dir,"/",filename,".pdf"))


## fig 1 overall
data <- data.frame(time=rep(mas$byw$time,2), IUw=10^5*c(mas$byw$Iw, mas$byw$Uw)/Npop,
                   State=rep(c("Clinic","Unasc"),each=length(mas$byw$time)))

p1 <- ggplot(data, aes(x=time)) + 
      geom_line(aes(y = IUw, group=State, color=State), lwd=0.8)  +
      theme(text=element_text(size=10),
            legend.key.size = unit(2, 'mm'),
            plot.title = element_text(size = 12),
            axis.text.y = element_text(color=1),
            axis.text.x = element_text(color=1)) +
      labs(y = "Infectious inc./100k/week", x = "Day", color = "State") +
	  ylim(0,IUw_novacc) +
      ggtitle(parsum$Disease) #+ theme(aspect.ratio=ar)

print(p1)
if (pset$platform=="repo" & pars$Disease=="RSV-illness") p1R<-p1
if (pset$platform=="repo" & pars$Disease=="Influenza")   p1F<-p1
if (pset$platform=="repo" & pars$Disease=="COVID-19")    p1C<-p1


## fig 2 by ses
# Iw_s/Uw_s are (nd|nw) x nimd matrices; flatten column-wise so col j -> IMD j
IUw_s_mat <- mas$byw$Iw_s + mas$byw$Uw_s
data <- data.frame(
  time = rep(mas$byw$time, nimd),
  IUw  = 10^5 * as.vector(sweep(IUw_s_mat,    2, Ns, "/")),
  Iw   = 10^5 * as.vector(sweep(mas$byw$Iw_s, 2, Ns, "/")),
  IMD  = rep(1:nimd, each = length(mas$byw$time)))
p2 <- ggplot(data, aes(x=time)) + 
      geom_line(aes(y = Iw, group=IMD, color=IMD), lwd=0.8)  +
      theme(text=element_text(size=10),
            legend.key.size = unit(2, 'mm'),
            plot.title = element_text(size = 12),
            axis.text.y = element_text(color=1),
            axis.text.x = element_text(color=1)) +
      labs(y = "Clinical infs. /100k/week", x = "Day", color = "IMD") +
	  ylim(0,Iw_imd_novacc) +
      ggtitle(parsum$Disease) #+ theme(aspect.ratio=ar)

print(p2)
if (pset$platform=="repo" & pars$Disease=="RSV-illness") p2R<-p2
if (pset$platform=="repo" & pars$Disease=="Influenza")   p2F<-p2
if (pset$platform=="repo" & pars$Disease=="COVID-19")    p2C<-p2


## fig 3 by age
# Iw_a/Uw_a are (nd|nw) x na matrices; flatten column-wise so col j -> AGE j
IUw_a_mat <- mas$byaw$Iw_a + mas$byaw$Uw_a
data <- data.frame(
  time = rep(mas$byw$time, na),
  IUw  = 10^5 * as.vector(sweep(IUw_a_mat,     2, Na, "/")),
  Iw   = 10^5 * as.vector(sweep(mas$byaw$Iw_a, 2, Na, "/")),
  AGE  = rep(1:na, each = length(mas$byw$time)))
p3 <- ggplot(data, aes(x=time)) + 
      geom_line(aes(y = Iw, group=AGE, color=AGE), lwd=0.8)  +
      theme(text=element_text(size=10),
        legend.key.size = unit(2, 'mm'),
        plot.title = element_text(size = 12),
        axis.text.y = element_text(color=1),
        axis.text.x = element_text(color=1)) +
      labs(y = "Clinical infs. /100k/week", x = "Day", color = "Age") +
	  ylim(0,Iw_age_novacc) +
      ggtitle(parsum$Disease)  #+ theme(aspect.ratio=ar)

print(p3)
dev.off()
if (pset$platform=="repo" & pars$Disease=="RSV-illness") p3R<-p3
if (pset$platform=="repo" & pars$Disease=="Influenza")   p3F<-p3
if (pset$platform=="repo" & pars$Disease=="COVID-19")    p3C<-p3

}##FIGURES


## Performance

if (pset$DIAGNOSTIC==1){
   filename=paste0(parsum$Disease,"_SEIRD_performance_",daily,pset$Namevacc,TODAY)
   if(pset$ncomparisons==1){  #1 #2
     niter= 100 #3000 #1000
     test <- bench::mark(model(parscpp45), min_iterations = niter)
     #test <- bench::mark(mas, min_iterations = niter)
     sink(file = paste0(output_dir,"/",filename,".txt"),append=FALSE,split=FALSE)
        print(paste0("Iterations = ", niter))
        print(test[1:11]) #cut last 2 columns
     sink()
     pdf(file=paste0(output_dir,"/",filename,".pdf"))
        print(plot(test))
        ## plot time vs memory allocation
        print(test %>% unnest(c(time, gc)) %>%
                 filter(gc == "none") %>%
                 mutate(expression = as.character(expression)) %>%
                 ggplot(aes(x = mem_alloc, y = time, color = expression)) + geom_point())
     dev.off()
     
     } else { #For comparisons
     niter= 3000 #1000
     test <- bench::mark(model(parscpp45), model2(parscpp45), min_iterations = 3000) #1000)
     sink(file = paste0(output_dir,"/",filename,".txt"),append=FALSE,split=FALSE)
        print(paste0("Iterations = ", niter))
        summary(test, relative = TRUE)
        summary(test, relative = FALSE)
     sink()
     pdf(file=paste0(output_dir,"/",filename,".pdf"))
        plot(test)
        ## plot time vs memory allocation
        test %>% unnest(c(time, gc)) %>%
                 filter(gc == "none") %>%
                 mutate(expression = as.character(expression)) %>%
                ggplot(aes(x = mem_alloc, y = time, color = expression)) + geom_point()
     dev.off()
     
  }#comparisons
  
  
}#Diagnostic



## Text summary

if (pset$SUMMARY==1){

filename=paste0(parsum$Disease,"_SEIRD_parameters_",daily,pset$Namevacc,TODAY)
sink(file = paste0(output_dir,"/",filename,".txt"),append=FALSE,split=FALSE)

cat("\n")

cat("\n Iw peak \n")
print(paste0("Peak:  ", round(Iwpeakval,3) ," million at ", Iwpeakloc, " days"))

cat("\n Study \n")
print(paste0("Disease:    ", parsum$Disease))
print(paste0("Population: ", Npop))
print(paste0("Age groups: ", parsum$na))
print(paste0("SE  groups: ", parsum$nimd))
print(paste0("All groups: ", parsum$na*parsum$nimd))
print(paste0("Age distribution: ")); print(pa*Npop)
print(paste0("Age proportions:  ")); print(round(pa,4))
print(paste0("Ages:             ")); print(parsum$ages)
print(paste0("Age (median):     ")); print(parsum$age)
print(paste0("Reporting rate:   ")); print(parsum$rrep)

cat("\n Natural history \n")
print(paste0("Assuming R0 = ", parsum$R0))
print(paste0("  then beta = (1/day) ", round(betanew,4) ))
print(paste0("latent period (tEI, tEU),   days:    ", 1/parsum$rEI))
print(paste0("infectious period (rIR),    days:    ", 1/parsum$rIR))
print(paste0("infectious period (rUR),    days:    ", 1/parsum$rUR))
print(paste0("infectious preclin (rI1I2), days:    ", 1/parsum$rI1I2))
print(paste0("infectious clinical (rI2R), days:    ", 1/parsum$rI2R))
print(paste0("relative subclinical infectiousness: ", parsum$f))
print(paste0("susceptibility by age     : ")); print(parsum$u)
print(paste0("clinical  fraction by age : ")); print(parsum$y)
print(paste0("mortality fraction (overall, given clinical) by age : ")); print(parsum$m)
print(paste0("hospitalisation fraction (given clinical) by age    : ")); print(parsum$h)
print(paste0("mortality fraction (given hospitalised)    by age   : ")); print(parsum$mH)
print(paste0("hospital stay length, days                          : ", round(1/parsum$rH,2)))

cat("\n Initial condition \n");
print(paste0("Initial latent proportion pE1g0: ")); print(as.numeric(parsum$pE1g0))
print(paste0("Initial latent infections  E1g0: ")); print(as.numeric(parsum$pE1g0*(1/oNg)))
print(paste0("Initial susceptible         Sg0: ")); print(Sg0)

cat("\n Temporal \n")
print(paste0("Time range:       ", range(pars$times)))
print(paste0("Number of weeks:  ", pars$nw))
print(paste0("Number of points: ", pars$nt))
print(paste0("dt:               ", pars$dt))

cat("\n Contacts \n")
print(paste0("Contact data: Polymod 2005"))
print(paste0("Average contact rate of cm45: ", round(cav,3), "/day"))
print(paste0("Contact matrix: ", parsum$cmdim1, " x ", parsum$cmdim1))
print(paste0("Contact matrix: ")); #cm

if (pset$Vaccination==1){
cat("\n Vaccination \n")
print(paste0("Coverage (per age x IMD):     ")); print(parsum$vcov)
print(paste0("VE against infection:         ")); print(parsum$VE_inf)
print(paste0("VE against symptoms:          ")); print(parsum$VE_sym)
print(paste0("VE against hospitalisation:   ")); print(parsum$VE_hosp)
print(paste0("VE against mortality (in H):  ")); print(parsum$VE_mort)
print(paste0("Vaccination rate (per day):   ", round(parsum$rV,5)))
print(paste0("Vaccine waning rate (per day):", round(parsum$rW,5)))
print(paste0("Natural waning rate (per day):", round(parsum$rW_nat,5))) }


cat("\n")
sink()

cat("\n")

}##SUMMARY


