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


## Contact matrix (square, ng x ng where ng = na*nimd)
# Built from Reconnect base_matrix.csv by codes/prepare_model_inputs.R.
cm45<-(as.matrix(read.csv(paste0(input_dir,"/Mas50.csv"),header=F))) # removes name of columns
cm45dim1 = dim(cm45)[1]


## Parameters
if(pset$Vaccination==0){
   if(pset$Disease=="COVID-19")    source(paste0(source_dir,"/parsC_.r"))
   if(pset$Disease=="Influenza")   source(paste0(source_dir,"/parsF_.r"))
   if(pset$Disease=="RSV-illness") source(paste0(source_dir,"/parsR_.r"))
}else{
   if(pset$Disease=="COVID-19")    source(paste0(source_dir,"/parsCv_.r"))
   if(pset$Disease=="Influenza")   source(paste0(source_dir,"/parsFv_.r"))
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


## Demography
# Built from ONS LSOA SYA + IoD 2025 by codes/prepare_model_inputs.R.
demog2021 <- read.csv(paste0(input_dir,"/demographics_10age.csv"),header=T)
# number of age groups
na   = pars$na
# number of SES
nimd = pars$nimd
# number of groups
ng   = na*nimd
# proportion by age (summed across IMD)
pa<-vector(); for (i in 1:na){pa[i]=sum(demog2021$Population[which(demog2021$Age==pars$ages[i])])/sum(demog2021$Population)}
# Assign back into pars so anything reading pars$ageons gets the CSV-derived value
pars$ageons <- pa


## Initial state: S, E1:2, I1:2, U1:2, R, D
oNg  <- vector();   # 1/Population
Sg0  <- vector();   # Susceptible - Initial population, unless there's acquired immunity
E1g0 <- rep(0,ng);  # Exposed     - seed of epidemic
E2g0 <- rep(0,ng);  # Exposed
U1g0 <- rep(0,ng);  # Pre-clinical cases
U2g0 <- rep(0,ng);  # Pre-clinical cases
I1g0 <- rep(0,ng);  # Sub-clinical cases
I2g0 <- rep(0,ng);  # clinical cases
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

  
# pars: age 30 to 39, imd=1, 1/100,000 latent infections
E1g0 = (1/oNg)*pars$pE1g0
Sg0  = Sg0 - E1g0


## R0 and average contacts
source(paste0(source_dir,"/R0_.r"))      #outputs av contact rate
betanew = R0(pars,as.numeric(pars$R0),0) #default 2.5
print(paste0("Assuming R0 = ", pars$R0 ,"... beta is ", round(betanew,4)) )


## Parameters
# Note: pars already contains h, mH, rH (and for vacc: vcov, VE_inf/sym/hosp/sev, rV, rW, rW_nat) - they flow through via within().
parscpp45 = within(parscpp45 <- pars, {
                 cm=as.vector(cm45); cmdim1=cm45dim1; beta=betanew;
                 Sg0=Sg0; E1g0=E1g0; E2g0=E2g0; I1g0=I1g0; I2g0=I2g0; U1g0=U1g0; U2g0=U2g0;
                 Rg0=Rg0; Dg0=Dg0; oNg=oNg })
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


## fig 4 - all diseases (currently DISABLED: only RSV runs in main_repo.R).
## To re-enable: uncomment the Influenza + COVID-19 blocks in main_repo.R AND
## this section below. References p1C/p2C/p3C/p1F/p2F/p3F that are only set
## when those diseases are also run.
# if (pset$platform=="repo" & pars$Disease=="RSV-illness"){
#   filename=paste0("All_diseases_",area,"_SEIRD_epidemics_",daily,pset$Namevacc,TODAY)
#   pdf(file=paste0(output_dir,"/",filename,".pdf"))
#      gridExtra::grid.arrange(p1C,p2C,p3C,p1F,p2F,p3F,p1R,p2R,p3R, nrow=3, ncol=3)
#   dev.off()
#
#   gridExtra::grid.arrange(p1C,p2C,p3C,p1F,p2F,p3F,p1R,p2R,p3R, nrow=3, ncol=3)
#
#   ggsave(paste0(output_dir,"/",filename,".png"),
#          gridExtra::grid.arrange(p1C,p2C,p3C,p1F,p2F,p3F,p1R,p2R,p3R, nrow=3, ncol=3),
#          device = "png", width = 8000, height = 3931, units = "px", dpi = 600)
# }

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


