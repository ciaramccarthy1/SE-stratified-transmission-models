#pipeline SEIRD Rcpp

library(bench)
library(magrittr)
library(ggplot2)
library(ggtext)
library(Rcpp)
library(tidyverse)


#TODO: vaccination, Risk groups


### folders
input_dir0 <- getwd()
input_dir  <- paste0(getwd(),"/data")
if (pset$platform=="repo"){
source_dir <- paste0(getwd(),"/codes")
output_dir <- paste0(getwd(),"/output")
TODAY      <- pset$TODAY
#temporarily
#TODAY      <- paste0("_",format(Sys.Date(), "%d-%m-%Y"))
} else {
source_dir <- paste0(getwd()) #,"/codes")
output_dir <- paste0(getwd()) #,"/output")
TODAY      <- format(Sys.Date(), "%d-%m-%Y")
}


#Disease choice
#source(paste0(source_dir,"/setup.r"))


## Contact matrix 45x45
cm45<-(as.matrix(read.csv(paste0(input_dir,"/Mas45_urban.csv"),header=F))) # removes name of columns
cm45dim1 = dim(cm45)[1]


## Parameters
if(pset$Disease=="COVID-19")    source(paste0(source_dir,"/parsC_.r"))
if(pset$Disease=="Influenza")   source(paste0(source_dir,"/parsF_.r"))
if(pset$Disease=="RSV-illness") source(paste0(source_dir,"/parsR_.r"))
print(paste0("Disease: ", pars$Disease))


## Demography
demog2021 <- read.csv(paste0(input_dir,"/demographics2021.csv"),header=T)
# number of age groups
na   = pars$na
# number of SES
nimd = pars$nimd
# number of groups
ng   = na*nimd
# area - region
urb  = pars$urban #urban (T), rural (F)
if(pars$urb==T){ area = "Urban"} else { area ="Rural"}
print(paste0("Area: ", area))
# proportion by age
pa<-vector(); for (i in 1:na){pa[i]=sum(demog2021$Population[which(demog2021$Age==pars$ages[i])])/sum(demog2021$Population)}
#  check:
#  round(pa,4)          [1] 0.0573 0.0873 0.0693 0.1500 0.1337 0.1258 0.1351 0.1058 0.1358
#  round(pars$ageons,4) [1] 0.0471 0.0882 0.0700 0.1516 0.1351 0.1272 0.1366 0.1069 0.1373


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
for (is in 1:nimd) {
  Sg0[(is-1)*na + 1:na] =   demog2021$Population[1:na + na*(is-1) + na*nimd*(1-urb)] # whole urb population
  oNg[(is-1)*na + 1:na] = 1/demog2021$Population[1:na + na*(is-1) + na*nimd*(1-urb)] # 1/number in each age group
}
## Population by age group (over SES), by SES (over age), overall
Na<-rep(0,na)
Ns<-rep(0,nimd)
for (ia in 1:na) { for (is in 1:nimd) {
    Na[ia] = Na[ia] + 1/oNg[(is-1)*na + ia] 
    Ns[is] = Ns[is] + 1/oNg[(is-1)*na + ia] }}
Npop = sum(1/oNg);

  
# Checks:
#  sum(demog2021$Population)                                    #[1] 56550138
#  sum(demog2021$Population[which(demog2021$rural=="Rural")])   #[1]  9683314
#  sum(demog2021$Population[which(demog2021$rural=="Urban")])   #[1] 46866824
#  sum(1/oNa)                                                   #[1] 46866824
# pars: imd=1, age 30 to 39", 1/100,000 latent infections
E1g0 = (1/oNg)*pars$pE1g0
Sg0  = Sg0 - E1g0


## R0 and average contacts
source(paste0(source_dir,"/R0_.r"))      #outputs cav
betanew = R0(pars,as.numeric(pars$R0),0) #default 2.5
print(paste0("Assuming R0 = ", pars$R0 ,"... beta is ", round(betanew,4), "/day"))


## Parameters
parscpp45 = within(parscpp45 <- pars, {
                 cm=as.vector(cm45); cmdim1=cm45dim1; mI=pars$m; beta=betanew;
                 Sg0=Sg0; E1g0=E1g0; I1g0=I1g0; I2g0=I2g0; U1g0=U1g0; U2g0=U2g0; 
                 Rg0=Rg0; Dg0=Dg0; oNg=oNg })
#  for output
parsum = parscpp45;
#  remove what's not needed for Rcpp
parscpp45 <- parscpp45 %>% magrittr::inset(c('age', 'ages', 'ageons', 'm'), NULL)  #parscpp45[['age']] <- NULL; etc


## Model output (for the proposed parameters)
if (pset$COMPILE==1) sourceCpp(file = paste0(source_dir,"/","SEIRDas_.cpp"))
mas <- model(parscpp45)
Iwpeakval = max(mas$byw$Iw)*10^(-6)
Iwpeakloc = mas$byw$time[which(mas$byw$Iw==max(mas$byw$Iw))]
print(paste0("Peak:  ", round(Iwpeakval,3) ," (million) at ", Iwpeakloc, " days"))
cat("\n")


## Plotting overall
#  (for comparison with past model runs)
#plot(mas$byw$time,  (mas$byw$It)/(7/pars$dt),type="l",lty="solid",col=1,ylab="Infectious",xlab="Day")
#lines(mas$byw$time, (mas$byw$Ut)/(7/pars$dt),type="l",lty="solid",col=3,ylab="Infectious",xlab="Day")
#lines(mas$byw$time, (mas$byw$Iw)/7,type="l",lty="dashed",col=1)
#lines(mas$byw$time, (mas$byw$Uw)/7,type="l",lty="dashed",col=3)



## fig 1 overall
data <- data.frame(time=rep(mas$byw$time,2), IUw=10^5*c(mas$byw$Iw, mas$byw$Uw)/Npop,
                   State=rep(c("Clini","Unasc"),each=length(mas$byw$time)))

p1 <- ggplot(data, aes(x=time)) + 
      geom_line(aes(y = IUw, group=State, color=State), lwd=0.8)  +
      theme(text=element_text(size=10),
      legend.key.size = unit(3, 'mm'),
      plot.title = element_text(size = 13),
      axis.text.y = element_text(color=1),
      axis.text.x = element_text(color=1)) +
      labs(y = "Infectious incid. /100k/week", x = "Day", color = "State") +
      ggtitle(paste0(parsum$Disease,", ",area," all age & SE strata")) 

if (pset$platform=="repo" & pars$Disease=="RSV-illness") p1R<-p1
if (pset$platform=="repo" & pars$Disease=="Influenza")   p1F<-p1
if (pset$platform=="repo" & pars$Disease=="COVID-19")    p1C<-p1

#filename=paste0(parsum$Disease,"_",area,"_SEIRD_Iw_Uw_overall_",TODAY)
filename=paste0(parsum$Disease,"_",area,"_SEIRD_Infectious_incidence_",TODAY)
pdf(file=paste0(output_dir,"/",filename,".pdf"))
print(p1)
dev.off()


## fig 2 by ses
data <- data.frame(time=rep(mas$byw$time,5), 
      IUw=10^5*c(mas$byw$IUw_s1/Ns[1], mas$byw$IUw_s2/Ns[2], mas$byw$IUw_s3/Ns[3], mas$byw$IUw_s4/Ns[4], mas$byw$IUw_s5/Ns[5]),
      Iw =10^5*c(mas$byw$Iw_s1/Ns[1],  mas$byw$Iw_s2/Ns[2],  mas$byw$Iw_s3/Ns[3],  mas$byw$Iw_s4/Ns[4],  mas$byw$Iw_s5/Ns[4]),
                   IMD=rep(1:5,each=length(mas$byw$time)))
p2 <- ggplot(data, aes(x=time)) + 
      #geom_line(aes(y = IUw/10^6, group=IMD, color=IMD), lwd=0.8)  +
      geom_line(aes(y = Iw, group=IMD, color=IMD), lwd=0.8)  +
      theme(text=element_text(size=10),
      legend.key.size = unit(3, 'mm'),
      plot.title = element_text(size = 13),
      axis.text.y = element_text(color=1),
      axis.text.x = element_text(color=1)) +
      #labs(y = "Infectious clinical & subcl. incidence (10^6/week)", x = "Day", color = "IMD") +
      labs(y = "Clinical infs. /100k/week", x = "Day", color = "IMD") +
      ggtitle(paste0(parsum$Disease,", ",area," by SE strata"))

if (pset$platform=="repo" & pars$Disease=="RSV-illness") p2R<-p2
if (pset$platform=="repo" & pars$Disease=="Influenza")   p2F<-p2
if (pset$platform=="repo" & pars$Disease=="COVID-19")    p2C<-p2

#filename=paste0(parsum$Disease,"_",area,"_SEIRD_Iw+Uw_by_SES_",TODAY)
filename=paste0(parsum$Disease,"_",area,"_SEIRD_Clinical_indicence_by_SES_",TODAY)
pdf(file=paste0(output_dir,"/",filename,".pdf"))
print(p2)
dev.off()


## fig 3 by age
data <- data.frame(time=rep(mas$byw$time,9), 
       IUw=10^5*c(mas$byaw$IUw_a1/Na[1], mas$byaw$IUw_a2/Na[2], mas$byaw$IUw_a3/Na[3], mas$byaw$IUw_a4/Na[4], 
                  mas$byaw$IUw_a5/Na[5], mas$byaw$IUw_a6/Na[6], mas$byaw$IUw_a7/Na[7], mas$byaw$IUw_a8/Na[8], 
                  mas$byaw$IUw_a9/Na[9]),
        Iw=10^5*c(mas$byaw$Iw_a1/Na[1],  mas$byaw$Iw_a2/Na[2],  mas$byaw$Iw_a3/Na[3],  mas$byaw$Iw_a4/Na[4],  
                  mas$byaw$Iw_a5/Na[5],  mas$byaw$Iw_a6/Na[6],  mas$byaw$Iw_a7/Na[7],  mas$byaw$Iw_a8/Na[8],
                  mas$byaw$Iw_a9/Na[9]),
                   AGE=rep(1:9,each=length(mas$byw$time)))
p3 <- ggplot(data, aes(x=time)) + 
  #geom_line(aes(y = IUw/10^6, group=AGE, color=AGE), lwd=0.8)  +
  geom_line(aes(y = Iw, group=AGE, color=AGE), lwd=0.8)  +
  theme(text=element_text(size=10),
        legend.key.size = unit(3, 'mm'),
        plot.title = element_text(size = 13),
        axis.text.y = element_text(color=1),
        axis.text.x = element_text(color=1)) +
  #labs(y = "Infectious clinical & subcl. incidence (10^6/week)", x = "Day", color = "Age") +
  labs(y = "Clinical infs. /100k/week", x = "Day", color = "Age") +
  ggtitle(paste0(parsum$Disease,", ",area," by age group"))

if (pset$platform=="repo" & pars$Disease=="RSV-illness") p3R<-p3
if (pset$platform=="repo" & pars$Disease=="Influenza")   p3F<-p3
if (pset$platform=="repo" & pars$Disease=="COVID-19")    p3C<-p3

#filename=paste0(parsum$Disease,"_",area,"_SEIRD_Iw+Uw_by_Age_",TODAY)
filename=paste0(parsum$Disease,"_",area,"_SEIRD_Clinical_incidence_by_Age_",TODAY)
pdf(file=paste0(output_dir,"/",filename,".pdf"))
print(p3)
dev.off()


## fig 4 - all diseases

if (pset$platform=="repo" & pars$Disease=="RSV-illness"){
  filename=paste0("All_diseases_",area,"_SEIRD_all_indicators_",TODAY)
  #par(mfrow = c(3, 3))
  pdf(file=paste0(output_dir,"/",filename,".pdf"))
  #print(p1C)
  #print(p2C)
  #print(p3C)
  #print(p1F)
  #print(p2F)
  #print(p3F)
  #print(p1R)
  #print(p2R)
  #print(p3R)
  gridExtra::grid.arrange(p1C,p2C,p3C,p1F,p2F,p3F,p1R,p2R,p3R)
    dev.off()
}
  
  
## Performance
DIAGNOSTIC=pset$DIAGNOSTIC #1
ncomparisons=pset$ncomparisons #1 #2

if (DIAGNOSTIC==1){
   filename=paste0(parsum$Disease,"_SEIRD_performance_",TODAY)
   if(ncomparisons==1){
     niter= 100 #3000 #1000
     test <- bench::mark(model(parscpp45), min_iterations = niter)
     sink(file = paste0(output_dir,"/",filename,".txt"),append=FALSE,split=FALSE)
        print(paste0("Iterations = ", niter))
        print(test[1:11]) #cut last 2 columns
     sink()
     pdf(file=paste0(output_dir,"/",filename,".pdf"))
        print(plot(test))
        ## plot time vs memory allocation
        #  https://bench.r-lib.org/
        #  library(tidyr)
        print(test %>% unnest(c(time, gc)) %>%
                 filter(gc == "none") %>%
                 mutate(expression = as.character(expression)) %>%
                 ggplot(aes(x = mem_alloc, y = time, color = expression)) + geom_point())
     dev.off()
     
     } else { #For comparisons
     #sourceCpp(file = paste0(source_dir,"/","SEIURDasv_.cpp"))
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
     
     #if (require(ggplot2) && require(tidyr) && require(ggbeeswarm)) {
     #      # Beeswarm plot
     #        autoplot(test)
     #      # ridge (joyplot)
     #        autoplot(test, "ridge")
     #      # If you want to have the plots ordered by execution time you can do so by
     #      # ordering factor levels in the expressions.
     #          if (require(dplyr) && require(forcats)) {
     #                test %>%
     #                      mutate(expression = forcats::fct_reorder(as.character(expression), min, .desc = TRUE)) %>%
     #                      as_bench_mark() %>%
     #                      autoplot("violin")
     #          }}
     
     }#comparisons
  
  ##Profiling
  #sourceCpp(file = paste0(source_dir,"/","SEIRDa.cpp"))
  #profvis::profvis(SEIRDa(parscpp))
  #=> "Error in parse_rprof_lines(lines, expr_source) : 
  #    No parsing data available. Maybe your function was too fast?"
  
}#Diagnostic



## Text summary

sink(file = paste0(output_dir,"/",parsum$Disease,"_",area,"_SEIRD_parameters_",TODAY,".txt"),append=FALSE,split=FALSE)

cat("\n")
print(paste0("Disease:    ", parsum$Disease))
print(paste0("Area:       ", area))
print(paste0("Population: ", Npop))
print(paste0("Age groups: ", parsum$na))
print(paste0("SE  groups: ", parsum$nimd))
print(paste0("All groups: ", parsum$na*parsum$nimd))
print(paste0("Age distribution: ")); pa*Npop
print(paste0("Age proportions:  ")); round(pa,4)
print(paste0("Ages:             ")); parsum$ages
print(paste0("Age (median):     ")); parsum$age
print(paste0("Reporting rate:   ")); parsum$rrep

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
print(paste0("mortality fraction by age : ")); print(parsum$mI)

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
print(paste0("Average contact rate of cm45: ", round(cav,5)))
print(paste0("Contact matrix: ", parsum$cmdim1, " x ", parsum$cmdim1))
print(paste0("Contact matrix: ")); #cm

cat("\n Iw peak \n")
print(paste0("Peak:  ", round(Iwpeakval,3) ," (million) at ", Iwpeakloc, " days"))

cat("\n")
sink()

cat("\n")

#filename=paste0(output_dir,"/",parsum$Disease,"_",area,"_SEIRD_parameters_",TODAY,".txt")
#write.table(t(as.data.frame(out)), file = filename, row.names = FALSE, col.names = FALSE)

