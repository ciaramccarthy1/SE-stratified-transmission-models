################################################################################
# Parameter-uncertainty ensemble: propagate David's JOINT posterior through our
# model. For each posterior draw we derive every posterior-dependent parameter
# (beta=qp, f=alpha_i, seasonal b1/phi/psi, waning 1/om, latency 1/si, and the
# state-weighted u / y / rIR / rUR), rebuild the contact matrix from that draw's
# qc, set the initial conditions from that draw's post-burn-in state, run the
# cpp model, and collect weekly symptomatic incidence (total and by age).
# Output: median + 95% credible-interval ribbons vs David's target.
#
# Usage:  Rscript codes/uncertainty_ensemble.R [n_draws]   (default: all draws)
################################################################################
suppressPackageStartupMessages({library(here);library(Rcpp);library(magrittr);library(tidyverse)})
setwd(here())

## --- inputs -----------------------------------------------------------------
post   <- read.csv("data/david_posteriors.csv")
states <- read.csv("data/david_states_by_draw.csv")
# NB: named cnt_p/cnt_c (not cp/cc) — R0_.r overwrites globals `cp` mid-run.
cnt_p  <- as.matrix(read.csv("data/david_cnt_p25.csv")); dimnames(cnt_p)<-NULL
cnt_c  <- as.matrix(read.csv("data/david_cnt_c25.csv")); dimnames(cnt_c)<-NULL
pop25  <- read.csv("data/david_pop25.csv")$pop
ons    <- read.csv("data/demographics_9age.csv")               # ONS IMD shares for contact expansion
ages9  <- c("0 to 4","5 to 14","15 to 24","25 to 34","35 to 44","45 to 54","55 to 64","65 to 74","75+")
na<-9; nimd<-5; ng<-na*nimd

## band geometry (25 -> 9) --------------------------------------------------
d_lo<-c((0:11)/12,1,2,3,4,5,10,15,25,35,45,55,65,75); d_hi<-c((1:12)/12,2,3,4,5,10,15,25,35,45,55,65,75,90)
o_lo<-c(0,5,15,25,35,45,55,65,75); o_hi<-c(5,15,25,35,45,55,65,75,90)
ovl <-function(b,j) max(0,min(o_hi[j],d_hi[b])-max(o_lo[j],d_lo[b]))          # overlap in years
fr  <-outer(1:25,1:9,Vectorize(function(b,j) ovl(b,j)/(d_hi[b]-d_lo[b])))     # width-frac of band b in target j
## contact-aggregation weights (fixed across draws; only C25 varies with qc)
P    <-outer(1:25,1:9,Vectorize(function(b,i) pop25[b]*ovl(b,i)/(d_hi[b]-d_lo[b])))
Nmat <-sapply(1:na,function(i) sapply(1:nimd,function(s) ons$Population[ons$IMD==s & ons$Age==ages9[i]]))
share<-sweep(Nmat,2,colSums(Nmat),"/")                                        # 5 x 9 IMD share by age
colP <-colSums(P)

## --- derive scenario_overrides for one posterior draw ----------------------
derive <- function(d){
  pp <- post[post$draw==d,]
  st <- states[states$draw==d,]; st <- st[order(st$age_group),]
  S<-as.matrix(st[,c("S0","S1","S2","S3")]); E<-as.matrix(st[,c("E0","E1","E2","E3")])
  A<-as.matrix(st[,c("A0","A1","A2","A3")]); I<-as.matrix(st[,c("I0","I1","I2","I3")]); R<-as.matrix(st[,c("R0","R1","R2","R3")])
  sig <- c(1, pp$d1, pp$d1*pp$d2, pp$d1*pp$d2*pp$d3)                          # exposure-group susceptibility
  Ssum<-rowSums(S); Isum<-rowSums(I); Asum<-rowSums(A); Ntot<-rowSums(S+E+A+I+R)
  u25 <- as.numeric((S%*%sig)/Ssum)
  y25 <- 1 - c(rep(pp$pA1,12),rep(pp$pA2,4),rep(pp$pA3,2),rep(pp$pA4,7))       # 1-pA by age band
  rate<- 1/(pp$g0 * c(1, pp$g1, pp$g1*pp$g2, pp$g1*pp$g2))                     # recovery rate by exposure group (grp3=grp2)
  # aggregate 25 -> 9, each with its epidemiologically-correct weight
  wsum <- function(x,w){ v<-colSums(fr*(w*x)); d<-colSums(fr*w); ifelse(d>0, v/d, 0) }
  u9  <- wsum(u25, Ssum)                                                       # susceptible-weighted
  y9  <- wsum(y25, Isum+Asum)                                                  # infection(prevalence)-weighted
  # rIR/rUR: weight group rates by infectious counts, aggregated over (band,group) in each 9-band
  wrate<-function(W){ num<-colSums(fr * as.numeric(W%*%rate)); den<-colSums(fr*rowSums(W)); ifelse(den>0,num/den,mean(rate)) }
  rIR9<-wrate(I); rUR9<-wrate(A)
  # per-draw contact matrix: C25 = cnt_p + qc*cnt_c, aggregate to 9, IMD-expand
  C25 <- cnt_p + pp$qc*cnt_c
  Cna <- sweep(t(P)%*%C25%*%fr, 1, colP, "/")
  Cout<- matrix(0,ng,ng)
  for(s1 in 1:nimd) for(i in 1:na) for(s2 in 1:nimd) for(j in 1:na) Cout[(s1-1)*na+i,(s2-1)*na+j]<-Cna[i,j]*share[s2,j]
  # initial-condition states (25-band df in the init_conditions_allgroups format)
  ic <- data.frame(age_group=1:25, Ntot=Ntot,
                   frac_R=rowSums(R)/Ntot, frac_E=rowSums(E)/Ntot,
                   frac_A=rowSums(A)/Ntot, frac_I=rowSums(I)/Ntot)
  list(beta_override=pp$qp, f=pp$alpha_i, b1=pp$b1, phi=pp$phi, psi=pp$psi,
       rW_nat=1/pp$om, rEI=1/pp$si, u=u9, y45=rep(y9,nimd), rIR=rIR9, rUR=rUR9,
       cm_override=Cout, ic_states=ic,
       times=seq(0,364), nt=3641, nd=365)
}

## --- draw-1 sanity check (should be near parsR medians) ---------------------
cat("=== derived params, draw 1 (should be near parsR medians) ===\n")
ov1<-derive(post$draw[1])
cat("u  :",round(ov1$u,3),"\n"); cat("y  :",round(ov1$y45[1:9],3),"\n")
cat("rIR:",round(ov1$rIR,3),"\n"); cat(sprintf("beta=%.4f f=%.3f b1=%.3f phi=%.3f psi=%.3f rW=1/%.0f\n\n",
    ov1$beta_override,ov1$f,ov1$b1,ov1$phi,ov1$psi,1/ov1$rW_nat))

## --- run the ensemble -------------------------------------------------------
args<-commandArgs(trailingOnly=TRUE); ND<-if(length(args)) as.integer(args[1]) else nrow(post)
draws<-post$draw[1:ND]
source("codes/setup.r")
pset$Disease<-"RSV-illness";pset$Vaccination<-0;pset$DailyIncidence<-1;pset$Incidence<-"Daily"
pset$Namevacc<-"";pset$TODAY<-"";pset$FIGURES<-0;pset$DIAGNOSTIC<-0;pset$SUMMARY<-0;pset$COMPILE<-1
pset$DavidContacts<-TRUE; pset$DavidDemog<-TRUE                    # cm/ic overridden per draw; demog = David pop

wk_of<-function(v){ w<-tapply(v, ceiling(seq_along(v)/7), sum); as.numeric(w[as.integer(names(w))<=52]) }
tot_mat<-matrix(NA,52,ND); age_arr<-array(NA,dim=c(52,na,ND))
for(k in seq_len(ND)){
  scenario_overrides <- derive(draws[k])
  ok<-tryCatch({invisible(capture.output(source("codes/modelrun.r"))); TRUE}, error=function(e){cat("draw",draws[k],"failed:",conditionMessage(e),"\n");FALSE})
  if(!ok) next
  tot_mat[,k]<-wk_of(as.numeric(mas$byw$Iw))
  for(a in 1:na) age_arr[,a,k]<-wk_of(as.numeric(mas$byaw$Iw_a[,a]))
  if(k%%10==0) cat("ran",k,"/",ND,"draws\n")
}
rm(scenario_overrides)

## --- David: all sample runs -> ribbons (fallback to mean CSV as 1 run) -------
W9<-outer(1:25,1:9,Vectorize(function(b,j) ovl(b,j)/(d_hi[b]-d_lo[b])))
qs<-function(m) apply(matrix(m,nrow=52),1,quantile,probs=c(.025,.5,.975),na.rm=TRUE)
lab<-c("0-4","5-14","15-24","25-34","35-44","45-54","55-64","65-74","75+")
mk<-function(q,src,band=NA) data.frame(week=1:52,lo=q[1,],med=q[2,],hi=q[3,],source=src,band=band)

dav_file<-"data/david_outcomes_by_sample.csv"
if(file.exists(dav_file)){
  dvs<-read.csv(dav_file); if("outcome"%in%names(dvs)) dvs<-dvs[dvs$outcome=="symptomatic",]
  dvs$week<-dvs$week_no-104; dvs<-dvs[dvs$week>=1 & dvs$week<=52,]
  ss<-sort(unique(dvs$s)); DT<-matrix(0,52,length(ss)); DA<-array(0,dim=c(52,9,length(ss)))
  for(i in seq_along(ss)){ d<-dvs[dvs$s==ss[i],]; M<-matrix(0,52,25); M[cbind(d$week,d$age_group)]<-d$cases
    A<-M%*%W9; DA[,,i]<-A; DT[,i]<-rowSums(A) }
  cat(sprintf("David: %d sample runs loaded\n",length(ss)))
} else {
  dd<-read.csv("data/no_vacc_weekly_by_age_outcome.csv")|>filter(outcome=="symptomatic")|>
    select(age_group,week_no,cases)|>pivot_wider(names_from=age_group,values_from=cases)|>arrange(week_no)|>mutate(week_no=week_no-104)
  d0<-as.matrix(dd[,as.character(1:25)])%*%W9; DT<-matrix(rowSums(d0),52,1); DA<-array(d0,dim=c(52,9,1))
  cat("David: per-sample file not found; using mean CSV as a single run (points-equivalent)\n")
}

## --- our model = 95% CrI band; David = every rsvie run as a line ------------
linesdf<-function(M,band=NA){ M<-matrix(M,nrow=52); df<-as.data.frame(M)
  names(df)<-paste0("r",seq_len(ncol(df))); df$week<-1:52
  pivot_longer(df,-week,names_to="rep",values_to="cases")|>mutate(band=band) }
ribdf <-function(M,band=NA){ q<-qs(M); data.frame(week=1:52,lo=q[1,],med=q[2,],hi=q[3,],band=band) }

Td<-linesdf(DT); Tr<-ribdf(tot_mat)
pT<-ggplot()+
  geom_line(data=Td,aes(week,cases,group=rep),colour="#5e81ac",alpha=.30,linewidth=.32)+
  geom_ribbon(data=Tr,aes(week,ymin=lo,ymax=hi),fill="grey55",alpha=.45)+
  geom_line(data=Tr,aes(week,med),colour="grey10",linewidth=.9)+
  labs(x="Week",y="Weekly total symptomatic",
    title="Total symptomatic: our model 95% CrI (grey band) vs David's rsvie runs (blue lines)",
    subtitle=sprintf("%d model posterior draws | %d David sample runs",ncol(tot_mat),ncol(DT)))+
  theme_minimal(base_size=11)
ggsave("output/uncertainty_total.png",pT,width=8,height=5,dpi=140)

ageL<-do.call(rbind,lapply(1:na,function(a) linesdf(DA[,a,],lab[a])))
ageR<-do.call(rbind,lapply(1:na,function(a) ribdf(age_arr[,a,],lab[a])))
ageL$band<-factor(ageL$band,levels=lab); ageR$band<-factor(ageR$band,levels=lab)
pA<-ggplot()+
  geom_line(data=ageL,aes(week,cases,group=interaction(band,rep)),colour="#5e81ac",alpha=.30,linewidth=.26)+
  geom_ribbon(data=ageR,aes(week,ymin=lo,ymax=hi),fill="grey55",alpha=.45)+
  geom_line(data=ageR,aes(week,med),colour="grey10",linewidth=.6)+
  facet_wrap(~band,scales="free_y",ncol=3)+
  labs(x="Week",y="Weekly symptomatic",
    title="By age: our model 95% CrI (grey band) vs David's rsvie runs (blue lines), David demographics")+
  theme_minimal(base_size=10)
ggsave("output/uncertainty_byage.png",pA,width=10,height=8.5,dpi=140)

cat("\n=== season-total by band: ours vs David, median [2.5%, 97.5%] (thousands) ===\n")
for(a in 1:na){
  o<-colSums(matrix(age_arr[,a,],nrow=52)); d<-colSums(matrix(DA[,a,],nrow=52))
  ovr<-if(quantile(o,.975)<quantile(d,.025)) "David ABOVE (disjoint)" else if(quantile(o,.025)>quantile(d,.975)) "ours ABOVE (disjoint)" else "overlap"
  cat(sprintf("%-6s ours %5.0f [%5.0f,%5.0f] | David %5.0f [%5.0f,%5.0f]  -> %s\n", lab[a],
    median(o)/1e3, quantile(o,.025)/1e3, quantile(o,.975)/1e3,
    median(d)/1e3, quantile(d,.025)/1e3, quantile(d,.975)/1e3, ovr))
}
cat(sprintf("\nSeason-total symptomatic: ours %.2fM [%.2f, %.2f] | David %.2fM [%.2f, %.2f]\n",
  median(colSums(tot_mat),na.rm=TRUE)/1e6, quantile(colSums(tot_mat),.025,na.rm=TRUE)/1e6, quantile(colSums(tot_mat),.975,na.rm=TRUE)/1e6,
  median(colSums(DT))/1e6, quantile(colSums(DT),.025)/1e6, quantile(colSums(DT),.975)/1e6))
cat("Wrote output/uncertainty_total.png and output/uncertainty_byage.png\n")
