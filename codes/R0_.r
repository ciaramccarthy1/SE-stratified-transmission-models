# R0 and beta
# COVID-19, Influenza
#
# Requires cm45, pars


## Demography
na   <- pars$na
nimd <- pars$nimd



## Average contact rate of cm45
#   contact rate of participant across contacts
cp <- vector(); for (i in 1:(na*nimd)){ cp[i]=sum(cm45[i,])}  #cp #[1] 8.053193 14.369821 15.171069 11.534715 10.880512 10.948698  8.581716  7.393931  5.250907  8.945841 ...
#   population proportion across age x imd strata
pa0 <- vector()
for (is in 1:nimd) {
  pa0[(is-1)*na + 1:na] =   demog2021$Proportion[1:na + na*(is-1) + na*nimd*(1-urb)] } #length(pa0) #[1] 45
pa0 <- pa0/sum(pa0)
#   Check:
#     sum(pa0 - demog2021$Proportion[which(demog2021$rural=="Urban")]/sum(demog2021$Proportion[which(demog2021$rural=="Urban")]))
#     [1] 0 
#   average contact rate over participants
cav = sum(pa0*cp)
print(paste0("Average contact rate of cm45: ", round(cav,4), "/day")) #[1] 10.80047



## R0 and beta #################################################################

# Applied directly to cm45
#   eigen(cm45)$values       #[1] 46.054457390 25.056689207 21.525676962 10.116110636  9.056574541  7.461609085  5.418768558  4.596350633
#                            [41] 0.053084981  0.038487281  0.020231916  0.011700634  0.008343733
#   max(eigen(cm45)$values) # [1] 46.05446

## NGM
#Erlang has no effect on R0 without vital dynamics (mu=0)
#https://www.biorxiv.org/content/10.1101/319574v1.full.pdf

R0 <- function(pars,R0assumed=2.5,printout=1){
  
ng    = na*nimd
u45   = rep(pars$u,nimd)
y45   = pars$y45
orIR  = 1/pars$rIR          
orUR  = 1/pars$rUR          
beta0 = pars$beta            
fu    = pars$f              
ngm   = cm45

for (k in 1:ng){ 
  y_k=y45[k]
  for (j in 1:ng) {
    ngm[j,k] = beta0*u45[j]*cm45[j,k]*( y_k*orIR + fu*(1-y_k)*orUR) }}

# max EV
EVs = eigen(ngm)$values
R00 = max(Re(EVs[which(Im(EVs)==0)]))

# beta given R0assumed
beta = R0assumed/(R00/beta0)
if(printout==1) print(paste0("Assuming R0 = ", R0assumed, ", then beta = ", round(beta,4), " /day"))

return(beta)
}


#beta = R0(pars)

