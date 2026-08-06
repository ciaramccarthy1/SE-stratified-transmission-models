################################################################################
# Build a David-rsvie contact matrix for our 9-age x 5-IMD model.
#
# Inputs (exported from rsvie into data/):
#   david_cnt_p25.csv, david_cnt_c25.csv : 25x25 David TRANSMISSION matrices
#                                          cnt_matrix_p / cnt_matrix_c (household
#                                          infant contacts included) 
#   david_pop25.csv                       : populationAgeGroup (25 bands)
# Steps:
#   1. Effective contacts C25 = cnt_p + qc*cnt_c  (qp = overall scale -> our beta0)
#   2. Aggregate 25 -> our 9 bands: contact side = sum, participant side =
#      population-weighted mean. Our bands are exact unions of David's (no straddling).
#   3. Expand na -> na*5, HOMOGENEOUS across IMD:
#        Cout[(s1,i),(s2,j)] = Cna[i,j] * N[s2,j]/N[j]    (our demographics)
#      i.e. age-i contacts split across IMD by population share, independent of
#      the participant's own IMD. In the FOI this collapses to age-only mixing.
# Output: data/Mas45_david.csv  (use via pset$DavidContacts <- TRUE)
################################################################################

qc  <- 0.999182                       # median posterior conversational weight
# David's TRANSMISSION matrices cnt_matrix_p/c (household infant contacts included),
# NOT the raw survey contactMatrixPhy/Con (~0 infant contacts, inflated adult contacts).
cp  <- as.matrix(read.csv("data/david_cnt_p25.csv"))
cc  <- as.matrix(read.csv("data/david_cnt_c25.csv"))
pop <- read.csv("data/david_pop25.csv")$pop
stopifnot(dim(cp) == c(25, 25), dim(cc) == c(25, 25), length(pop) == 25,
          is.numeric(cp), is.numeric(cc), is.numeric(pop))
C25 <- cp + qc*cc; dimnames(C25) <- NULL   # David's raw transmission matrix (no symmetrisation)

## --- Age-band edges (years) -------------------------------------------------
## David has 25 bands: 12 monthly (0-1yr), yearly 1-4, 5-10, 10-15, then decadal.
## We have 9. Our bands are EXACT UNIONS of David's (the edges coincide), so each
## David band lies wholly inside exactly one of our bands - none straddle a
## boundary. This makes the collapse exact (no interpolation).
na   <- 9
d_lo <- c((0:11)/12, 1,2,3,4, 5,10, 15,25,35,45,55,65,75)   # David 25 band lower edges
d_hi <- c((1:12)/12, 2,3,4,5, 10,15, 25,35,45,55,65,75,90)  # David 25 band upper edges
o_lo <- c(0,5,15,25,35,45,55,65,75); o_hi <- c(5,15,25,35,45,55,65,75,90)  # our 9 band edges

## --- Collapse the 25x25 matrix to our 9x9 -----------------------------------
## A contact-matrix entry C[a,b] = mean contacts a person in band a makes with
## band b (per person, per unit time). The two axes aggregate DIFFERENTLY:
##   * contacted axis -> SUM: one person's contacts with a merged group equal the
##       sum of their contacts with each of its sub-bands.
##   * participant axis -> POPULATION-WEIGHTED MEAN: the average contacts of a
##       merged group is the mean over its sub-bands, weighted by their sizes.

# david_bands_in[[i]] = indices of David's 25 bands that fall inside our band i.
# (Our bands are exact unions of David's, so each David band belongs to exactly one.)
david_bands_in <- lapply(1:na, function(i) which(d_lo >= o_lo[i] & d_hi <= o_hi[i]))
stopifnot(length(unlist(david_bands_in)) == 25,     # every David band is used...
          !anyDuplicated(unlist(david_bands_in)))   # ...exactly once (bands nest cleanly)

Cna <- matrix(0, na, na)
for (i in 1:na) {                      # our PARTICIPANT band
  A <- david_bands_in[[i]]             #   the David bands making it up
  for (j in 1:na) {                    # our CONTACT band
    B <- david_bands_in[[j]]           #   the David bands making it up
    # for each participant sub-band A: SUM its contacts across the contact sub-bands B
    contacts_per_A <- rowSums(C25[A, B, drop = FALSE])
    # then take the population-weighted MEAN over the participant sub-bands
    Cna[i, j] <- sum(pop[A] * contacts_per_A) / sum(pop[A])
  }
}

## expand 9-band Cna -> 45x45, homogeneous across IMD (IMD-major: g=(s-1)*na+i).
## Read the 10-band ONS demographics; fold the top two (75-84, 85+) into a 9-band 75+
## for the expansion, and keep their split for the 75+ contact-band split below.
demog  <- read.csv("data/demographics_10age.csv")
ages10 <- c("0 to 4","5 to 14","15 to 24","25 to 34","35 to 44","45 to 54","55 to 64","65 to 74","75 to 84","85+")
N10 <- matrix(0,5,10); for(s in 1:5) for(i in 1:10) N10[s,i] <- demog$Population[demog$IMD==s & demog$Age==ages10[i]]
N9  <- cbind(N10[,1:8], N10[,9] + N10[,10])       # 5 x 9 (band 9 = 75-84 + 85+)
share9 <- sweep(N9, 2, colSums(N9), "/")          # IMD share within each 9-band
ng9  <- na*5
Cout <- matrix(0,ng9,ng9)
for(s1 in 1:5) for(i in 1:na) for(s2 in 1:5) for(j in 1:na)
  Cout[(s1-1)*na+i, (s2-1)*na+j] <- Cna[i,j]*share9[s2,j]

## Split the 75+ contact band (9) into 75-84 (9) & 85+ (10) -> 50x50, as for Reconnect:
## participant side duplicates the 75+ row; contact side splits the 75+ column by the
## ONS 75-84 vs 85+ population share (per contact IMD). David has no data finer than 75+.
sh9  <- N10[,9]  / (N10[,9] + N10[,10])            # share of 75+ that is 75-84, per IMD
sh10 <- N10[,10] / (N10[,9] + N10[,10])
na10 <- 10; ng10 <- na10*5
Cd10 <- matrix(0, ng10, ng10)
for(s1 in 1:5) for(i in 1:na10) for(s2 in 1:5) for(j in 1:na10) {
  i9 <- min(i,9); j9 <- min(j,9)                  # bands 9,10 -> old 75+ (band 9)
  fac <- if (j==9) sh9[s2] else if (j==10) sh10[s2] else 1
  Cd10[(s1-1)*na10+i, (s2-1)*na10+j] <- Cout[(s1-1)*na+i9, (s2-1)*na+j9] * fac
}
write.table(Cd10, "data/Mas50_david.csv", sep=",", row.names=FALSE, col.names=FALSE)
cat("Wrote data/Mas50_david.csv (", nrow(Cd10), "x", ncol(Cd10), ")\n")
cat("Cna row sums (contacts/person by 9 bands, pre-split):\n"); print(round(rowSums(Cna),2))

## David demographics at our 10 bands: David's age totals (his populationAgeGroup),
## the 75+ total split into 75-84/85+ by the overall ONS share, then ONS IMD split.
## For REPRODUCING David's results only (pset$DavidDemog <- TRUE).
david_age9  <- sapply(david_bands_in, function(idx) sum(pop[idx]))  # David pop per 9-band
ons_sh9     <- sum(N10[,9]) / (sum(N10[,9]) + sum(N10[,10]))        # overall ONS 75-84 share
david_age10 <- c(david_age9[1:8], david_age9[9]*ons_sh9, david_age9[9]*(1-ons_sh9))
share10 <- sweep(N10, 2, colSums(N10), "/")
dpop    <- sweep(share10, 2, david_age10, "*")    # 5 x 10: David age totals, ONS IMD split
tot     <- sum(dpop)
david_demog <- do.call(rbind, lapply(1:5, function(s) data.frame(
  Age = ages10, IMD = s, Population = as.numeric(dpop[s, ]),
  tot_pop = tot, Proportion = as.numeric(dpop[s, ])/tot)))
write.csv(david_demog, "data/demographics_10age_david.csv", row.names=FALSE)
cat("Wrote data/demographics_10age_david.csv (David pops, ONS IMD split)\n")
cat("David pop by 10-band:\n"); print(round(david_age10))
