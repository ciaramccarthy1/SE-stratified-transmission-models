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

## expand na -> na*5, homogeneous across IMD (IMD-major: g=(s-1)*na+i)
demog <- read.csv("data/demographics_9age.csv")
ages  <- c("0 to 4","5 to 14","15 to 24","25 to 34","35 to 44","45 to 54","55 to 64","65 to 74","75+")
Nmat  <- matrix(0,5,na); for(s in 1:5) for(i in 1:na) Nmat[s,i] <- demog$Population[demog$IMD==s & demog$Age==ages[i]]
share <- sweep(Nmat, 2, colSums(Nmat), "/")
ng  <- na*5
Cout <- matrix(0,ng,ng)
for(s1 in 1:5) for(i in 1:na) for(s2 in 1:5) for(j in 1:na)
  Cout[(s1-1)*na+i, (s2-1)*na+j] <- Cna[i,j]*share[s2,j]

write.table(Cout, "data/Mas45_david.csv", sep=",", row.names=FALSE, col.names=FALSE)
cat("Wrote data/Mas45_david.csv (", nrow(Cout), "x", ncol(Cout), ")\n")
cat("Cna row sums (contacts/person by", na, "bands):\n"); print(round(rowSums(Cna),2))

## David demographics at our na bands (his stationary populationAgeGroup, from
## david_pop25), split across IMD by ONS shares so the age TOTALS are David's while
## the model keeps its IMD structure. For REPRODUCING David's results only:
## set pset$DavidDemog <- TRUE (default FALSE uses real ONS demographics_9age.csv).
david_age <- sapply(david_bands_in, function(idx) sum(pop[idx]))  # David pop per our band
dpop      <- sweep(share, 2, david_age, "*")  # 5 x na: David age totals, ONS IMD split
tot       <- sum(dpop)
david_demog <- do.call(rbind, lapply(1:5, function(s) data.frame(
  Age = ages, IMD = s, Population = as.numeric(dpop[s, ]),
  tot_pop = tot, Proportion = as.numeric(dpop[s, ])/tot)))
write.csv(david_demog, "data/demographics_9age_david.csv", row.names=FALSE)
cat("Wrote data/demographics_9age_david.csv (David pops, ONS IMD split)\n")
cat("Pop by band - David vs ONS:\n"); print(round(rbind(David=david_age, ONS=colSums(Nmat))))
