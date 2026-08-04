################################################################################
# Prepare derived 9-age-band model inputs from raw downloads.
# Pairs with codes/fetch_data.R (which downloads the raw files).
#
# Model age bands (aligned to David rsvie adult bands so his 25 bands and
# Reconnect's 5-year bands both nest cleanly):
#   0-4, 5-14, 15-24, 25-34, 35-44, 45-54, 55-64, 65-74, 75+
#
# Inputs (all sourced by fetch_data.R):
#   - data/base_matrix.csv               : Reconnect long-form 5-IMD x 16-age
#                                          (5 x 5 IMD pairs x 16 x 16 = 6400 rows)
#   - data/ons_lsoa_syoa_2022-2024.xlsx  : ONS LSOA population by SYA (Mid-2024 sheet)
#   - data/iod2025_lsoa_ranks_deciles.csv: IoD 2025 LSOA -> IMD decile lookup
#
# Outputs:
#   - data/Mas45.csv             : 45x45 wide contact matrix
#                                  (rows = participant ig = is*na + ia)
#   - data/demographics_9age.csv : 9-band x 5-IMD-quintile demographics
#
# Aggregation notes for the contact matrix:
#   Participant-side merge of {a1, a2} -> A: population-weighted average
#     (one participant is in exactly one sub-band, so we average rates).
#   Contact-side merge of {b1, b2} -> B: SUM
#     (a participant in A makes contacts with both sub-bands; counts add).
#   TODO: currently using UNIFORM weights on the participant side
#     (i.e. simple mean over sub-bands). Replace with population-weighted
#     means using ONS LSOA SYA data (already fetched by fetch_data.R).
################################################################################

input_dir  <- file.path("data")
output_dir <- file.path("data")

# Reconnect 16 bands 
reconnect_ages <- c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29",
               "30-34", "35-39", "40-44", "45-49", "50-54", "55-59",
               "60-64", "65-69", "70-74", "75+")
# 9 bands (aligned with David's) - each is an exact aggregation of Reconnect bands
new_ages  <- c("0-4", "5-14", "15-24", "25-34", "35-44",
               "45-54", "55-64", "65-74", "75+")
# Mapping: new_age -> the reconnect_ages that compose it
new_to_reconnect <- list(
  "0-4"   = "0-4",
  "5-14"  = c("5-9", "10-14"),
  "15-24" = c("15-19", "20-24"),
  "25-34" = c("25-29", "30-34"),
  "35-44" = c("35-39", "40-44"),
  "45-54" = c("45-49", "50-54"),
  "55-64" = c("55-59", "60-64"),
  "65-74" = c("65-69", "70-74"),
  "75+"   = "75+")

## --- 0. Population by IMD quintile (shared input) ----------------------------
## Load ONS single-year-of-age (SYOA) LSOA populations and the IoD 2025 decile
## lookup ONCE, join, and aggregate to (quintile x SYOA). Both the contact-matrix
## weights (pop_recon, 16-band) and the demographics (9-band, added in a later
## step) derive from this, so the ONS file is read a single time.
##   Decile -> Quintile: quintile = ceiling(decile / 2)  (1-2 -> Q1 most deprived)

# Raw inputs are fetched by codes/fetch_data.R.
source(file.path("codes", "fetch_data.R"))
suppressPackageStartupMessages(library(readxl))

# ONS: row 4 of the Mid-2024 sheet is the header (rows 1-3 are notes / merged title).
ons <- read_excel(file.path(input_dir, "ons_lsoa_syoa_2022-2024.xlsx"),
                  sheet = "Mid-2024 LSOA 2021", skip = 3)
syoa_ages <- 0:90                                    # F0..F90 / M0..M90 (90 = 90+)
pop_syoa  <- sapply(syoa_ages, function(a)
  as.numeric(ons[[paste0("F", a)]]) + as.numeric(ons[[paste0("M", a)]]))
colnames(pop_syoa) <- as.character(syoa_ages)

iod <- read.csv(file.path(input_dir, "iod2025_lsoa_ranks_deciles.csv"),
                header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
decile_col <- grep("^Index of Multiple Deprivation \\(IMD\\) Decile",
                   names(iod), value = TRUE)[1]
lsoa_imd <- data.frame(lsoa     = iod[["LSOA code (2021)"]],
                       decile   = iod[[decile_col]],
                       quintile = ceiling(iod[[decile_col]] / 2),
                       stringsAsFactors = FALSE)

# Inner join SYOA populations to IMD (England-only LSOAs present in both)
joined <- merge(data.frame(lsoa = ons[["LSOA 2021 Code"]], pop_syoa, check.names = FALSE),
                lsoa_imd, by = "lsoa")
cat(sprintf("Joined %d LSOAs (England) on (LSOA 2021 code)\n", nrow(joined)))

syoa_mat   <- as.matrix(joined[, as.character(syoa_ages)])
pop_q_syoa <- rowsum(syoa_mat, joined$quintile)      # 5 x 91  (quintile x SYOA)
pop_d_syoa <- rowsum(syoa_mat, joined$decile)        # 10 x 91 (decile x SYOA)

# Aggregate single years into a set of bands (each band = a run of single years)
agg_bands <- function(pmat, def)
  sapply(def, function(a) rowSums(pmat[, as.character(a), drop = FALSE]))

# 16 Reconnect bands -> participant-side contact weights
recon_def <- list(
  "0-4"=0:4, "5-9"=5:9, "10-14"=10:14, "15-19"=15:19, "20-24"=20:24,
  "25-29"=25:29, "30-34"=30:34, "35-39"=35:39, "40-44"=40:44, "45-49"=45:49,
  "50-54"=50:54, "55-59"=55:59, "60-64"=60:64, "65-69"=65:69, "70-74"=70:74,
  "75+"=75:90)
pop_recon <- agg_bands(pop_q_syoa, recon_def)        # 5 quintile x 16 Reconnect band

# 9 model bands -> demographics
band_def <- list(
  "0 to 4"=0:4, "5 to 14"=5:14, "15 to 24"=15:24, "25 to 34"=25:34,
  "35 to 44"=35:44, "45 to 54"=45:54, "55 to 64"=55:64, "65 to 74"=65:74,
  "75+"=75:90)
band_names <- names(band_def)
pop_band9  <- agg_bands(pop_q_syoa, band_def)        # 5 quintile x 9 model band

## --- 1. Contact matrix ------------------------------------------------------

# Raw inputs are fetched by codes/fetch_data.R.
source(file.path("codes", "fetch_data.R"))

reconnect_path <- file.path(input_dir, "base_matrix.csv")
cm_long <- read.csv(reconnect_path, header = TRUE, stringsAsFactors = FALSE)
stopifnot(nrow(cm_long) == 5 * 5 * length(reconnect_ages) * length(reconnect_ages))

n_imd_p <- length(unique(cm_long$Participant_IMD))
n_imd_c <- length(unique(cm_long$Contact_IMD))
stopifnot(n_imd_p == 5, n_imd_c == 5)
na_new  <- length(new_ages)   # 9
nimd    <- 5
ng_new  <- na_new * nimd      # 45

key <- function(pi, ci, pa, ca) paste(pi, ci, pa, ca, sep = "|")
mean_lookup <- setNames(cm_long$mean,
                        with(cm_long, key(Participant_IMD, Contact_IMD,
                                          Participant_age_group, Contact_age_group)))

cm50 <- matrix(0, ng_new, ng_new)
for (is_p in 1:nimd) {            # participant IMD
  for (is_c in 1:nimd) {          # contact IMD
    for (ia_p in 1:na_new) {      # new participant age
      for (ia_c in 1:na_new) {    # new contact age
        sub_p <- new_to_reconnect[[new_ages[ia_p]]]
        sub_c <- new_to_reconnect[[new_ages[ia_c]]]
        # contact-side: SUM (a participant makes contacts with both sub-bands)
        sum_over_contacts_per_subp <- vapply(
          sub_p, function(pa) {
            sum(vapply(sub_c, function(ca)
              mean_lookup[[key(is_p, is_c, pa, ca)]],
              numeric(1)))
          }, numeric(1))
        # participant-side: POPULATION-WEIGHTED mean over participant sub-bands,
        # using this IMD's sub-band populations (pop_recon[is_p, sub_p]).
        w <- pop_recon[is_p, sub_p]
        cm50[(is_p - 1)*na_new + ia_p,
             (is_c - 1)*na_new + ia_c] <- sum(w * sum_over_contacts_per_subp) / sum(w)
      }
    }
  }
}

write.table(cm50, file.path(output_dir, "Mas45.csv"),
            sep = ",", row.names = FALSE, col.names = FALSE)
cat(sprintf("Wrote %s: %d x %d (from Reconnect base_matrix.csv)\n",
            "data/Mas45.csv", nrow(cm50), ncol(cm50)))

## --- 2. Demographics (9-band x IMD quintile) --------------------------------
## Built from the shared quintile population computed in Section 0.
demog10 <- do.call(rbind, lapply(1:5, function(q) {
  pops <- pop_band9[q, ]
  tot  <- sum(pops)
  data.frame(
    Age        = band_names,
    IMD        = q,
    Population = as.numeric(pops),
    tot_pop    = tot,
    Proportion = as.numeric(pops) / tot,
    stringsAsFactors = FALSE)
}))
rownames(demog10) <- NULL

write.csv(demog10, file.path(output_dir, "demographics_9age.csv"),
          row.names = FALSE)
cat(sprintf("Wrote %s: %d rows (expected %d)\n",
            "data/demographics_9age.csv",
            nrow(demog10), 5 * length(band_names)))

## --- 3. RSV uptake by IMD quintile (population-weighted) -------------------
## Collapse decile -> quintile using the 75+ population in each decile as
## weights, since the programme targets 75+. NOTE: will need to edit if looking at 65-74yos

uptake_path <- file.path(input_dir, "rsv_uptake_by_imd_decile.csv")
uptake_dec  <- read.csv(uptake_path, stringsAsFactors = FALSE)
uptake_dec  <- uptake_dec[order(uptake_dec$decile), ]
stopifnot(all(uptake_dec$decile == 1:10))

# 75+ population by decile (sum of single years 75..90, from the shared Section 0 pop)
pop_75_by_decile <- setNames(rowSums(pop_d_syoa[, as.character(75:90)]),
                             rownames(pop_d_syoa))
stopifnot(length(pop_75_by_decile) == 10,
          all(as.integer(names(pop_75_by_decile)) == 1:10))

# Population-weighted collapse: paired deciles (2k-1, 2k) -> quintile k
uptake_dec$pop_75    <- as.numeric(pop_75_by_decile[as.character(uptake_dec$decile)])
uptake_dec$quintile  <- ceiling(uptake_dec$decile / 2)
uptake_quintile <- with(uptake_dec,
                        tapply(uptake_pct * pop_75, quintile, sum) /
                        tapply(pop_75,              quintile, sum))

uptake_qntl_df <- data.frame(
  quintile   = 1:5,
  uptake_pct = as.numeric(uptake_quintile))
write.csv(uptake_qntl_df,
          file.path(output_dir, "rsv_uptake_by_imd_quintile.csv"),
          row.names = FALSE)
cat(sprintf("Wrote %s: %d quintiles (75+-population-weighted from decile)\n",
            "data/rsv_uptake_by_imd_quintile.csv", nrow(uptake_qntl_df)))


