################################################################################
# Prepare derived 10-age-band model inputs from raw downloads.
# Pairs with codes/fetch_data.R (which downloads the raw files).
#
# Model age bands: 0-4, 5-14, 15-19, 20-29, 30-39, 40-49, 50-59, 60-64, 65-74, 75+
#
# Inputs (all sourced by fetch_data.R):
#   - data/base_matrix.csv               : Reconnect long-form 5-IMD x 16-age
#                                          (5 x 5 IMD pairs x 16 x 16 = 6400 rows)
#   - data/ons_lsoa_syoa_2022-2024.xlsx  : ONS LSOA population by SYA (Mid-2024 sheet)
#   - data/iod2025_lsoa_ranks_deciles.csv: IoD 2025 LSOA -> IMD decile lookup
#
# Outputs:
#   - data/Mas50.csv             : 50x50 wide contact matrix
#                                  (rows = participant ig = is*na + ia)
#   - data/demographics_10age.csv: 10-band x 5-IMD-quintile demographics
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

# Reconnect 16 bands (sorted in epidemiologically meaningful order)
reconnect_ages <- c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29",
               "30-34", "35-39", "40-44", "45-49", "50-54", "55-59",
               "60-64", "65-69", "70-74", "75+")
# New 10 bands - each is an exact aggregation of one or two of Reconnect bands
new_ages  <- c("0-4", "5-14", "15-19", "20-29", "30-39", "40-49",
               "50-59", "60-64", "65-74", "75+")
# Mapping: new_age -> the reconnect_ages that compose it
new_to_reconnect <- list(
  "0-4"   = "0-4",
  "5-14"  = c("5-9", "10-14"),
  "15-19" = "15-19",
  "20-29" = c("20-24", "25-29"),
  "30-39" = c("30-34", "35-39"),
  "40-49" = c("40-44", "45-49"),
  "50-59" = c("50-54", "55-59"),
  "60-64" = "60-64",
  "65-74" = c("65-69", "70-74"),
  "75+"   = "75+")

## --- 1. Contact matrix ------------------------------------------------------

# Raw inputs are fetched (and cached) by codes/fetch_data.R.
# Sourcing it here makes scaffold idempotent: re-running this script "just works"
# whether you've fetched before or not.
source(file.path("codes", "fetch_data.R"))

reconnect_path <- file.path(input_dir, "base_matrix.csv")
cm_long <- read.csv(reconnect_path, header = TRUE, stringsAsFactors = FALSE)
stopifnot(nrow(cm_long) == 5 * 5 * length(reconnect_ages) * length(reconnect_ages))

n_imd_p <- length(unique(cm_long$Participant_IMD))
n_imd_c <- length(unique(cm_long$Contact_IMD))
stopifnot(n_imd_p == 5, n_imd_c == 5)
na_new  <- length(new_ages)   # 10
nimd    <- 5
ng_new  <- na_new * nimd      # 50

# Build a fast lookup: for each (participant_IMD, contact_IMD, p_age, c_age) -> mean rate
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
        # participant-side: MEAN (TODO population-weighted; uniform for now)
        sum_over_contacts_per_subp <- vapply(
          sub_p, function(pa) {
            sum(vapply(sub_c, function(ca)
              mean_lookup[[key(is_p, is_c, pa, ca)]],
              numeric(1)))
          }, numeric(1))
        cm50[(is_p - 1)*na_new + ia_p,
             (is_c - 1)*na_new + ia_c] <- mean(sum_over_contacts_per_subp)
      }
    }
  }
}

write.table(cm50, file.path(output_dir, "Mas50.csv"),
            sep = ",", row.names = FALSE, col.names = FALSE)
cat(sprintf("Wrote %s: %d x %d (from Reconnect base_matrix.csv)\n",
            "data/Mas50.csv", nrow(cm50), ncol(cm50)))

## --- 2. Demographics (from ONS LSOA SYA + IoD 2025) -------------------------
## Build a (IMD quintile x model age band) population table by joining the
## ONS LSOA-level single-year-of-age data to the IoD 2025 LSOA -> decile
## lookup, then aggregating.
##   Decile -> Quintile: quintile = ceiling(decile / 2)
##   (deciles 1-2 -> quintile 1 most deprived, ..., 9-10 -> quintile 5 least)

suppressPackageStartupMessages(library(readxl))

ons_path <- file.path(input_dir, "ons_lsoa_syoa_2022-2024.xlsx")
iod_path <- file.path(input_dir, "iod2025_lsoa_ranks_deciles.csv")

# ONS: row 4 of the Mid-2024 sheet is the header (rows 1-3 are notes / merged title)
ons <- read_excel(ons_path, sheet = "Mid-2024 LSOA 2021", skip = 3)
# Columns: LAD 2023 Code, LAD 2023 Name, LSOA 2021 Code, LSOA 2021 Name, Total,
# then F0, F1, ..., F90 and M0, M1, ..., M90 (where F90/M90 are 90 and over).

# Sum female + male for each single year of age
sya_ages <- 0:90
pop_sya  <- sapply(sya_ages, function(a)
  as.numeric(ons[[paste0("F", a)]]) + as.numeric(ons[[paste0("M", a)]]))
colnames(pop_sya) <- sya_ages

# Aggregate single-year ages into the 10 model bands
band_def <- list(
  "0 to 4"   = 0:4,
  "5 to 14"  = 5:14,
  "15 to 19" = 15:19,
  "20 to 29" = 20:29,
  "30 to 39" = 30:39,
  "40 to 49" = 40:49,
  "50 to 59" = 50:59,
  "60 to 64" = 60:64,
  "65 to 74" = 65:74,
  "75+"      = 75:90)
band_names <- names(band_def)
pop_band <- sapply(band_def, function(ages)
  rowSums(pop_sya[, as.character(ages), drop = FALSE]))

# LSOA-level pop-by-band
lsoa_demog <- data.frame(
  lsoa = ons[["LSOA 2021 Code"]],
  pop_band,
  check.names = FALSE,
  stringsAsFactors = FALSE)

# IoD 2025: LSOA -> IMD decile
iod <- read.csv(iod_path, header = TRUE, check.names = FALSE,
                stringsAsFactors = FALSE)
decile_col <- grep("^Index of Multiple Deprivation \\(IMD\\) Decile",
                   names(iod), value = TRUE)[1]
lsoa_imd <- data.frame(
  lsoa     = iod[["LSOA code (2021)"]],
  quintile = ceiling(iod[[decile_col]] / 2),  # decile 1-2 -> quintile 1
  stringsAsFactors = FALSE)

# Inner join: keep only LSOAs present in both (England only, since IoD is England)
joined <- merge(lsoa_demog, lsoa_imd, by = "lsoa")
cat(sprintf("Joined %d LSOAs (England) on (LSOA 2021 code)\n", nrow(joined)))
rm(ons, pop_sya, pop_band, lsoa_demog, iod, lsoa_imd)

# Aggregate by (IMD quintile, age band)
demog10 <- do.call(rbind, lapply(1:5, function(q) {
  pops <- colSums(joined[joined$quintile == q, band_names, drop = FALSE])
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

write.csv(demog10, file.path(output_dir, "demographics_10age.csv"),
          row.names = FALSE)
cat(sprintf("Wrote %s: %d rows (expected %d)\n",
            "data/demographics_10age.csv",
            nrow(demog10), 5 * length(band_names)))
