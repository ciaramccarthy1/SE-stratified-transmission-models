################################################################################
# Produces 10-age-band model inputs aligned with Reconnect survey
# imd_matrices contact matrices (https://github.com/lucy-gf/imd_matrices).
#
# New 10-band model layout (cuts aligned with 5-year bands so no
# fractional within-band splitting is needed):
#   0-4, 5-14, 15-19, 20-29, 30-39, 40-49, 50-59, 60-64, 65-74, 75+
#
# Inputs:
#   - data/base_matrix.csv         : Reconnect long-form 5-IMD x 16-age matrix
#                                    (5 x 5 IMD pairs x 16 x 16 age cells = 6400 rows)
#   - data/demographics2021.csv    : original 9-band ONS demographics
#
# Outputs:
#   - data/Mas50_urban.csv         : 50x50 wide contact matrix
#                                    (rows = participant ig = is*na + ia)
#   - data/demographics2021_10age.csv : 10-band demographics
#                                    (TODO: still scaffolded from 9-band source;
#                                    replace with real ONS 10-band data)
#
# Aggregation notes for the contact matrix:
#   Participant-side merge of {a1, a2} -> A is a population-weighted average
#     (one participant is in exactly one sub-band, so we average rates).
#   Contact-side merge of {b1, b2} -> B is a SUM
#     (a participant in A makes contacts with both sub-bands; counts add).
#   TODO: currently using UNIFORM weights on the participant side
#     (i.e. simple mean over sub-bands). Replace with population-weighted
#     means using real ONS sub-band populations when available.
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

# Source: IMD matrices repo (Reconnect survey, balanced).
# Downloaded on demand; cached locally so subsequent runs don't re-fetch.
reconnect_url  <- "https://raw.githubusercontent.com/lucy-gf/imd_matrices/main/matrices/base_matrix.csv"
reconnect_path <- file.path(input_dir, "base_matrix.csv")
if (!file.exists(reconnect_path)) {
  message("Downloading base_matrix.csv from ", reconnect_url)
  download.file(reconnect_url, reconnect_path, mode = "wb", quiet = TRUE)
}
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

write.table(cm50, file.path(output_dir, "Mas50_urban.csv"),
            sep = ",", row.names = FALSE, col.names = FALSE)
cat(sprintf("Wrote %s: %d x %d (from Reconnect base_matrix.csv)\n",
            "data/Mas50_urban.csv", nrow(cm50), ncol(cm50)))

## --- 2. Demographics --------------------------------------------------------
# The original demographics2021.csv has 9 bands with different cut points,
# so re-cutting to the new 10 bands needs splitting/merging assumptions.
# Old bands: 0-4, 5-11, 12-17, 18-29, 30-39, 40-49, 50-59, 60-69, 70+
# New bands: 0-4, 5-14, 15-19, 20-29, 30-39, 40-49, 50-59, 60-64, 65-74, 75+
#
# Linear-uniform splits (TODO: replace with real ONS 10-band data):
#   new "5-14"  = old "5-11" (7 years) + 3/6 of old "12-17"
#   new "15-19" = 3/6 of old "12-17"   + 2/12 of old "18-29"
#   new "20-29" = 10/12 of old "18-29"
#   new "60-64" = 5/10 of old "60-69"
#   new "65-74" = 5/10 of old "60-69"  + 5/20 of old "70+"
#   new "75+"   = 15/20 of old "70+"
#
# Where 70+ is assumed to span 70-89 effectively (width 20).

old_ages <- c("0 to 4","5 to 11","12 to 17","18 to 29","30 to 39",
              "40 to 49","50 to 59","60 to 69","70+")
new_age_labels <- c("0 to 4","5 to 14","15 to 19","20 to 29","30 to 39",
                    "40 to 49","50 to 59","60 to 64","65 to 74","75+")

demog <- read.csv(file.path(input_dir, "demographics2021.csv"),
                  header = TRUE, stringsAsFactors = FALSE)
blocks <- split(demog, list(demog$IMD, demog$rural), drop = TRUE)

expand_block <- function(b) {
  b <- b[match(old_ages, b$Age), ]
  pop <- setNames(b$Population, b$Age)
  new_pop <- c(
    pop["0 to 4"],
    pop["5 to 11"]   + 0.5 * pop["12 to 17"],     # new 5-14
    0.5 * pop["12 to 17"] + (2/12) * pop["18 to 29"],   # new 15-19
    (10/12) * pop["18 to 29"],                    # new 20-29
    pop["30 to 39"],
    pop["40 to 49"],
    pop["50 to 59"],
    0.5 * pop["60 to 69"],                        # new 60-64
    0.5 * pop["60 to 69"] + (5/20) * pop["70+"],  # new 65-74
    (15/20) * pop["70+"]                          # new 75+
  )
  data.frame(
    Age        = new_age_labels,
    IMD        = b$IMD[1],
    rural      = b$rural[1],
    Population = as.numeric(new_pop),
    tot_pop    = sum(new_pop),
    Proportion = as.numeric(new_pop) / sum(new_pop),
    stringsAsFactors = FALSE)
}

demog10 <- do.call(rbind, lapply(blocks, expand_block))
demog10 <- demog10[order(demog10$IMD,
                         factor(demog10$rural, levels = c("Urban","Rural")),
                         match(demog10$Age, new_age_labels)), ]
rownames(demog10) <- NULL

write.csv(demog10, file.path(output_dir, "demographics2021_10age.csv"),
          row.names = FALSE)
cat(sprintf("Wrote %s: %d rows (expected %d)\n",
            "data/demographics2021_10age.csv",
            nrow(demog10),
            length(new_age_labels) * length(unique(demog$IMD)) *
              length(unique(demog$rural))))
