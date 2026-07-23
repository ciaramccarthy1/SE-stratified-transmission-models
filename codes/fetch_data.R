################################################################################
# Download raw input data files for the SE-stratified transmission model.
#
# Files downloaded by this script (see data/SOURCES.md):
#   data/base_matrix.csv                : age x IMD contact matrix
#                                         (Reconnect survey via lucy-gf/imd_matrices)
#   data/ons_lsoa_syoa_2022-2024.xlsx   : LSOA population by single year of age + sex
#                                         (ONS mid-2022 to mid-2024)
#   data/iod2025_lsoa_ranks_deciles.csv : LSOA -> IMD 2025 lookup
#                                         (GOV.UK English Indices of Deprivation 2025)

#
# codes/prepare_model_inputs.R: joins ons_lsoa_syoa with iod2025_lsoa_ranks on LSOA code, 
# aggregates to IMD quintile x model age bands, and writes data/demographics_9age.csv.
#
################################################################################

input_dir <- file.path("data")
if (!dir.exists(input_dir)) dir.create(input_dir, recursive = TRUE)

fetch_if_missing <- function(url, path, expected_bytes = NULL, label = path) {
  if (file.exists(path)) {
    message(sprintf("Using cached %s (delete to force re-download)", path))
    return(invisible(path))
  }
  message(sprintf("Downloading %s from %s", label, url))
  download.file(url, path, mode = "wb", quiet = TRUE)
  if (!is.null(expected_bytes)) {
    got <- file.info(path)$size
    if (got < expected_bytes) {
      stop(sprintf("Downloaded %s is only %d bytes (expected at least %d). Possible network truncation?",
                   path, got, expected_bytes))
    }
  }
  invisible(path)
}

## --- Reconnect contact matrix ----------------------------------------------
fetch_if_missing(
  url = "https://raw.githubusercontent.com/lucy-gf/imd_matrices/main/matrices/base_matrix.csv",
  path = file.path(input_dir, "base_matrix.csv"),
  label = "Reconnect base_matrix.csv")

# Sanity-check Reconnect matrix shape:
# 5 participant_IMD x 5 contact_IMD x 16 age_p x 16 age_c = 6400 rows + header
.cm <- read.csv(file.path(input_dir, "base_matrix.csv"),
                header = TRUE, stringsAsFactors = FALSE, nrows = 6500)
stopifnot(nrow(.cm) == 6400)
rm(.cm)

## Population by age and IMD ##

## --- ONS LSOA population by single year of age and sex --------------------
## Most recent ONS release: mid-2022 revised through mid-2024 (Nov 2025).
## Excel workbook with one sheet per year. 

fetch_if_missing(
  url   = paste0("https://www.ons.gov.uk/file?uri=",
                 "/peoplepopulationandcommunity/populationandmigration/",
                 "populationestimates/datasets/",
                 "lowersuperoutputareamidyearpopulationestimates/",
                 "mid2022revisednov2025tomid2024/sapelsoasyoa20222024.xlsx"),
  path  = file.path(input_dir, "ons_lsoa_syoa_2022-2024.xlsx"),
  expected_bytes = 1e6,
  label = "ONS LSOA mid-2022 to mid-2024 population by SYA and sex")


## --- LSOA -> IMD 2025 lookup (English Indices of Deprivation 2025) --------
## File 7 from gov.uk/government/statistics/english-indices-of-deprivation-2025
## Single CSV containing for every English LSOA: all 7 IMD domain ranks +
## scores + deciles, plus population denominators. Convert decile -> quintile
## downstream (quintile = ceiling(decile / 2)).
fetch_if_missing(
  url   = paste0("https://assets.publishing.service.gov.uk/media/",
                 "691ded56d140bbbaa59a2a7d/",
                 "File_7_IoD2025_All_Ranks_Scores_Deciles_",
                 "Population_Denominators.csv"),
  path  = file.path(input_dir, "iod2025_lsoa_ranks_deciles.csv"),
  expected_bytes = 1e6,
  label = "IoD 2025 LSOA ranks/scores/deciles + population denominators")

## --- Source-paper per-age parameter tables (RSV) ---------------------------
## Feeds: u, y, m, h, mH, rrep per age in parsR_.r / parsRv_.r.

## RSV - Hodgson 2020 (Lancet ID) - get from GitHub
# fetch_if_missing(
#   url   = "TODO Hodgson 2020 supp / GitHub: per-age u, y, IFR, hospitalisation",
#   path  = file.path(input_dir, "hodgson2020_rsv.csv"),
#   label = "Hodgson 2020 RSV per-age params")

## --- RSV vaccine efficacy from Phase 3 trials ------------------------------
## Feeds: VE_inf, VE_sym, VE_hosp, VE_mort in parsRv_.r.
## Trial papers report composite endpoints (typically VE vs MA-LRTI and
## VE vs severe LRTD)
##


## --- England RSV vaccination uptake by IMD decile ---------------------
## Source: UKHSA "RSV older adults vaccination coverage in England" report.
## Updated periodically; URL below points to the January 2026 report.
## TO DO: published table includes routine + catch-up cohorts. Routine cohort
## uptake is currently lower than the all-cohort figure - need to
## update with cohort-specific values when separable.
## Feeds: vcov in parsRv_.r (currently uniform 100% in band 10 / 75+ only).
rsv_uptake_url  <- paste0("https://www.gov.uk/government/statistics/",
                          "respiratory-syncytial-virus-rsv-older-adults-vaccination-coverage-in-england/",
                          "respiratory-syncytial-virus-rsv-older-adults-vaccination-coverage-in-england-january-2026-report")
rsv_uptake_html <- file.path(input_dir, "ukhsa_rsv_uptake_jan2026.html")
fetch_if_missing(
  url   = rsv_uptake_url,
  path  = rsv_uptake_html,
  expected_bytes = 1e5,
  label = "UKHSA RSV uptake (Jan 2026 report)")

rsv_uptake_csv <- file.path(input_dir, "rsv_uptake_by_imd_decile.csv")
if (!file.exists(rsv_uptake_csv)) {
  .rsv_html <- paste(readLines(rsv_uptake_html, warn = FALSE), collapse = "\n")
  # Extract all <table>...</table> blocks individually, then pick the one containing "Deprivation deciles".
  .all_tables <- regmatches(.rsv_html,
                            gregexpr("(?s)<table>.*?</table>", .rsv_html, perl = TRUE))[[1]]
  .tbl_block  <- .all_tables[grepl("Deprivation deciles", .all_tables, fixed = TRUE)]
  stopifnot(length(.tbl_block) == 1)
  # 10 rows x 2 cells = 20 <td> contents, alternating decile / uptake_pct
  .cells <- regmatches(.tbl_block,
                       gregexpr("<td[^>]*>([^<]+)</td>", .tbl_block, perl = TRUE))[[1]]
  .cells <- trimws(gsub("<[^>]+>", "", .cells))
  stopifnot(length(.cells) == 20)
  .rsv_uptake_df <- data.frame(
    decile     = as.integer(sub(" .*", "", .cells[seq(1, 20, 2)])),
    uptake_pct = as.numeric(.cells[seq(2, 20, 2)]))
  write.csv(.rsv_uptake_df, rsv_uptake_csv, row.names = FALSE)
  cat(sprintf("Parsed RSV uptake table from HTML: %d deciles -> %s\n",
              nrow(.rsv_uptake_df), rsv_uptake_csv))
  rm(.rsv_html, .tbl_block, .cells, .rsv_uptake_df)
} else {
  message(sprintf("Using cached %s (delete to re-parse from HTML)", rsv_uptake_csv))
}

## TO THINK ABOUT: Need uptake by quintile -> use population weighting once have popn/quintile (although they should be similar)

## --- RSV hospital length-of-stay -------------------------------------------
## Feeds: rH in parsR_.r / parsRv_.r (currently 1/6 placeholder).

## --- Population-weighted contact aggregation -------------------------------
## prepare_model_inputs.R currently aggregates contact rates participant-side
## with uniform mean (TO DO -> population-weighted). 

