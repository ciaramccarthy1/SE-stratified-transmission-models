################################################################################
# Download raw input data files for the SE-stratified transmission model.
#
# Files downloaded by this script (see data/SOURCES.md for full provenance):
#   data/base_matrix.csv : age x IMD contact matrix (Reconnect survey,
#                          16 five-year bands, 5 IMD quintiles, balanced)
#                          Source: lucy-gf/imd_matrices repo
#
# Files NOT auto-fetched (require manual sourcing - see data/SOURCES.md):
#   data/demographics2021.csv  - ONS population, inherited from upstream;
#                                rebuild method TODO
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

# fetch_if_missing(
#   url   = "TODO ONS: 2011 rural-urban classification by LSOA (still current)",
#   path  = file.path(input_dir, "lsoa_ruc.csv"),
#   label = "LSOA rural-urban classification")
##
## After downloading, a separate prep step (NOT this script) should join +
## aggregate them into data/demographics2021_real.csv with the 10 model bands.

## --- Source-paper per-age parameter tables (RSV) ---------------------------
## Feeds: u, y, m, h, mH, rrep per age in parsR_.r / parsRv_.r.

## RSV - Hodgson 2020 (Lancet ID); model code may live on GitHub
# fetch_if_missing(
#   url   = "TODO Hodgson 2020 supp / GitHub: per-age u, y, IFR, hospitalisation",
#   path  = file.path(input_dir, "hodgson2020_rsv.csv"),
#   label = "Hodgson 2020 RSV per-age params")

## --- RSV vaccine efficacy from Phase 3 trials ------------------------------
## Feeds: VE_inf, VE_sym, VE_hosp, VE_mort in parsRv_.r.
## Trial papers report composite endpoints (typically VE vs MA-LRTI and
## VE vs severe LRTD); decomposition into the four model VEs requires a
## modelling choice. Useful to record raw trial-reported values for reference.
# fetch_if_missing(
#   url   = "TODO GSK Arexvy (RSVPreF3) Phase 3 VE estimates by endpoint+season",
#   path  = file.path(input_dir, "arexvy_ve.csv"),
#   label = "GSK Arexvy VE")
#
# fetch_if_missing(
#   url   = "TODO Pfizer Abrysvo Phase 3 VE estimates by endpoint+season",
#   path  = file.path(input_dir, "abrysvo_ve.csv"),
#   label = "Pfizer Abrysvo VE")
#
# fetch_if_missing(
#   url   = "TODO Moderna mRESVIA Phase 3 VE estimates by endpoint+season",
#   path  = file.path(input_dir, "mresvia_ve.csv"),
#   label = "Moderna mRESVIA VE")


## --- UK vaccination programme uptake data ----------------------------------
## Feeds: vcov in pars*v_.r files.
## UKHSA / NHS England publish uptake periodically (Green Book + dashboards).
# fetch_if_missing(
#   url   = "TODO UKHSA/NHS England: RSV 75+ programme uptake (when reported)",
#   path  = file.path(input_dir, "rsv_75plus_uptake.csv"),
#   label = "RSV 75+ uptake by IMD/region")
#
# fetch_if_missing(
#   url   = "TODO UKHSA annual flu vaccine uptake (clinical risk groups + 65+)",
#   path  = file.path(input_dir, "flu_uptake.csv"),
#   label = "Flu uptake by age")
#
# fetch_if_missing(
#   url   = "TODO UKHSA COVID booster uptake by age + IMD",
#   path  = file.path(input_dir, "covid_uptake.csv"),
#   label = "COVID uptake by age")


## --- Hospital length-of-stay data ------------------------------------------
## Feeds: rH per disease in pars*_.r files.
## HES (Hospital Episode Statistics) is restricted; published reports give
## aggregate LOS by ICD-10 code, age, and sometimes IMD.
# fetch_if_missing(
#   url   = "TODO NHS Digital HES respiratory admissions: LOS by age, ICD-10",
#   path  = file.path(input_dir, "hes_resp_los.csv"),
#   label = "HES respiratory LOS by age")


## --- Sub-band populations for population-weighted contact aggregation ------
## Feeds: scaffold_10age_data.R participant-side aggregation (currently
## uniform mean - TODO population-weighted).
## If ONS demographics above are added at 5-year-band granularity, this
## becomes a derived calculation rather than a separate fetch.
# (no new fetch; uses ons_pop_lsoa.csv aggregated by 5-year bands)


message("fetch_data.R: all raw inputs present and validated.")
