################################################################################
# Build (age x IMD quintile x year) demographic rate schedules for the ONS
# (MatchDavid = FALSE) 2026-2036 cohort-component run, from ONS 2024-based
# national population projections plus deprivation gradients/splits.
#
# Rates are TIME-VARYING (year-by-year, following the projection).
#
# Inputs (fetched/extracted by codes/fetch_data.R):
#   data/en_ppp_machine_readable.xlsx  : ONS 2024-based PRINCIPAL projection,
#       England, national, single year of age x sex x year. Sheets used:
#         Population            - mid-year population by age x sex x year
#         Mortality_assumptions - mortality rate per 100,000 by age x sex x year
#         Deaths                - death counts (used only to verify rate units)
#         Total_migration       - In/Out/Net migration by age x sex x year
#         Births                - births by mother's age x baby sex x year
#   [Bit B] ONS 2009-2020 age x IMD-decile mortality  -> mortality gradient
#   [Bit C] births by deprivation                     -> birth IMD split
#
# Outputs (data/):
#   [A] nat_mortality_9band.csv    national annual mortality rate, 9-band x year
#   [A] nat_netmigration_9band.csv national net migration (persons/yr), 9-band x yr
#   [A] nat_births_total.csv       national total births per year
#   [B/C] to follow (age x IMD x year schedules)
#
# Model age bands: 0-4, 5-14, 15-24, 25-34, 35-44, 45-54, 55-64, 65-74, 75+
#   (75+ = all ages >= 75). Sex-combined (model is not sex-stratified),
#   population-weighted for rates. Rates kept as annual per-person values here;
#   conversion to a daily hazard happens where the cpp is wired.
################################################################################

suppressPackageStartupMessages({library(here); library(readxl); library(dplyr); library(tidyr)})

ppp <- file.path(here("data"), "en_ppp_machine_readable.xlsx")

## --- A. National schedules from the ONS principal projection ----------------

ystart <- 2024:2036                                   # projection years to extract
rng    <- paste(ystart, ystart + 1, sep = " - ")      # rate/flow cols: "2024 - 2025"
popc   <- as.character(ystart)                        # population cols: "2024"

# assign a single year of age to one of our 9 model bands (75+ = all >= 75)
band_of <- function(a) {as.integer(cut(a, c(0,5,15,25,35,45,55,65,75,Inf), right = FALSE, labels = 1:9))}

# read a sheet, coerce Age -> integer (drops "Birth"/"All ages"/text to NA)
read_sheet <- function(sheet) {
  d <- suppressMessages(read_excel(ppp, sheet = sheet))
  d$a <- suppressWarnings(as.integer(d$Age))
  d
}
# long-form the year columns of a sheet, keeping single-year ages only
longify <- function(d, cols, value) {
  d |> select(Sex, a, any_of(cols)) |> filter(!is.na(a)) |>
    pivot_longer(any_of(cols), names_to = "k", values_to = value)
}

pop <- longify(read_sheet("Population"), popc, "pop") |> mutate(year = as.integer(k))
mort <- longify(read_sheet("Mortality_assumptions"), rng, "rate") |>
  mutate(year = as.integer(substr(k, 1, 4)))
deaths <- longify(read_sheet("Deaths"), rng, "deaths") |>
  mutate(year = as.integer(substr(k, 1, 4)))
migrate <- read_sheet("Total_migration") |> select(Sex, a, Flow, any_of(rng)) |>
  filter(!is.na(a)) |> pivot_longer(any_of(rng), names_to = "k", values_to = "n") |>
  mutate(year = as.integer(substr(k, 1, 4)))
births <- longify(read_sheet("Births"), rng, "births") |>
  mutate(year = as.integer(substr(k, 1, 4)))

# Unit check: the mortality "rate" should equal deaths/population * 100,000.
check <- mort |> inner_join(deaths, by = c("Sex","a","year", "k")) |>
  inner_join(pop, by = c("Sex","a","year")) |> filter(year == 2026) |>
  mutate(ratio = rate / (deaths / pop))
stopifnot(abs(median(check$ratio, na.rm = TRUE) - 1e5) < 100)   # confirms per-100,000
cat(sprintf("Mortality rates confirmed per 100,000 (median rate/(deaths/pop) = %.0f)\n",
            median(check$ratio, na.rm = TRUE)))

# Mortality -> annual per-person rate by band (sex-combined, population-weighted)
nat_mortality <- mort |> inner_join(pop, by = c("Sex","a","year")) |>
  mutate(band = band_of(a)) |> filter(!is.na(band)) |>
  group_by(band, year) |>
  summarise(mort_rate = sum(rate / 1e5 * pop) / sum(pop), .groups = "drop")

# Net migration (persons/yr) by band: the "Net" flow rows, summed over sex & age
nat_netmigrate <- migrate |> filter(Flow == "Net") |> mutate(band = band_of(a)) |>
  filter(!is.na(band)) |> group_by(band, year) |>
  summarise(net = sum(n), .groups = "drop")

# Total births per year (sum over mother's age and baby sex)
nat_births <- births |> group_by(year) |> summarise(births = sum(births), .groups = "drop")

write.csv(nat_mortality, file.path(here("data"), "nat_mortality_9band.csv"), row.names = FALSE)
write.csv(nat_netmigrate, file.path(here("data"), "nat_netmigration_9band.csv"), row.names = FALSE)
write.csv(nat_births,    file.path(here("data"), "nat_births_total.csv"), row.names = FALSE)
cat(sprintf("Wrote national schedules for %d-%d: mortality, net migration, births\n",
            min(ystart), max(ystart)))

## --- B. IMD mortality gradient -> mortality[band, IMD quintile, year] --------
## Split the national mortality schedule (A) across IMD quintiles using the relative
## deprivation gradient from ONS 2009-2020 age x decile mortality (Table 1, per 100k).
## Choices:
##  - gradient years: 2016-2020 (IMD 2019-based, closest to current deprivation);
##  - aggregate 5yr age groups -> our 9 bands and combine sexes with national
##    population weights (by sex x 5yr group); average the 2 deciles in a quintile
##    and the gradient years with equal weight;
##  - normalise the gradient so its population-weighted mean over quintiles = 1,
##    which preserves the national projected total after the split.

dep_path   <- file.path(here("data"), "dep_asmr_2009_2020.xlsx")
grad_years <- 2016:2020

# dep_mort 5yr age group -> lower edge (single year) -> our band
grp_lo <- c("<1"=0,"01-04"=1,"05-09"=5,"10-14"=10,"15-19"=15,"20-24"=20,"25-29"=25,
            "30-34"=30,"35-39"=35,"40-44"=40,"45-49"=45,"50-54"=50,"55-59"=55,
            "60-64"=60,"65-69"=65,"70-74"=70,"75-79"=75,"80-84"=80,"85-89"=85,"90+"=90)

dm <- suppressMessages(read_excel(dep_path, sheet = "Table 1", skip = 1))
names(dm) <- c("year","decile","sex","agegrp","rate")
dm <- dm |> mutate(year = as.integer(year), decile = as.integer(decile),
                   rate = as.numeric(rate), sex = ifelse(grepl("^M", sex), "M", "F"),
                   quint = ceiling(decile / 2), band = band_of(grp_lo[agegrp])) |>
  filter(year %in% grad_years, !is.na(rate), !is.na(band))

# national population weights by (sex, 5yr group), representative year 2024
grp_of <- function(a) names(grp_lo)[findInterval(a, grp_lo)]
popw <- pop |> filter(year == 2024) |>
  mutate(sex = ifelse(grepl("^M", Sex), "M", "F"), agegrp = grp_of(a)) |>
  group_by(sex, agegrp) |> summarise(w = sum(pop), .groups = "drop")

# rate[band, quintile]: pop-weighted over (sex x 5yr group), equal over decile & year
rate_agequint <- dm |> left_join(popw, by = c("sex","agegrp")) |>
  group_by(band, quint) |> summarise(rate = sum(rate * w) / sum(w), .groups = "drop")

# quintile population by band (from demographics_9age.csv) to normalise the gradient
ages9 <- c("0 to 4","5 to 14","15 to 24","25 to 34","35 to 44",
           "45 to 54","55 to 64","65 to 74","75+")
dg   <- read.csv(file.path(here("data"), "demographics_9age.csv"))
popq <- data.frame(band = match(dg$Age, ages9), quint = dg$IMD, popq = dg$Population)

grad <- rate_agequint |> left_join(popq, by = c("band","quint")) |> group_by(band) |>
  mutate(gradient = rate / (sum(rate * popq) / sum(popq))) |> ungroup()

# apply gradient to the national projected mortality schedule -> imd x band x year
imd_mortality <- nat_mortality |>
  left_join(grad |> select(band, quint, gradient), by = "band",
            relationship = "many-to-many") |>   # band-year x band-quintile -> band-imd-year
  transmute(band, imd = quint, year, mort_rate = mort_rate * gradient)
write.csv(imd_mortality, file.path(here("data"), "imd_mortality_9band.csv"), row.names = FALSE)

# sanity: gradient Q1(most deprived) vs Q5, and the normalisation
cat("\nMortality gradient (Q1 most-deprived / Q5 least-deprived) by band:\n")
gchk <- grad |> select(band, quint, gradient) |> tidyr::pivot_wider(names_from = quint,
          values_from = gradient, names_prefix = "Q")
print(as.data.frame(gchk |> mutate(across(-band, ~round(., 3)))))
nrm <- grad |> group_by(band) |> summarise(wmean = sum(gradient * popq) / sum(popq))
cat("normalisation check (pop-weighted mean of gradient per band, should all be 1):",
    paste(round(nrm$wmean, 3), collapse = " "), "\n")

