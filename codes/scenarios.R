################################################################################
# Run RSV 75+ vaccination scenarios and compare hospitalisation outcomes
# by IMD quintile.
#
# Each scenario defines a `vcov` vector (length na*nimd, IMD-major) which
# replaces parsRv_.r's default before the cpp model runs.
################################################################################

suppressPackageStartupMessages({
  library(here)
  library(Rcpp)
  library(magrittr)
  library(ggplot2)
  library(tidyverse)
  library(MetBrewer)
})

setwd(here())

## --- pset for all scenarios -------------------------------------------------
source("codes/setup.r")
pset$Disease        <- "RSV-illness"
pset$Vaccination    <- 1
pset$DailyIncidence <- 1
pset$Incidence      <- "Daily"
pset$Namevacc       <- "vaccine_"
pset$TODAY          <- ""
pset$FIGURES        <- 0
pset$DIAGNOSTIC     <- 0
pset$SUMMARY        <- 0
pset$COMPILE        <- 1


## --- Build the comparator scenarios ----------------------------------------
na   <- 9
nimd <- 5

build_vcov <- function(uptake_per_quintile, na, nimd) {
  stopifnot(length(uptake_per_quintile) == nimd)
  vc <- matrix(0, na, nimd)
  vc[na, ] <- uptake_per_quintile        # band 9 = last band (75+) only
  as.vector(vc)                          # IMD-major
}

# Status quo: UKHSA uptake by IMD quintile 
uptake_status_quo <- read.csv("data/rsv_uptake_by_imd_quintile.csv")
uptake_status_quo <- uptake_status_quo$uptake_pct[
  order(uptake_status_quo$quintile)] / 100

# Population-weighted mean of current uptake = dose-neutral equal-coverage counterfactual
demog        <- read.csv("data/demographics_9age.csv")
pop_75_by_q  <- demog$Population[demog$Age == "75+"]
stopifnot(length(pop_75_by_q) == nimd)
equal_uptake <- sum(uptake_status_quo * pop_75_by_q) / sum(pop_75_by_q)

scenarios <- list(
  status_quo = build_vcov(uptake_status_quo,         na, nimd),
  equal_avg  = build_vcov(rep(equal_uptake, nimd),   na, nimd))

cat("Status-quo uptake by quintile:",
    paste0(round(uptake_status_quo * 100, 1), "%", collapse = ", "), "\n")
cat("Equal uptake (pop-weighted mean):",
    paste0(round(equal_uptake * 100, 1), "%"), "\n\n")

## --- Run each scenario ------------------------------------------------------
results <- list()
for (nm in names(scenarios)) {
  cat(sprintf("=== Running scenario: %s ===\n", nm))
  scenario_overrides <- list(vcov = scenarios[[nm]])
  source("codes/modelrun.r")
  results[[nm]] <- list(
    mas = mas,
    Npop = Npop, Na = Na, Ns = Ns,
    parsum = parsum)
  rm(mas)
}
rm(scenario_overrides)


## --- Build comparison data frame --------------------------------------------
# Hw_ag  = cpp output: daily H by (age x IMD) [nd x ng], unvacc chain
# Hwv_ag = same, vacc-breakthrough chain
# Slice to 75+ band (row na of the reshaped (age x IMD) matrix) per quintile.
df <- do.call(rbind, lapply(names(results), function(nm) {
  r          <- results[[nm]]
  totH_g     <- colSums(r$mas$byaw$Hw_ag) +
                colSums(r$mas$byaw$Hwv_ag)        # length ng, IMD-major
  totH_mat   <- matrix(totH_g, na, nimd)           # rows=age, cols=IMD
  cumH_75pl  <- totH_mat[na, ]                     # length nimd: 75+ per quintile
  data.frame(
    scenario   = nm,
    quintile   = factor(1:nimd,
                        labels = c("1 (most dep.)", "2", "3", "4",
                                   "5 (least dep.)")),
    cumH_75pl = cumH_75pl,
    cumH_per_100k = cumH_75pl / pop_75_by_q * 1e5,
    stringsAsFactors = FALSE)
}))


## --- Plot -------------------------------------------------------------------
# Total hospitalisations
p_tot <- ggplot(df, aes(x = quintile, y = cumH_75pl, fill = scenario)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  scale_fill_manual(
    values = setNames(met.brewer("Derain", n = length(scenarios)),
                      names(scenarios)),
    labels = c(status_quo = "Status quo uptake",
               equal_avg  = "Equal coverage (population-weighted mean)")) +
  labs(x = "IMD quintile",
       y = "Total RSV hospitalisations in 75+ age group",
       fill = NULL,
       title = "RSV hospitalisations in 75+ age group by IMD quintile") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

p_100k <- ggplot(df, aes(x = quintile, y = cumH_per_100k, fill = scenario)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  scale_fill_manual(
    values = setNames(met.brewer("Derain", n = length(scenarios)),
                      names(scenarios)),
    labels = c(status_quo = "Status quo uptake",
               equal_avg  = "Equal coverage (population-weighted mean)")) +
  labs(x = "IMD quintile",
       y = "Total RSV hospitalisations in 75+\nper 100k population",
       fill = NULL,
       title = "RSV hospitalisations in 75+ age group by IMD quintile") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")
  
  out_pdf <- "output/scenarios_RSV_hosp_tot_IMD.pdf"
  out_png <- "output/scenarios_RSV_hosp_tot_IMD.png"
  ggsave(out_pdf, plot = p_tot, width = 7, height = 5)
  ggsave(out_png, plot = p_tot, width = 7, height = 5, dpi = 200)
  
  out_pdf <- "output/scenarios_RSV_hosp_100k_IMD.pdf"
  out_png <- "output/scenarios_RSV_hosp_100k_IMD.png"
  ggsave(out_pdf, plot = p_100k, width = 7, height = 5)
  ggsave(out_png, plot = p_100k, width = 7, height = 5, dpi = 200)
  
## --- Print summary numbers --------------------------------------------------
print_summary <- function(value_col, label) {
  cat(sprintf("\n%s, by IMD quintile:\n", label))
  print(df %>%
          pivot_wider(id_cols = quintile,
                      names_from = scenario,
                      values_from = {{ value_col }}) %>%
          mutate(diff = status_quo - equal_avg,
                 pct_change_vs_equal = 100 * diff / equal_avg))
}

print_summary(cumH_75pl,     "Total cumulative hospitalisations in 75+")
print_summary(cumH_per_100k, "Cumulative hospitalisations in 75+ per 100k 75+ pop")
