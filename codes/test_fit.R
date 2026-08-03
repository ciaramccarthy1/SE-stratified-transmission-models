################################################################################
# Eyeball-fit beta to David's (rsvie) weekly total symptomatic cases.
#
# Target : data/no_vacc_weekly_by_age_outcome.csv  (David's "no vaccination" run)
#          -> filter outcome == "symptomatic", sum over age -> weekly total.
# Model  : our SEIRD (no vaccination, weekly incidence), run over a full 52-week
#          season, with beta swept over a few candidate values (beta_override
#          hook in modelrun.r bypasses the R0-derived beta).
#
# Just eyeball which beta curve best matches the target, then set it in parsR_.r.
################################################################################

library(here)
library(Rcpp)
library(magrittr)
library(ggplot2)
library(tidyverse)
library(rsvie)

RSVempty@post |> select(l1, l2, pA1, pA2, pA3, pA4, b1, phi, psi) |>
  apply(2, median)

## --- Target: David's weekly total symptomatic cases (across ages) -----------
no_vac_weekly <- read.csv(here("data", "no_vacc_weekly_by_age_outcome.csv"))
target <- no_vac_weekly |>
  filter(outcome == "symptomatic") |>
  group_by(week_no) |>
  summarise(cases_total = sum(cases), .groups = "drop")  #|>
 # mutate(week_no = week_no - 104)

cat(sprintf("Target (David, symptomatic): peak %.0f at week %d, season total %.0f\n",
            max(target$cases_total),
            target$week_no[which.max(target$cases_total)],
            sum(target$cases_total)))

## --- pset: RSV, no vaccination, weekly incidence ----------------------------
source("codes/setup.r")
pset$Disease        <- "RSV-illness"
pset$Vaccination    <- 0
pset$DailyIncidence <- 1            # daily incidence -> SEIRDasday_.cpp (aggregated to weeks below)
pset$Incidence      <- "Daily"
pset$Namevacc       <- ""
pset$TODAY          <- ""
pset$FIGURES        <- 0
pset$DIAGNOSTIC     <- 0
pset$SUMMARY        <- 0
pset$COMPILE        <- 1

# demographics follow setup.r's DavidDemog flag (single source of truth)
davidDemog <- isTRUE(pset$DavidDemog)

## --- Run the model over a full 52-week season for each candidate beta -------
betas <- c(0.098)

weeks_total <- 520
Tmax        <- weeks_total * 7      # 364 days
dt_fit      <- 0.1                  # must match parsR_.r dt

# NB: source() evaluates modelrun.r in the global env, so scenario_overrides
# must be a global (top-level for-loop), not local to a function.
model_curves <- list()
for (b in betas) {
  scenario_overrides <- list(
    beta_override = b,
    times = seq(0, Tmax),
    nt    = Tmax / dt_fit + 1,
    nd    = Tmax + 1
  )
  source("codes/modelrun.r")       # builds inputs, runs daily model -> `mas`
  # Daily incident symptomatic -> aggregate to weeks (week = ceiling(day/7)).
  daily <- as.numeric(mas$byw$Iw)                # length nd = Tmax+1 (days 1..365)
  wk    <- ceiling(seq_along(daily) / 7)         # day d -> week ceiling(d/7)
  weekly <- tapply(daily, wk, sum)
  weekly <- weekly[as.integer(names(weekly)) <= weeks_total]   # keep weeks 1-52
  model_curves[[sprintf("%.3f", b)]] <- data.frame(
    week_no = as.integer(names(weekly)),
    cases   = as.numeric(weekly),
    beta    = sprintf("beta = %.3f", b))
}
rm(scenario_overrides)
model_df <- do.call(rbind, model_curves)

cat(sprintf("\nModel population (Npop): %.0f\n", Npop))
model_df |>
  group_by(beta) |>
  summarise(peak = max(cases), peak_week = week_no[which.max(cases)],
            season_total = sum(cases), .groups = "drop") |>
  print()

## --- Overlay plot -----------------------------------------------------------
p2 <- ggplot() +
  geom_line(data = model_df,
            aes(week_no, cases, colour = beta), linewidth = 0.8) +
  #geom_line(data = target,
  #          aes(week_no, cases_total), colour = "black", linewidth = 1.1) +
  geom_point(data = target,
             aes(week_no, cases_total), colour = "black", size = 1) +
  labs(x = "Week", y = "Weekly total symptomatic cases",
       colour = NULL,
       title = sprintf("Total symptomatic: model (beta sweep) vs David, black (%s demographics)",
                       if (davidDemog) "David's stationary" else "real ONS")) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

p2

tag <- if (davidDemog) "daviddemog" else "ons"
out <- sprintf("output/test_fit_symptomatic_%s.pdf", tag)
ggsave(out, p2, width = 8, height = 5)
ggsave(sub("\\.pdf$", ".png", out), p2, width = 8, height = 5, dpi = 150)
cat(sprintf("\nWrote %s (and .png)\n", out))
