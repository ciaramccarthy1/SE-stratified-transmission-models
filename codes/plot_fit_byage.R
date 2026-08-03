################################################################################
# Plot weekly symptomatic incidence by age band: model vs David (rsvie).
#
# David reports 25 age bands; our model has 9. We aggregate David's 25 -> our
# 9 by age-range overlap (uniform-within-band assumption for his 10-year adult
# bands; the 0-4 / 65-74 / 75+ overlaps are exact). Facet by our 9 bands.
#
# Run from the repo root:  Rscript codes/plot_fit_byage.R
################################################################################

suppressPackageStartupMessages({
  library(here); library(Rcpp); library(magrittr); library(ggplot2); library(tidyverse)
})
setwd(here())
source("codes/setup.r")
pset$Disease<-"RSV-illness"; pset$Vaccination<-0; pset$DailyIncidence<-1
pset$Incidence<-"Daily"; pset$Namevacc<-""; pset$TODAY<-""
pset$FIGURES<-0; pset$DIAGNOSTIC<-0; pset$SUMMARY<-0; pset$COMPILE<-1

# Demographics follow the DavidDemog flag in codes/setup.r 
# TRUE = David's stationary populations (reproduce his output), FALSE = real ONS.
# Change it in setup.r
davidDemog <- isTRUE(pset$DavidDemog)

# baseline transmission to plot (eyeball fit)
beta0 <- 0.098

## our 9 bands (aligned with David's)
o_lo <- c(0,5,15,25,35,45,55,65,75); o_hi <- c(5,15,25,35,45,55,65,75,90)
o_nm <- factor(c("0-4","5-14","15-24","25-34","35-44","45-54","55-64","65-74","75+"),
               levels=c("0-4","5-14","15-24","25-34","35-44","45-54","55-64","65-74","75+"))

## --- model: symptomatic incidence by age (Iw_a), full 52-week season ---
scenario_overrides <- list(beta_override=beta0, times=seq(0,364), nt=3641, nd=365)
invisible(capture.output(source("codes/modelrun.r")))
Iw_a <- mas$byaw$Iw_a                       # nd x 9 (daily symptomatic by age)
wk   <- ceiling(seq_len(nrow(Iw_a))/7)
mod_wk <- apply(Iw_a, 2, function(col) tapply(col, wk, sum))   # 53 x 9
mod_wk <- mod_wk[as.integer(rownames(mod_wk))<=52, ]
mod_df <- as.data.frame(mod_wk) |> setNames(as.character(o_nm)) |>
  mutate(week_no=1:52) |> pivot_longer(-week_no, names_to="band", values_to="cases") |>
  mutate(band=factor(band, levels=levels(o_nm)),
         src=sprintf("Model (beta0=%.3f)", beta0))

## --- David: symptomatic by 25 bands -> aggregate to our 9 by age overlap ---
d_lo <- c((0:11)/12, 1,2,3,4, 5,10, 15,25,35,45,55,65,75)
d_hi <- c((1:12)/12, 2,3,4,5, 10,15, 25,35,45,55,65,75,90)
W <- outer(1:25, seq_along(o_lo), Vectorize(function(b,j)
  max(0, min(o_hi[j],d_hi[b]) - max(o_lo[j],d_lo[b])) / (d_hi[b]-d_lo[b])))  # 25 x na
dd <- read.csv("data/no_vacc_weekly_by_age_outcome.csv") |>
  filter(outcome=="symptomatic") |>
  select(age_group, week_no, cases) |>
  pivot_wider(names_from=age_group, values_from=cases) |> arrange(week_no) # |>
 # mutate(week_no = week_no - 104)
Mmat <- as.matrix(dd[, as.character(1:25)])              # 52 x 25
dav  <- Mmat %*% W                                       # 52 x na
dav_df <- as.data.frame(dav) |> setNames(as.character(o_nm)) |>
  mutate(week_no=dd$week_no) |> pivot_longer(-week_no, names_to="band", values_to="cases") |>
  mutate(band=factor(band, levels=levels(o_nm)), src="David (aggregated)")

## --- plot ---
both <- bind_rows(mod_df, dav_df)
p <- ggplot(both, aes(week_no, cases, colour=src)) +
  geom_line(linewidth=0.8) +
  facet_wrap(~band, scales="free_y", ncol=2) +
  scale_colour_manual(values=setNames(c("black","#d1495b"),
                                       c("David (aggregated)", sprintf("Model (beta0=%.3f)", beta0)))) +
  labs(x="Week", y="Weekly symptomatic cases", colour=NULL,
       title=sprintf("Weekly symptomatic cases by age: model vs David (%s demographics)",
                     if (davidDemog) "David's stationary" else "real ONS")) +
  theme_minimal(base_size=10) + theme(legend.position="bottom")
p
outfile <- if (davidDemog) "output/test_fit_byage_daviddemog.png" else "output/test_fit_byage_current.png"
ggsave(outfile, p, width=9, height=11, dpi=130)
cat("Wrote", outfile, "\n")
