## ============================================================
## extra_hazard_test.R
## Supplementary analysis: hazards vs probabilities and
## counterperturbation patterns (CR and CTMC)
##
## This script:
##  1. Reads the Lievre annual transition probabilities.
##  2. Shows why logm()-based CTMC hazards are problematic.
##  3. Constructs competing-risks (CR) hazards from the same probs
##     and shows they are demographically plausible.
##  4. Compares original probs vs probs implied by CR hazards.
##  5. Runs perturbation experiments:
##       - CR mapping hazards -> probs
##       - CTMC mapping hazards -> probs via expm(Q)
##     and computes fractions of counterperturbation
##     absorbed by self-transitions and within the origin state.
##  6. (Optional) Adds a Taylor-based approximation block.
##  7. (Template) Shows how one might compare DemoDecomp::horiuchi()
##     on CR hazards vs CR-derived probabilities.
##
## These analyses are *illustrative* and based on hazards that
## imply plausible probabilities, but they do **not** reproduce
## the exact probabilities used in the main manuscript.
## ============================================================

## ------------------------------------------------------------
## 0) Setup
## ------------------------------------------------------------

library(tidyverse)
library(expm)
library(compositions)
library(DemoDecomp)
source("R/00_functions_classic.R")
source("R/00_sensitivity_functions.R")
source("R/extra_hazards_compare_functions.R")

## ------------------------------------------------------------


# Note: in the following, you could reproduce this using probs_monthly
# instead of probs, in which case you'll need to always set dt = 1/12 
# instead of dt = 1.

dt = 1; dt = 1/12

if (dt == 1){
  probs <- readr::read_csv("transitions_lievre2003_annual.csv") |>
    pivot_longer(
      cols      = -c(age, sex),
      names_to  = "transition",
      values_to = "p"
    ) |>
    mutate(transition = paste0("p_", transition)) |>
    pivot_wider(
      names_from  = "transition",
      values_from = "p"
    )
} else {
  probs <- readr::read_csv("transitions_lievre2003_monthly.csv") |>
    pivot_longer(
      cols      = -c(age, sex),
      names_to  = "transition",
      values_to = "p"
    ) |>
    mutate(transition = paste0("p_", transition)) |>
    pivot_wider(
      names_from  = "transition",
      values_from = "p"
    )
}

## ------------------------------------------------------------
##    Compare CTMC-logm hazards vs CR hazards
##    and original probs vs CR-derived probs
## ------------------------------------------------------------

## 3.1 Hazards: logm-based CTMC vs CR
haz_ctmc_logm <- probs_to_hazards_ctmc_logm(probs, dt = dt)
haz_cr        <- probs_to_hazards_cr(probs, dt = dt)

## Example: look at h_HD and h_HU by sex and age
haz_long <- bind_rows(
  haz_ctmc_logm |>
    mutate(variant = "CTMC_logm"),
  haz_cr |>
    mutate(variant = "CR")
) |>
  pivot_longer(
    cols      = c(h_HU, h_HD, h_UH, h_UD),
    names_to  = "hazard",
    values_to = "h"
  )

## Plot: hazard comparison (e.g., h_HD)
haz_long |> 
  mutate(sex = if_else(sex == "f","females","males")) |> 
ggplot(
  aes(x = age, y = h, colour = hazard, linetype = variant)
) +
  geom_line() +
  facet_wrap(~ sex) +
  theme_minimal() +
  labs(
    title = "Hazards from H: CTMC logm vs CR",
    y     = "hazard (log scale)",
    x     = "Age"
  ) +
  scale_y_log10()

## 3.2 Probabilities: original vs CR-CTMC-implied probabilities
# these are paired to haz_cr
probs_cr <- hazards_to_probs_ctmc(haz_cr, dt = dt)
# Note: you can get back the original probabilities exatly
# using hazards_to_probs_cr(), but then we can't do this 
# comparison!

probs_orig_long <- probs |>
  pivot_longer(
    cols      = starts_with("p_"),
    names_to  = "transition",
    values_to = "p_orig"
  )

probs_cr_long <- probs_cr |>
  pivot_longer(
    cols      = starts_with("p_"),
    names_to  = "transition",
    values_to = "p_cr"
  )

probs_compare <- probs_orig_long |>
  inner_join(
    probs_cr_long,
    by = c("sex", "age", "transition")
  ) |>
  mutate(diff = p_cr - p_orig)

## Plot example: p_HD original vs CR
probs_compare |> 
  filter(!transition%in%c("p_HH","p_UU"))  |> 
  select(-diff) |> 
  pivot_longer(c(p_orig, p_cr), names_to = "variant",values_to = "p") |> 
  mutate(variant = if_else(variant == "p_cr","CR","Original")) |> 
  ggplot(aes(x = age,color=transition,y=p, linetype=variant)) +
    geom_line() +
    facet_wrap(~ sex) +
    theme_minimal() +
    labs(
      title = "Original vs CR-implied probabilities",
      y     = "probability",
      x     = "Age"
    ) +
  scale_y_log10()

## In what follows, we treat haz_cr as our working hazards and
## probs_cr as the corresponding probabilities. These are
## demographically more plausible than the logm-based hazards,
## but they do not exactly match the original probabilities
## used in the main manuscript.

# but to clear up potential confusion of what probs_cr are, here it is:
# 1. take original probs
# 2. convert them to hazards using cr (competing risks) transform
# 3. convert back to probabilities using ctmc transform (expm) these are probs_cr!
# 4. we can toggle between haz_cr and probs_cr using expm(), logm(), respectively,
# so we now check to see what happens to probs_cr when we perturb haz_cr. We could
# not do this with the original data because logm() of those probabilities gives
# junk. But these tests should give valid feedback on how much self-transitions
# absorb such perturbations in practice.


## ------------------------------------------------------------
## 6) Run perturbation experiments: CR and CTMC (based on CR hazards)
## ------------------------------------------------------------

## Choose perturbation size and dt
eps_pert <- 1e-3

## 6.1 CR mapping hazards -> probs
perturb_results_cr <- run_hazard_perturbation_cr(
  haz_df = haz_cr,
  dt     = dt,
  eps    = eps_pert
)

## 6.2 CTMC mapping hazards -> probs via expm(Q)
perturb_results_ctmc <- run_hazard_perturbation_ctmc(
  haz_df = haz_cr,
  dt     = dt,
  eps    = eps_pert
)

## ------------------------------------------------------------
## 7) Summaries: fractions absorbed by self-transitions
##    and within the origin state
## ------------------------------------------------------------

## 7.1 CR: fraction of shrink absorbed by self-transitions (global)
cr_frac_global <- compute_counterpert_fraction_staying_global(perturb_results_cr)

## 7.2 CR: fraction of shrink absorbed by self-transitions within origin row
cr_frac_origin <- compute_counterpert_fraction_staying_origin_only(perturb_results_cr)

## 7.3 CR: fraction of total shuffle happening in origin row
cr_frac_shuffle_origin <- compute_shuffle_fraction_origin(perturb_results_cr)

## Quick summaries
cr_frac_origin_summary <- cr_frac_origin |>
  group_by(hazard_perturbed) |>
  summarise(
    median_frac_stay = median(frac_stay_row, na.rm = TRUE),
    min_frac_stay    = min(frac_stay_row, na.rm = TRUE),
    max_frac_stay    = max(frac_stay_row, na.rm = TRUE),
    .groups = "drop"
  )
print(cr_frac_origin_summary)

## 7.4 CTMC: same quantities
ctmc_frac_global <- compute_counterpert_fraction_staying_global(perturb_results_ctmc)
ctmc_frac_origin <- compute_counterpert_fraction_staying_origin_only(perturb_results_ctmc)
ctmc_frac_shuffle_origin <- compute_shuffle_fraction_origin(perturb_results_ctmc)

ctmc_frac_global <- ctmc_frac_origin |>
  group_by(hazard_perturbed) |>
  summarise(
    median_frac_stay = median(frac_stay_row, na.rm = TRUE),
    min_frac_stay    = min(frac_stay_row, na.rm = TRUE),
    max_frac_stay    = max(frac_stay_row, na.rm = TRUE),
    q.025_frac_stay  = quantile(frac_stay_row, probs = 0.025, na.rm = TRUE),
    q.975_frac_stay  = quantile(frac_stay_row, probs = 0.975, na.rm = TRUE),
    .groups = "drop"
  )
ctmc_frac_global
# The vast majority of total counterperturbation is taken by the self-transition
# located in the same origin state as the hazard perturbation.



ctmc_frac_origin_summary <- ctmc_frac_origin |>
  group_by(hazard_perturbed) |>
  summarise(
    median_frac_stay = median(frac_stay_row, na.rm = TRUE),
    min_frac_stay    = min(frac_stay_row, na.rm = TRUE),
    max_frac_stay    = max(frac_stay_row, na.rm = TRUE),
    q.025_frac_stay  = quantile(frac_stay_row, probs = 0.025, na.rm = TRUE),
    q.975_frac_stay  = quantile(frac_stay_row, probs = 0.975, na.rm = TRUE),
    .groups = "drop"
  )
# vast majority of counterperturbation among origin-state perturbations happens
# in the self-transition
ctmc_frac_origin_summary

ctmc_frac_shuffle_origin_summary <- ctmc_frac_shuffle_origin |>
  group_by(hazard_perturbed) |>
  summarise(
    median_frac_stay = median(frac_shuffle_origin, na.rm = TRUE),
    min_frac_stay    = min(frac_shuffle_origin, na.rm = TRUE),
    max_frac_stay    = max(frac_shuffle_origin, na.rm = TRUE),
    q.025_frac_stay  = quantile(frac_shuffle_origin, probs = 0.025, na.rm = TRUE),
    q.975_frac_stay  = quantile(frac_shuffle_origin, probs = 0.975, na.rm = TRUE),
    .groups = "drop"
  )
# vast majority of total perturbation is kept in the origin state
ctmc_frac_shuffle_origin_summary
## ------------------------------------------------------------
## 8) Plots: age patterns of self-absorption and origin-share
## ------------------------------------------------------------

## 8.1 CR: fraction of shrink in self-transition within origin
ggplot(cr_frac_origin,
       aes(x = age, y = frac_stay_row,
           colour = hazard_perturbed)) +
  geom_line(alpha = 0.7) +
  facet_wrap(~ sex) +
  theme_minimal() +
  labs(
    title = "CR: fraction of counterperturbation absorbed by self-transition\n(within origin state)",
    x     = "Age",
    y     = "Fraction absorbed by self-transition"
  )

## 8.2 CTMC: same plot
ggplot(ctmc_frac_origin,
       aes(x = age, y = frac_stay_row,
           colour = hazard_perturbed)) +
  geom_line(alpha = 0.7) +
  facet_wrap(~ sex) +
  theme_minimal() +
  labs(
    title = "CTMC: fraction of counterperturbation absorbed by self-transition\n(within origin state)",
    x     = "Age",
    y     = "Fraction absorbed by self-transition"
  )

## 8.3 CTMC: fraction of *total shuffle* that remains in origin state
ggplot(ctmc_frac_shuffle_origin,
       aes(x = age, y = frac_shuffle_origin,
           colour = hazard_perturbed)) +
  geom_line(alpha = 0.7) +
  facet_wrap(~ sex) +
  theme_minimal() +
  labs(
    title = "CTMC: fraction of total perturbation shuffle within origin state",
    x     = "Age",
    y     = "Fraction of |Δp| in origin state"
  )

## ------------------------------------------------------------
## 10) DemoDecomp::horiuchi() comparison
##     (CR hazards vs CR-derived probabilities)
## ------------------------------------------------------------
probs_dec <-
  probs_cr |> 
  pivot_longer(-c(sex,age),names_to = "transition", values_to = "p", names_prefix = "p_") |> 
  pivot_wider(names_from = sex, values_from = p) |> 
  filter(!transition %in% c("HH","UU"))
  

probs_dec$cc <- horiuchi(f2dd,
                         probs_dec$m,
                         probs_dec$f,
                         df = probs_dec[,c("age","transition")],
                         dt = dt,
                         N = 20)
probs_dec$variant = "P2 horiuchi"


probs_dec2 <-
  probs_cr |> 
  pivot_longer(-c(sex,age), values_to = "p", names_to = "transition", names_prefix = "p_") |> 
  pivot_wider(names_from = sex, values_from = p) |> 
  mutate(p = (m+f) / 2,
         delta = f-m) 

probs_dec2<-
  probs_dec2 |>   
  group_modify(~s2t(data = .x,
                    expectancy = "h", 
                    init = c(H=1,U=0),
                    interval=dt)) |> filter(transition != "init") |> 
  mutate(age = age + 50) |> 
  left_join(probs_dec2, by = join_by(age,transition)) |> 
  mutate(cc = delta * effect) |> 
  select(age, transition, cc) |> 
  mutate(variant="P2 analytic")


dec_p2_compare <- bind_rows(probs_dec, probs_dec2)
# this plot should not match the p2 decomp results in the manuscript because the
# probabilities it's based on are different
dec_p2_compare |> 
  ggplot(aes(x=age,y=cc,color=transition, linetype=  variant)) +
  geom_line() +
  theme_minimal() +
  labs(title = "compare analytic P2 with Horiuchi P2",
       subtitle = "small differences due to sensitivity being evaluated at exact midpoint")

haz_cr_dec <-
  haz_cr %>%
  dplyr::select(
    sex, age, starts_with("h_")
  ) %>%
  tidyr::pivot_longer(
    cols = -c(sex,age),
    names_to = "transition",
    values_to = "h"
  ) %>%
  tidyr::pivot_wider(names_from = sex, values_from = h) %>%
  dplyr::arrange(age, transition)

# Note on computational cost here!
# hazard-based decomp using horiuchi needs to run expm() ca 2.1*10^7 times,
# and it's running f2(), and pivoting over 1e5 times, so this takes a very very
# long time to execute. Point being, we can compare to see convergence,
# but then know that this approach really isn't required to do a good job.

if (dt == 1/12){
  stop("This code might take >5 hours to execute\nYou can choose to read in its results in the next line if you want\nOtherwise, manually execute the next two lines in this block")
}
# started 20:12 23-12-2025
haz_cr_dec$cc <- horiuchi(f2dd_haz_ctmc, haz_cr_dec$m,haz_cr_dec$f,df=haz_cr_dec |> select(transition, age),dt=dt,init=c(H=1,U=0),N=20)
haz_cr_dec$variant = "hazard"

# uncomment this line below if you want to skip the 
# haz_cr_dec <- read_csv("haz_decomp_monthly.csv.gz")

# now do P2 to probs_cr
probs_cr_dec <-
  probs_cr |> 
  pivot_longer(-c(sex,age), values_to = "p", names_to = "transition", names_prefix = "p_") |> 
  pivot_wider(names_from = sex, values_from = p) |> 
  mutate(p = (m+f) / 2,
         delta = f-m) 

probs_cr_dec<-
  probs_cr_dec |>   
  group_modify(~s2t(data = .x,
                    expectancy = "h", 
                    init = c(H=1,U=0), 
                    interval = dt)) |> 
  filter(transition != "init") |> 
  mutate(age = age + 50) |> 
  left_join(probs_cr_dec, by = join_by(age,transition)) |> 
  mutate(cc = delta * effect) |> 
  select(age, transition, m,f,cc) |> 
  mutate(variant="P2 analytic")

haz_prob_compare <-
  haz_cr_dec |> 
  mutate(transition = substr(transition,3,4)) |> 
  bind_rows(probs_cr_dec)
if (dt==1){
  write_csv(haz_prob_compare,"haz_prob_decomp_annual.csv.gz")
} else {
  write_csv(haz_prob_compare,"haz_prob_decomp_monthly.csv.gz")
}
haz_prob_compare |> 
  ggplot(aes(x=age,y=cc, color = transition, linetype = variant)) +
  geom_line() +
  theme_minimal() +
  labs(y = "contribution to gap (years)")


haz_prob_compare |> 
  group_by(transition, variant) |> 
  summarize(cc = sum(cc, na.rm = TRUE)) |> 
  pivot_wider(names_from = variant, values_from = cc)

haz_prob_compare |> 
  mutate(sen = cc/(f-m)) |> 
  ggplot(aes(x=age,y=sen,color=transition))+
  geom_line()+
  facet_wrap(~variant)

# haz_cr (due to end-of-interval assumption)
# HLE      males     females    diff
# monthly  24.91891  26.36728   1.44837
# annual   25.41989  26.87229   1.4524
# diff     0.50098   0.50501

# probs_cr
# HLE      males     females    diff
# monthly  24.91891  26.36728   1.44837
# annual   25.41989  26.87229   1.4524

# probs (due to end-of-interval assumption)

results_annual <- read_csv("haz_prob_decomp_annual.csv.gz")
results_monthly <- read_csv("haz_prob_decomp_monthly.csv.gz")


make_decomp_margin_table <- function(results_annual,
                                     results_monthly,
                                     prob_variant = "P2 analytic",
                                     haz_variant  = "hazard",
                                     diff_label   = "Prob - Hazard",
                                     digits = 3,
                                     caption = "Marginal decomposition sums by transition under analytic P2 versus hazard-based Horiuchi decomposition (CTMC), for annual and monthly time steps.",
                                     label = "tab:p2_vs_hazard_margins") {
  stopifnot(all(c("transition","variant","cc") %in% names(results_annual)))
  stopifnot(all(c("transition","variant","cc") %in% names(results_monthly)))
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(xtable)
  })
  
  # helper: summarize -> wide with prob/haz + diff + total row
  summarise_one <- function(df) {
    wide <- df %>%
      group_by(transition, variant) %>%
      summarise(cc = sum(cc, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = variant, values_from = cc)
    
    # ensure columns exist
    if (!prob_variant %in% names(wide)) wide[[prob_variant]] <- NA_real_
    if (!haz_variant  %in% names(wide)) wide[[haz_variant]]  <- NA_real_
    
    wide <- wide %>%
      mutate(
        !!diff_label := .data[[prob_variant]] - .data[[haz_variant]]
      )
    
    # total row
    total <- wide %>%
      summarise(
        transition = "Total",
        !!prob_variant := sum(.data[[prob_variant]], na.rm = TRUE),
        !!haz_variant  := sum(.data[[haz_variant]],  na.rm = TRUE),
        !!diff_label   := sum(.data[[diff_label]],   na.rm = TRUE)
      )
    
    bind_rows(wide, total)
  }
  
  ann <- summarise_one(results_annual)
  mon <- summarise_one(results_monthly)
  
  # Combine side-by-side
  tab <- ann %>%
    rename(
      Annual_Prob = !!prob_variant,
      Annual_Haz  = !!haz_variant,
      Annual_Diff = !!diff_label
    ) %>%
    left_join(
      mon %>%
        rename(
          Monthly_Prob = !!prob_variant,
          Monthly_Haz  = !!haz_variant,
          Monthly_Diff = !!diff_label
        ),
      by = "transition"
    )
  
  # Optional: order transitions, keep Total last
  tab <- tab %>%
    mutate(is_total = transition == "Total") %>%
    arrange(is_total, transition) %>%
    select(-is_total)
  
  # Build xtable
  xt <- xtable(tab, caption = caption, label = label, digits = digits)
  
  # Vertical divider between Annual and Monthly blocks:
  # transition | Annual Prob Annual Haz Annual Diff || Monthly Prob Monthly Haz Monthly Diff
  align(xt) <- c("l", "l", "r", "r", "r", "|", "r", "r", "r")
  
  # Add a multicolumn header row (Annual / Monthly) + underline rules
  add.to.row <- list()
  add.to.row$pos <- list(0)
  add.to.row$command <- paste0(
    "\\hline\n",
    " & \\multicolumn{3}{c}{Annual} & \\multicolumn{3}{c}{Monthly} \\\\\n",
    " \\cline{2-4} \\cline{5-7}\n"
  )
  
  # Slightly nicer column names
  colnames(tab) <- c("Transition",
                     "P2", "Hazard", "Diff",
                     "P2", "Hazard", "Diff")
  xt <- xtable(tab, caption = caption, label = label, digits = digits)
  align(xt) <- c("l", "l", "r", "r", "r", "|", "r", "r", "r")
  
  print(xt,
        include.rownames = FALSE,
        booktabs = FALSE,   # you can set TRUE if you prefer booktabs styling
        hline.after = c(-1, 0, nrow(tab)), # top, after header, bottom
        add.to.row = add.to.row,
        sanitize.text.function = identity)
}

# Example usage:
 make_decomp_margin_table(results_annual, results_monthly,
                          prob_variant = "P2 analytic",
                          haz_variant = "hazard",
                          digits = 3)


