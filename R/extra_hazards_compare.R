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
## ------------------------------------------------------------
## 1) Read annual transition probabilities
## ------------------------------------------------------------
## transitions_lievre2003_annual.csv must have columns:
##   age, sex, HH, HU, HD, UH, UU, UD
## as in the main manuscript example.

probs_annual <- readr::read_csv("transitions_lievre2003_annual.csv") |>
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

## Just to check structure:
# dplyr::glimpse(probs_annual)

## Optionally: reconstruct monthly probabilities from Lievre eq. 33
## This reproduces the original source pipeline, but the main
## supplementary analyses below use the *annual* probabilities.

tab2 <- tibble(
  from_to  = c("HU",  "HD",   "UH",   "UD"),
  intercept = c(-11.9, -11.9, -1.4,  -5.4),
  age       = c(0.0853, 0.0812, -0.0358, 0.0207),
  sex       = c(-2.4,   3.0,   -0.2,   -2.5),
  age_sex   = c(0.0316, -0.0472, 0.0010, 0.0249)
)

## Equation 33 (Lievre et al. 2003): ALR of transition probs
eq33 <- function(tab2, .from_to = "HU", ages = 70:90, sex = 1) {
  coefs <- tab2 |> filter(from_to == .from_to)
  coefs$intercept +
    ages * coefs$age +
    sex  * coefs$sex +
    ages * sex * coefs$age_sex
}

## Example: monthly probabilities (not used later, but kept for reference)
ages_monthly <- seq(50, 110, by = 1/12)

fromH_f <- cbind(
  HU = eq33(tab2, "HU", ages_monthly, 1),
  HD = eq33(tab2, "HD", ages_monthly, 1)
) |>
  alrInv() |>
  as_tibble() |>
  rename(HH = V3) |>
  mutate(age = ages_monthly, sex = "f")

fromU_f <- cbind(
  UH = eq33(tab2, "UH", ages_monthly, 1),
  UD = eq33(tab2, "UD", ages_monthly, 1)
) |>
  alrInv() |>
  as_tibble() |>
  rename(UU = V3) |>
  mutate(age = ages_monthly, sex = "f")

fromH_m <- cbind(
  HU = eq33(tab2, "HU", ages_monthly, 0),
  HD = eq33(tab2, "HD", ages_monthly, 0)
) |>
  alrInv() |>
  as_tibble() |>
  rename(HH = V3) |>
  mutate(age = ages_monthly, sex = "m")

fromU_m <- cbind(
  UH = eq33(tab2, "UH", ages_monthly, 0),
  UD = eq33(tab2, "UD", ages_monthly, 0)
) |>
  alrInv() |>
  as_tibble() |>
  rename(UU = V3) |>
  mutate(age = ages_monthly, sex = "m")

probs_monthly <- fromU_m |>
  left_join(fromH_m, by = join_by(age, sex)) |>
  bind_rows(
    fromU_f |>
      left_join(fromH_f, by = join_by(age, sex))
  ) |>
  relocate(age, .before = 1) |>
  relocate(sex, .before = 1) |>
  pivot_longer(
    cols      = -c(sex, age),
    names_to  = "transition",
    values_to = "p"
  ) |>
  mutate(transition = paste0("p_", transition)) |>
  pivot_wider(
    names_from  = "transition",
    values_from = "p"
  )
# Note: in the following, you could reproduce this using probs_monthly
# instead of probs, in which case you'll need to always set dt = 1/12 
# instead of dt = 1.


## ------------------------------------------------------------
## 2) Functions: CTMC (logm/expm) and CR hazard mappings
## ------------------------------------------------------------

## 2.1 CTMC hazards via matrix log
## probs: tibble/data.frame with columns
##   sex, age, p_HH, p_HU, p_HD, p_UH, p_UU, p_UD
probs_to_hazards_ctmc_logm <- function(probs, dt = 1) {
  prob_mat <- as.data.frame(
    probs[, c("p_HH", "p_HU", "p_HD", "p_UH", "p_UU", "p_UD")]
  )
  
  one_row_to_haz <- function(v) {
    p_HH <- v[1]; p_HU <- v[2]; p_HD <- v[3]
    p_UH <- v[4]; p_UU <- v[5]; p_UD <- v[6]
    
    ## 3x3 P (standard Markov orientation)
    P <- matrix(
      c(
        p_HH, p_HU, p_HD,
        p_UH, p_UU, p_UD,
        0,    0,    1
      ),
      nrow = 3, byrow = TRUE
    )
    
    ## matrix log -> generator Q
    L <- expm::logm(P)
    if (any(abs(Im(L)) > 1e-8)) {
      warning("Non-negligible imaginary parts in logm(P); check probabilities.")
    }
    Q <- Re(L) / dt
    
    ## force row sums 0 (numerical cleanup)
    rs <- rowSums(Q)
    diag(Q) <- diag(Q) - rs
    
    c(
      h_HU = Q[1, 2],
      h_HD = Q[1, 3],
      h_UH = Q[2, 1],
      h_UD = Q[2, 3]
    )
  }
  
  haz_mat <- t(apply(prob_mat, 1, one_row_to_haz))
  haz_df  <- as.data.frame(haz_mat)
  
  dplyr::bind_cols(
    probs |> select(sex, age),
    haz_df
  )
}

## 2.2 CTMC: hazards -> probs via expm(Q)
## haz: tibble/data.frame with columns
##   sex, age, h_HU, h_HD, h_UH, h_UD
hazards_to_probs_ctmc <- function(haz, dt = 1) {
  haz_mat <- as.data.frame(haz[, c("h_HU", "h_HD", "h_UH", "h_UD")])
  
  one_row_to_probs <- function(v) {
    h_HU <- v[1]; h_HD <- v[2]
    h_UH <- v[3]; h_UD <- v[4]
    
    Q <- matrix(
      c(
        -(h_HU + h_HD), h_HU,            h_HD,
        h_UH,           -(h_UH + h_UD),  h_UD,
        0,              0,               0
      ),
      nrow = 3, byrow = TRUE
    )
    
    P <- expm::expm(Q * dt)
    
    c(
      p_HH = P[1, 1],
      p_HU = P[1, 2],
      p_HD = P[1, 3],
      p_UH = P[2, 1],
      p_UU = P[2, 2],
      p_UD = P[2, 3]
    )
  }
  
  prob_mat <- t(apply(haz_mat, 1, one_row_to_probs))
  prob_df  <- as.data.frame(prob_mat)
  
  dplyr::bind_cols(
    haz |> select(sex, age),
    prob_df
  )
}

## 2.3 Competing-risks (CR) hazards from probabilities
## Rowwise: treat each origin state as a CR system with
## total exit rate r = -log(p_ii)/dt and hazard shares
probs_to_hazards_cr <- function(probs, dt = 1) {
  probs |>
    rowwise() |>
    mutate(
      ## H row
      r_H  = ifelse(p_HH > 0 & p_HH < 1, -log(p_HH) / dt, NA_real_),
      h_HU = ifelse(!is.na(r_H),
                    r_H * p_HU / (1 - p_HH),
                    NA_real_),
      h_HD = ifelse(!is.na(r_H),
                    r_H * p_HD / (1 - p_HH),
                    NA_real_),
      ## U row
      r_U  = ifelse(p_UU > 0 & p_UU < 1, -log(p_UU) / dt, NA_real_),
      h_UH = ifelse(!is.na(r_U),
                    r_U * p_UH / (1 - p_UU),
                    NA_real_),
      h_UD = ifelse(!is.na(r_U),
                    r_U * p_UD / (1 - p_UU),
                    NA_real_)
    ) |>
    ungroup() |>
    select(sex, age, h_HU, h_HD, h_UH, h_UD)
}

## 2.4 CR: hazards -> probs using standard CR formulas
hazards_to_probs_cr <- function(haz, dt = 1) {
  haz |>
    rowwise() |>
    mutate(
      ## H row
      r_H  = h_HU + h_HD,
      p_HH = exp(-r_H * dt),
      p_HU = ifelse(r_H > 0,
                    h_HU / r_H * (1 - p_HH),
                    0),
      p_HD = ifelse(r_H > 0,
                    h_HD / r_H * (1 - p_HH),
                    0),
      ## U row
      r_U  = h_UH + h_UD,
      p_UU = exp(-r_U * dt),
      p_UH = ifelse(r_U > 0,
                    h_UH / r_U * (1 - p_UU),
                    0),
      p_UD = ifelse(r_U > 0,
                    h_UD / r_U * (1 - p_UU),
                    0)
    ) |>
    ungroup() |>
    select(sex, age, p_HH, p_HU, p_HD, p_UH, p_UU, p_UD)
}

## ------------------------------------------------------------
## 3) Compare CTMC-logm hazards vs CR hazards
##    and original probs vs CR-derived probs
## ------------------------------------------------------------

## 3.1 Hazards: logm-based CTMC vs CR
haz_ctmc_logm <- probs_to_hazards_ctmc_logm(probs_annual, dt = 1)
haz_cr        <- probs_to_hazards_cr(probs_annual, dt = 1)

## Example: look at h_HD and h_HU by sex and age
haz_long <- bind_rows(
  haz_ctmc_logm |>
    mutate(variant = "CTMC_logm"),
  haz_cr |>
    mutate(variant = "CR_rowwise")
) |>
  pivot_longer(
    cols      = c(h_HU, h_HD, h_UH, h_UD),
    names_to  = "hazard",
    values_to = "h"
  )

## Plot: hazard comparison (e.g., h_HD)
ggplot(
  haz_long |> filter(hazard == "h_HD"),
  aes(x = age, y = h, colour = variant, linetype = variant)
) +
  geom_line() +
  facet_wrap(~ sex) +
  theme_minimal() +
  labs(
    title = "Death hazard from H (h_HD): CTMC logm vs CR rowwise",
    y     = "h_HD",
    x     = "Age",
    subtitle = "for females, CTMC using logm() simply does not give plausible results with this data"
  )

## 3.2 Probabilities: original vs CR-CTMC-implied probabilities
probs_cr <- hazards_to_probs_ctmc(haz_cr, dt = 1)
# Note: you can get back the original probabilities exatly
# using hazards_to_probs_cr(), but then we can't do this 
# comparison!

probs_orig_long <- probs_annual |>
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
  ggplot(aes(x = age,color=transition)) +
    geom_line(aes(y = p_orig)) +
    geom_line(aes(y = p_cr), linetype = "dashed") +
    facet_wrap(~ sex) +
    theme_minimal() +
    labs(
      title = "Original vs CR-implied probabilities",
      y     = "prob",
      x     = "Age",
      subtitle = "notable differences in probabilities, but CR probabilities are\nstill in the plausible range and usable for testing",
      caption = "(we take the probs->cr hazards->ctmc probs result as the point of comparison)"
    )

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
## 4) Generic perturbation helpers (hazard-level perturbation)
## ------------------------------------------------------------

## Perturb one hazard column by a factor (1 + eps).
# note, this perturbs the whole age pattern of the given hazard, and it's
# not the same as perturbing just one age of that hazard. We need to add this
# more specific check in a second iteration to disentangle direct and indirect
# effects?
perturb_one_hazard <- function(haz_df, hazard_name, eps = 1e-4) {
  haz_df[[hazard_name]] <- haz_df[[hazard_name]] * (1 + eps)
  haz_df
}

## Run perturbation experiment using a mapping
##   hazards -> probabilities
run_hazard_perturbation_generic <- function(haz_df,
                                            hazard_to_probs_fun,
                                            dt  = 1,
                                            eps = 1e-4) {
  ## Baseline probabilities implied by hazards
  probs_orig <- hazard_to_probs_fun(haz_df, dt = dt)
  
  probs_orig_long <- probs_orig |>
    pivot_longer(
      cols      = starts_with("p_"),
      names_to  = "transition",
      values_to = "p_orig"
    )
  
  hazards_to_perturb <- c("h_HU", "h_HD", "h_UH", "h_UD")
  
  out_list <- lapply(hazards_to_perturb, function(hname) {
    haz_pert <- perturb_one_hazard(haz_df, hname, eps = eps)
    
    probs_pert <- hazard_to_probs_fun(haz_pert, dt = dt)
    
    probs_pert_long <- probs_pert |>
      pivot_longer(
        cols      = starts_with("p_"),
        names_to  = "transition",
        values_to = "p_pert"
      ) |>
      mutate(hazard_perturbed = hname)
    
    ## Join with baseline p_orig
    probs_orig_long |>
      inner_join(
        probs_pert_long,
        by = c("sex", "age", "transition")
      )
  })
  
  bind_rows(out_list) |>
    mutate(
      dt  = dt,
      eps = eps
    ) |>
    select(sex, age, hazard_perturbed, transition, p_orig, p_pert, dt, eps)
}

## Convenience wrappers:
run_hazard_perturbation_cr <- function(haz_df, dt = 1, eps = 1e-4) {
  run_hazard_perturbation_generic(haz_df, hazards_to_probs_cr, dt = dt, eps = eps)
}

run_hazard_perturbation_ctmc <- function(haz_df, dt = 1, eps = 1e-4) {
  run_hazard_perturbation_generic(haz_df, hazards_to_probs_ctmc, dt = dt, eps = eps)
}

## ------------------------------------------------------------
## 5) Counterperturbation summaries
## ------------------------------------------------------------

## 5.1 Global: among all shrinking probabilities (both rows),
## what fraction of total shrink is in self-transitions p_HH + p_UU?
compute_counterpert_fraction_staying_global <- function(df) {
  df %>%
    mutate(
      delta = p_pert - p_orig,
      origin_state = dplyr::case_when(
        hazard_perturbed %in% c("h_HU", "h_HD") ~ "H",
        hazard_perturbed %in% c("h_UH", "h_UD") ~ "U",
        TRUE ~ NA_character_
      )
    ) %>%
    group_by(sex, age, hazard_perturbed) %>%  # full system for this perturbation
    summarise(
      total_shrink = sum(abs(delta[delta < 0]), na.rm = TRUE),
      self_shrink  = {
        os <- unique(origin_state)
        if (length(os) != 1 || is.na(os)) {
          NA_real_
        } else {
          stay_transition <- if (os == "H") "p_HH" else "p_UU"
          sum(
            abs(delta[transition == stay_transition & delta < 0]),
            na.rm = TRUE
          )
        }
      },
      frac_self_global = ifelse(total_shrink > 0,
                                self_shrink / total_shrink,
                                NA_real_),
      .groups = "drop"
    )
}

## 5.2 Within-row: only consider shrinkage in the origin row (H or U),
## and ask what fraction of that is in the self-transition (p_HH or p_UU)
compute_counterpert_fraction_staying_origin_only <- function(df) {
  if (!"hazard_perturbed" %in% names(df)) {
    stop("Data frame must contain a 'hazard_perturbed' column.")
  }
  
  df |>
    mutate(delta = p_pert - p_orig) |>
    group_by(sex, age, hazard_perturbed) |>
    group_modify(~{
      h <- .y$hazard_perturbed[[1]]
      
      ## map hazard -> origin row & staying transition
      origin_transitions <- if (h %in% c("h_HU", "h_HD")) {
        c("p_HH", "p_HU", "p_HD")  # from H
      } else if (h %in% c("h_UH", "h_UD")) {
        c("p_UH", "p_UU", "p_UD")  # from U
      } else {
        stop("Unknown hazard_perturbed: ", h)
      }
      
      stay_transition <- if (h %in% c("h_HU", "h_HD")) "p_HH" else "p_UU"
      
      row_df <- .x |> filter(transition %in% origin_transitions)
      
      total_shrink_row <- sum(abs(row_df$delta[row_df$delta < 0]), na.rm = TRUE)
      stay_shrink_row  <- sum(
        abs(row_df$delta[row_df$transition == stay_transition &
                           row_df$delta < 0]),
        na.rm = TRUE
      )
      
      tibble(
        total_shrink_row = total_shrink_row,
        stay_shrink_row  = stay_shrink_row,
        frac_stay_row    = if (total_shrink_row > 0)
          stay_shrink_row / total_shrink_row
        else NA_real_
      )
    }) |>
    ungroup()
}

## 5.3 Fraction of total "shuffle" that happens in the origin state
compute_shuffle_fraction_origin <- function(df) {
  if (!all(c("hazard_perturbed", "transition", "p_orig", "p_pert") %in% names(df))) {
    stop("Data frame must contain hazard_perturbed, transition, p_orig, p_pert.")
  }
  
  df |>
    mutate(
      delta = p_pert - p_orig,
      origin_row = dplyr::case_when(
        hazard_perturbed %in% c("h_HU", "h_HD") ~ "H",
        hazard_perturbed %in% c("h_UH", "h_UD") ~ "U",
        TRUE ~ NA_character_
      ),
      in_origin_row = dplyr::case_when(
        origin_row == "H" & transition %in% c("p_HH", "p_HU", "p_HD") ~ TRUE,
        origin_row == "U" & transition %in% c("p_UH", "p_UU", "p_UD") ~ TRUE,
        TRUE ~ FALSE
      )
    ) |>
    group_by(sex, age, hazard_perturbed) |>
    summarise(
      total_shuffle       = sum(abs(delta), na.rm = TRUE),
      origin_shuffle      = sum(abs(delta[in_origin_row]), na.rm = TRUE),
      frac_shuffle_origin = ifelse(total_shuffle > 0,
                                   origin_shuffle / total_shuffle,
                                   NA_real_),
      .groups = "drop"
    )
}

## ------------------------------------------------------------
## 6) Run perturbation experiments: CR and CTMC (based on CR hazards)
## ------------------------------------------------------------

## We use CR-based hazards haz_cr as our starting point.
## These imply probabilities probs_cr that are demographically
## plausible but do not exactly match probs_annual.

haz_rowwise <- haz_cr

## Choose perturbation size and dt
eps_pert <- 1e-3
dt_pert  <- 1

## 6.1 CR mapping hazards -> probs
perturb_results_cr <- run_hazard_perturbation_cr(
  haz_df = haz_rowwise,
  dt     = dt_pert,
  eps    = eps_pert
)

## 6.2 CTMC mapping hazards -> probs via expm(Q)
perturb_results_ctmc <- run_hazard_perturbation_ctmc(
  haz_df = haz_rowwise,
  dt     = dt_pert,
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
# The vast majority of total counterperturbation is taken by the self-transition
# located in the same origin state as the hazard perturbation.
ctmc_frac_global

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
## 9) Taylor-based approximation addendum
## ------------------------------------------------------------

## Matrix exponential via truncated Taylor series:
##   expm(A t) ≈ I + A t + (A t)^2 / 2! + ... + (A t)^n / n!
expm_taylor <- function(A, t = 1, n = 5) {
  if (!is.matrix(A) || nrow(A) != ncol(A)) {
    stop("A must be square.")
  }
  d  <- nrow(A)
  At <- A * t
  out  <- diag(1, d)  # k = 0 term
  term <- diag(1, d)  # accumulator
  
  if (n == 0) return(out)
  
  for (k in 1:n) {
    term <- term %*% At / k  # (At)^k / k!
    out  <- out + term
  }
  out
}

## CTMC hazards -> probs via Taylor
hazards_to_probs_ctmc_taylor <- function(haz, dt = 1, n_terms = 5) {
  haz_mat <- as.data.frame(haz[, c("h_HU", "h_HD", "h_UH", "h_UD")])
  
  one_row_to_probs <- function(v) {
    h_HU <- v[1]; h_HD <- v[2]
    h_UH <- v[3]; h_UD <- v[4]
    
    Q <- matrix(
      c(
        -(h_HU + h_HD), h_HU,            h_HD,
        h_UH,           -(h_UH + h_UD),  h_UD,
        0,              0,               0
      ),
      nrow = 3, byrow = TRUE
    )
    
    P <- expm_taylor(Q, t = dt, n = n_terms)
    
    c(
      p_HH = P[1, 1],
      p_HU = P[1, 2],
      p_HD = P[1, 3],
      p_UH = P[2, 1],
      p_UU = P[2, 2],
      p_UD = P[2, 3]
    )
  }
  
  prob_mat <- t(apply(haz_mat, 1, one_row_to_probs))
  prob_df  <- as.data.frame(prob_mat)
  
  dplyr::bind_cols(
    haz |> select(sex, age),
    prob_df
  )
}

run_hazard_perturbation_ctmc_taylor <- function(haz_df,
                                                dt      = 1,
                                                eps     = 1e-4,
                                                n_terms = 5) {
  hazard_to_probs_fun <- function(h, dt) {
    hazards_to_probs_ctmc_taylor(h, dt = dt, n_terms = n_terms)
  }
  run_hazard_perturbation_generic(haz_df, hazard_to_probs_fun, dt = dt, eps = eps)
}

## This addendum isn't fully worked out here, but the 
# above diagnostics could be run under different n. It's 
# always interesting to see how results approach
# expm() as n increases. Low order n is just CR, whereas
# high order n approaches expm(). We take expm() results
# as the canonical ones here.

## ------------------------------------------------------------
## 10) DemoDecomp::horiuchi() comparison
##     (CR hazards vs CR-derived probabilities)
## ------------------------------------------------------------
probs_annual_dec <-
  probs_annual |> 
  pivot_longer(-c(sex,age),names_to = "transition", values_to = "p", names_prefix = "p_") |> 
  pivot_wider(names_from = sex, values_from = p) |> 
  filter(!transition %in% c("HH","UU"))
  
f2dd <- function(probs_vec, df, init = c(H=1,U=0), expectancy = "h"){

  df$p <- probs_vec
  df   <- pivot_wider(df, names_from = transition, values_from = p)
  f2(hd = df$HD,
     hu = df$HU,
     ud = df$UD,
     uh = df$UH,
     init = init,
     expectancy = expectancy)
}

f2dd(probs_vec = probs_annual_dec$f,
                 df=probs_annual_dec[,c("age","transition")])
probs_annual_dec$cc <- horiuchi(f2dd,
                                probs_annual_dec$m,
                                probs_annual_dec$f,
                                df = probs_annual_dec[,c("age","transition")],
                                N = 20)
probs_annual_dec$variant = "P2 horiuchi"


probs_annual_dec2 <-
probs_annual |> 
  pivot_longer(-c(sex,age), values_to = "p", names_to = "transition", names_prefix = "p_") |> 
  pivot_wider(names_from = sex, values_from = p) |> 
  mutate(p = (m+f) / 2,
         delta = f-m) 

probs_annual_dec2<-
  probs_annual_dec2 |>   
  group_modify(~s2t(data = .x,expectancy = "h", init = c(H=1,U=0))) |> 
  filter(transition != "init") |> 
  mutate(age = age + 50) |> 
  left_join(probs_annual_dec2, by = join_by(age,transition)) |> 
  mutate(cc = delta * effect) |> 
  select(age, transition, cc) |> 
  mutate(variant="P2 analytic")


dec_p2_compare <- bind_rows(probs_annual_dec, probs_annual_dec2)
# this plot should match the p2 decomp results in the manuscript
dec_p2_compare |> 
  ggplot(aes(x=age,y=cc,color=transition, linetype=  variant)) +
  geom_line() +
  theme_minimal() +
  labs(title = "compare analytic P2 with Horiuchi P2",
       subtitle = "small differences due to sensitivity being evaluated at exact midpoint")

# ------------------------------------------------------------------ #
# now have established that Horiuchi can match our analytic results  #
# this was done using the manuscript probabilities; Now we can       #
# compare hazard-based Horiuchi with P2 probability-based Horiuchi   #
# for the case of haz_cr and probs_cr (which correspond)             #
## ----------------------------------------------------------------- #
## End of script
## ------------------------------------------------------------
