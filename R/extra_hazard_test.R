## -----------------------------------------------------------
## Setup
## -----------------------------------------------------------
library(expm)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(tidyverse)
library(compositions)

# read annual probs
probs <- read_csv("transitions_lievre2003_annual.csv") |> 
  pivot_longer(-c(age,sex), names_to = "transition",values_to = "p") |> 
  mutate(transition = paste0("p_",transition)) |> 
  pivot_wider(names_from = transition, values_from = p)

# derive monthly probs
tab2 <- tibble(from_to = c("HU","HD","UH","UD"),
               intercept = c(-11.9,-11.9,-1.4,-5.4),
               age = c(0.0853,0.0812,-0.0358,0.0207),
               sex = c(-2.4,3,-.2,-2.5),
               age_sex = c(0.0316,-0.0472,0.0010,0.0249))

# a custom function that implements equation 33 from Lievre et al 2003
# this converts coefficients to the ALR of HD, HU, UH, UD, using HH and UU, 
# respectively in the denominator.
eq33 <- function(tab2, .from_to = "HU", ages = 70:90, sex = 1){
  coefs <- tab2 |>
    filter(from_to == .from_to)
  coefs$intercept + ages * coefs$age + sex * coefs$sex + ages * sex *  coefs$age_sex
}


# get monthly transitions, one sex and one origin state at a time.
# eq33, which returns ALR transforms of HU and HD (UH and UD) such that
# we can retrieve all six transitions using the alrInv function.
# these will be monthly transitions.
ages <- seq(50,110,by=1/12)

# for women
fromHf <- 
  cbind(HU = eq33(tab2, "HU", ages, 1),
        HD = eq33(tab2, "HD", ages, 1)) %>% 
  alrInv() %>% 
  as_tibble() %>% 
  rename(HH = V3) %>% 
  mutate(age = ages,
         sex = "f") 

fromUf <- 
  cbind(UH = eq33(tab2, "UH", ages, 1),
        UD = eq33(tab2, "UD", ages, 1)) %>% 
  alrInv() %>% 
  as_tibble() %>% 
  rename(UU = V3) %>% 
  mutate(age = ages,
         sex = "f")

fromHm <- 
  cbind(HU = eq33(tab2, "HU", ages, 0),
        HD = eq33(tab2, "HD", ages, 0)) %>% 
  alrInv() %>% 
  as_tibble() %>% 
  rename(HH = V3) %>% 
  mutate(age = ages,
         sex = "m") 

fromUm <- 
  cbind(UH = eq33(tab2, "UH", ages, 0),
        UD = eq33(tab2, "UD", ages, 0)) %>% 
  alrInv() %>% 
  as_tibble() %>% 
  rename(UU = V3) %>% 
  mutate(age = ages,
         sex = "m")

probs_monthly <-
  fromUm |> 
  left_join(fromHm,by = join_by(age, sex)) |> 
  bind_rows(fromUf |> 
              left_join(fromHf,by = join_by(age, sex))) |> 
  relocate(age, .before = 1) |> 
  relocate(sex, .before = 1) |> 
  pivot_longer(-c(sex,age), values_to = "p", names_to = "transition") |> 
  mutate(transition = paste0("p_",transition)) |> 
  pivot_wider(names_from = transition, values_from = p)
# --------------------------------------------------

# probs: tibble/data.frame with columns
# sex, age, p_HH, p_HU, p_HD, p_UH, p_UU, p_UD
probs_to_hazards_ctmc_logm <- function(probs, dt = 1) {
  
  prob_mat <- as.data.frame(probs[, c("p_HH", "p_HU", "p_HD",
                                      "p_UH", "p_UU", "p_UD")])
  
  one_row_to_haz <- function(v) {
    p_HH <- v[1]; p_HU <- v[2]; p_HD <- v[3]
    p_UH <- v[4]; p_UU <- v[5]; p_UD <- v[6]
    
    # 3x3 P
    P <- matrix(
      c(
        p_HH, p_HU, p_HD,
        p_UH, p_UU, p_UD,
        0,    0,    1
      ),
      nrow = 3, byrow = TRUE
    )
    
    # matrix log -> Q
    L <- expm::logm(P)
    if (any(abs(Im(L)) > 1e-8)) {
      warning("Non-negligible imaginary parts in logm(P); check probabilities.")
    }
    Q <- Re(L) / dt
    
    # enforce row sums 0 via diagonals
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
  
  bind_cols(
    probs %>% select(sex, age),
    haz_df
  )
}

# haz: tibble/data.frame with columns
# sex, age, h_HU, h_HD, h_UH, h_UD
hazards_to_probs_ctmc <- function(haz, dt = 1) {
  
  haz_mat <- as.data.frame(haz[, c("h_HU", "h_HD", "h_UH", "h_UD")])
  
  one_row_to_probs <- function(v) {
    h_HU <- v[1]; h_HD <- v[2]
    h_UH <- v[3]; h_UD <- v[4]
    
    Q <- matrix(
      c(
        -(h_HU + h_HD), h_HU,            h_HD,
        h_UH,          -(h_UH + h_UD),  h_UD,
        0,              0,              0
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
  
  bind_cols(
    haz %>% select(sex, age),
    prob_df
  )
}


# probs: tibble with sex, age, p_HH, p_HU, p_HD, p_UH, p_UU, p_UD
probs_to_hazards_cr <- function(probs, dt = 1) {
  probs %>%
    rowwise() %>%
    mutate(
      # H row
      r_H  = ifelse(p_HH > 0 & p_HH < 1,
                    -log(p_HH) / dt,
                    NA_real_),
      h_HU = ifelse(!is.na(r_H),
                    r_H * p_HU / (1 - p_HH),
                    NA_real_),
      h_HD = ifelse(!is.na(r_H),
                    r_H * p_HD / (1 - p_HH),
                    NA_real_),
      # U row
      r_U  = ifelse(p_UU > 0 & p_UU < 1,
                    -log(p_UU) / dt,
                    NA_real_),
      h_UH = ifelse(!is.na(r_U),
                    r_U * p_UH / (1 - p_UU),
                    NA_real_),
      h_UD = ifelse(!is.na(r_U),
                    r_U * p_UD / (1 - p_UU),
                    NA_real_)
    ) %>%
    ungroup() %>%
    select(sex, age, h_HU, h_HD, h_UH, h_UD)
}

hazards_to_probs_cr <- function(haz, dt = 1) {
  haz %>%
    rowwise() %>%
    mutate(
      # H row
      r_H  = h_HU + h_HD,
      p_HH = exp(-r_H * dt),
      p_HU = ifelse(r_H > 0,
                    h_HU / r_H * (1 - p_HH),
                    0),
      p_HD = ifelse(r_H > 0,
                    h_HD / r_H * (1 - p_HH),
                    0),
      # U row
      r_U  = h_UH + h_UD,
      p_UU = exp(-r_U * dt),
      p_UH = ifelse(r_U > 0,
                    h_UH / r_U * (1 - p_UU),
                    0),
      p_UD = ifelse(r_U > 0,
                    h_UD / r_U * (1 - p_UU),
                    0)
    ) %>%
    ungroup() %>%
    select(sex, age, p_HH, p_HU, p_HD, p_UH, p_UU, p_UD)
}
# Perturb one hazard column by (1 + eps)
perturb_one_hazard <- function(haz_df, hazard_name, eps = 1e-4) {
  haz_df %>%
    mutate(
      !!hazard_name := !!sym(hazard_name) * (1 + eps)
    )
}
run_hazard_perturbation_ctmc <- function(haz_df, dt = 1, eps = 1e-4) {
  
  # 1) Baseline CTMC probabilities implied by hazards
  probs_orig <- hazards_to_probs_ctmc(haz_df, dt = dt)
  
  probs_orig_long <- probs_orig %>%
    pivot_longer(
      cols      = starts_with("p_"),
      names_to  = "transition",
      values_to = "p_orig"
    )
  
  # 2) For each hazard, perturb and recompute probs
  hazards_to_perturb <- c("h_HU", "h_HD", "h_UH", "h_UD")
  
  out_list <- lapply(hazards_to_perturb, function(hname) {
    haz_pert <- perturb_one_hazard(haz_df, hname, eps = eps)
    
    probs_pert <- hazards_to_probs_ctmc(haz_pert, dt = dt)
    
    probs_pert_long <- probs_pert %>%
      pivot_longer(
        cols      = starts_with("p_"),
        names_to  = "transition",
        values_to = "p_pert"
      ) %>%
      mutate(hazard_perturbed = hname)
    
    # join with baseline p_orig
    probs_orig_long %>%
      inner_join(
        probs_pert_long,
        by = c("sex", "age", "transition")
      )
  })
  
  bind_rows(out_list) %>%
    mutate(
      dt  = dt,
      eps = eps
    ) %>%
    select(sex, age, hazard_perturbed, transition, p_orig, p_pert, dt, eps)
}
# Global: among all shrinking probabilities (both rows),
# what fraction of total shrink is in p_HH + p_UU?
compute_counterpert_fraction_staying_global <- function(df) {
  df %>%
    mutate(delta = p_pert - p_orig) %>%
    group_by(sex, age, hazard_perturbed) %>%
    summarise(
      total_shrink = sum(abs(delta[delta < 0]), na.rm = TRUE),
      stay_shrink  = sum(abs(delta[transition %in% c("p_HH","p_UU") & delta < 0]),
                         na.rm = TRUE),
      frac_stay    = ifelse(total_shrink > 0, stay_shrink / total_shrink, NA_real_),
      .groups = "drop"
    )
}

# Within-row: only consider shrinkage in the origin row (H or U),
# and ask what fraction of that is in p_HH (if from H) or p_UU (if from U)
compute_counterpert_fraction_staying_origin_only <- function(df) {
  # sanity check
  if (!"hazard_perturbed" %in% names(df)) {
    stop("Data frame must contain a 'hazard_perturbed' column.")
  }
  
  df %>%
    mutate(delta = p_pert - p_orig) %>%
    group_by(sex, age, hazard_perturbed) %>%
    group_modify(~{
      # .x = rows in this group
      # .y = one-row tibble with group keys
      h <- .y$hazard_perturbed[[1]]
      
      # map hazard -> origin row & staying transition
      origin_transitions <- if (h %in% c("h_HU", "h_HD")) {
        c("p_HH", "p_HU", "p_HD")   # from H
      } else if (h %in% c("h_UH", "h_UD")) {
        c("p_UH", "p_UU", "p_UD")   # from U
      } else {
        stop("Unknown hazard_perturbed: ", h)
      }
      
      stay_transition <- if (h %in% c("h_HU", "h_HD")) "p_HH" else "p_UU"
      
      row_df <- .x %>% filter(transition %in% origin_transitions)
      
      total_shrink_row <- sum(abs(row_df$delta[row_df$delta < 0]), na.rm = TRUE)
      stay_shrink_row  <- sum(abs(row_df$delta[row_df$transition == stay_transition &
                                                 row_df$delta < 0]),
                              na.rm = TRUE)
      
      tibble(
        total_shrink_row = total_shrink_row,
        stay_shrink_row  = stay_shrink_row,
        frac_stay_row    = if (total_shrink_row > 0)
          stay_shrink_row / total_shrink_row else NA_real_
      )
    }) %>%
    ungroup()
}
# how much or total perturbation happens in the origin?
compute_shuffle_fraction_origin <- function(df) {
  # df: output of run_hazard_perturbation_ctmc()
  # columns: sex, age, hazard_perturbed, transition, p_orig, p_pert, dt, eps
  
  if (!all(c("hazard_perturbed", "transition", "p_orig", "p_pert") %in% names(df))) {
    stop("Data frame must contain hazard_perturbed, transition, p_orig, p_pert.")
  }
  
  df %>%
    mutate(
      delta = p_pert - p_orig,
      
      # identify origin row for each perturbed hazard
      origin_row = dplyr::case_when(
        hazard_perturbed %in% c("h_HU", "h_HD") ~ "H",
        hazard_perturbed %in% c("h_UH", "h_UD") ~ "U",
        TRUE ~ NA_character_
      ),
      
      # flag transitions that belong to the origin row
      in_origin_row = dplyr::case_when(
        origin_row == "H" & transition %in% c("p_HH", "p_HU", "p_HD") ~ TRUE,
        origin_row == "U" & transition %in% c("p_UH", "p_UU", "p_UD") ~ TRUE,
        TRUE ~ FALSE
      )
    ) %>%
    group_by(sex, age, hazard_perturbed) %>%
    summarise(
      total_shuffle       = sum(abs(delta), na.rm = TRUE),
      origin_shuffle      = sum(abs(delta[in_origin_row]), na.rm = TRUE),
      frac_shuffle_origin = ifelse(total_shuffle > 0,
                                   origin_shuffle / total_shuffle,
                                   NA_real_),
      .groups = "drop"
    )
}

# 1) Get CR-based hazards from probs
haz_rowwise <- probs_to_hazards_cr(probs, dt = 1)

# 2) Run CTMC perturbation from those hazards
perturb_results <- run_hazard_perturbation_ctmc(
  haz_df = haz_rowwise,
  dt     = 1,
  eps    = 1e-3  # or 1e-4 etc.
)

# 3) Global fraction: where do negative changes live?
global_frac <- compute_counterpert_fraction_staying_global(perturb_results)

# 4) Within-origin-row fraction
row_frac <- compute_counterpert_fraction_staying_origin_only(perturb_results)

origin_frac <- compute_shuffle_fraction_origin(perturb_results)



# 5) Summaries for the response letter, e.g.:
row_frac %>%
  group_by(hazard_perturbed) %>%
  summarise(
    min_frac_row    = min(frac_stay_row, na.rm = TRUE),
    median_frac_row = median(frac_stay_row, na.rm = TRUE),
    max_frac_row    = max(frac_stay_row, na.rm = TRUE)
  )
global_frac %>%
  group_by(hazard_perturbed) %>%
  summarise(
    min_frac_gl    = min(frac_stay, na.rm = TRUE),
    median_frac_gl = median(frac_stay, na.rm = TRUE),
    max_frac_gl    = max(frac_stay, na.rm = TRUE)
  )



row_frac |> 
  ggplot(aes(x=age, y = frac_stay_row, color = hazard_perturbed)) +
  geom_line() +
  facet_wrap(~sex)

global_frac |> 
  ggplot(aes(x=age, y = frac_stay, color = hazard_perturbed)) +
  geom_line() +
  facet_wrap(~sex)

origin_frac |> 
  ggplot(aes(x=age, y = frac_shuffle_origin, color = hazard_perturbed)) +
  geom_line() +
  facet_wrap(~sex)


expm_taylor <- function(A, t = 1, n = 10) {
  # A: square matrix
  # t: time step
  # n: highest order term in the Taylor series
  #
  # returns: sum_{k=0}^n (A t)^k / k!
  
  if (!is.matrix(A) || nrow(A) != ncol(A)) {
    stop("A must be a square matrix")
  }
  
  d <- nrow(A)
  At <- A * t
  
  # start with k = 0 term: I
  out  <- diag(1, d)
  term <- diag(1, d)
  
  if (n == 0) return(out)
  
  for (k in 1:n) {
    term <- term %*% At / k   # (At)^k / k!
    out  <- out + term
  }
  
  out
}
P_true   <- expm::expm(Q * dt)
P_approx <- expm_taylor(Q, t = dt, n = 3)  # for example

max(abs(P_true - P_approx))


test_expm_taylor_single <- function(Q, dt = 1, n = 5) {
  P_true   <- expm::expm(Q * dt)
  P_taylor <- expm_taylor(Q, t = dt, n = n)
  
  tibble::tibble(
    max_abs_err = max(abs(P_true - P_taylor)),
    max_rel_err = max(abs((P_true - P_taylor) / P_true), na.rm = TRUE)
  )
}


test_expm_taylor_on_hazards <- function(haz_df, dt = 1, n = 5) {
  
  build_Q <- function(h) {
    h_HU <- h["h_HU"]; h_HD <- h["h_HD"]
    h_UH <- h["h_UH"]; h_UD <- h["h_UD"]
    
    Q <- matrix(
      c(
        -(h_HU + h_HD), h_HU,            h_HD,
        h_UH,          -(h_UH + h_UD),   h_UD,
        0,              0,               0
      ),
      3, 3, byrow = TRUE
    )
    
    Q
  }
  
  purrr::pmap_dfr(
    list(
      sex  = haz_df$sex,
      age  = haz_df$age,
      hHU  = haz_df$h_HU,
      hHD  = haz_df$h_HD,
      hUH  = haz_df$h_UH,
      hUD  = haz_df$h_UD
    ),
    function(sex, age, hHU, hHD, hUH, hUD) {
      Q <- build_Q(c(h_HU = hHU, h_HD = hHD, h_UH = hUH, h_UD = hUD))
      res <- test_expm_taylor_single(Q, dt = dt, n = n)
      
      tibble::tibble(
        sex   = sex,
        age   = age,
        n     = n,
        max_abs_err = res$max_abs_err,
        max_rel_err = res$max_rel_err
      )
    }
  )
}

# res5  <- test_expm_taylor_on_hazards(haz_rowwise, dt = 1, n = 5)
# res10 <- test_expm_taylor_on_hazards(haz_rowwise, dt = 1, n = 10)
# ggplot(res10, aes(age, max_rel_err, color = sex)) +
#   geom_line() +
#   facet_wrap(~sex) +
#   scale_y_log10() +
#   theme_minimal() +
#   ggtitle("Absolute error of 10-term Taylor approx to expm(Q)")
# ggplot(res5, aes(age, max_abs_err, color = sex)) +
#   geom_line() +
#   facet_wrap(~sex) +
#   scale_y_log10() +
#   theme_minimal() +
#   ggtitle("Absolute error of 5-term Taylor approx to expm(Q)")

# assuming expm_taylor(A, t, n) already defined

hazards_to_probs_ctmc_taylor <- function(haz, dt = 1, n_terms = 4) {
  
  haz_mat <- as.data.frame(haz[, c("h_HU", "h_HD", "h_UH", "h_UD")])
  
  one_row_to_probs <- function(v) {
    h_HU <- v[1]; h_HD <- v[2]
    h_UH <- v[3]; h_UD <- v[4]
    
    Q <- matrix(
      c(
        -(h_HU + h_HD), h_HU,            h_HD,
        h_UH,          -(h_UH + h_UD),  h_UD,
        0,              0,              0
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
    haz %>% dplyr::select(sex, age),
    prob_df
  )
}

run_hazard_perturbation_ctmc_taylor <- function(haz_df, dt = 1, eps = 1e-4, n_terms = 4) {
  
  # Baseline probs via Taylor
  probs_orig <- hazards_to_probs_ctmc_taylor(haz_df, dt = dt, n_terms = n_terms)
  
  probs_orig_long <- probs_orig %>%
    tidyr::pivot_longer(
      cols      = dplyr::starts_with("p_"),
      names_to  = "transition",
      values_to = "p_orig"
    )
  
  hazards_to_perturb <- c("h_HU", "h_HD", "h_UH", "h_UD")
  
  out_list <- lapply(hazards_to_perturb, function(hname) {
    haz_pert <- perturb_one_hazard(haz_df, hname, eps = eps)
    
    probs_pert <- hazards_to_probs_ctmc_taylor(haz_pert, dt = dt, n_terms = n_terms)
    
    probs_pert_long <- probs_pert %>%
      tidyr::pivot_longer(
        cols      = dplyr::starts_with("p_"),
        names_to  = "transition",
        values_to = "p_pert"
      ) %>%
      dplyr::mutate(hazard_perturbed = hname)
    
    probs_orig_long %>%
      dplyr::inner_join(
        probs_pert_long,
        by = c("sex", "age", "transition")
      )
  })
  
  dplyr::bind_rows(out_list) %>%
    dplyr::mutate(
      dt      = dt,
      eps     = eps,
      n_terms = n_terms
    ) %>%
    dplyr::select(sex, age, hazard_perturbed, transition,
                  p_orig, p_pert, dt, eps, n_terms)
}
compare_expm_vs_taylor <- function(perturb_expm, perturb_taylor) {
  df <- dplyr::inner_join(
    perturb_expm,
    perturb_taylor,
    by = c("sex", "age", "hazard_perturbed", "transition", "dt", "eps"),
    suffix = c("_true", "_taylor")
  ) %>%
    dplyr::mutate(
      delta_true   = p_pert_true   - p_orig_true,
      delta_taylor = p_pert_taylor - p_orig_taylor,
      sign_true    = sign(delta_true),
      sign_taylor  = sign(delta_taylor),
      sign_match   = sign_true == sign_taylor,
      rel_err_delta = dplyr::if_else(
        delta_true != 0,
        (delta_taylor - delta_true) / delta_true,
        NA_real_
      )
    )
  
  df
}
summarise_sign_agreement <- function(comp_df) {
  comp_df %>%
    dplyr::group_by(hazard_perturbed) %>%
    dplyr::summarise(
      n = dplyr::n(),
      sign_match_frac = mean(sign_match, na.rm = TRUE),
      max_abs_rel_err = max(abs(rel_err_delta), na.rm = TRUE),
      median_abs_rel_err = median(abs(rel_err_delta), na.rm = TRUE),
      .groups = "drop"
    )
}
summarise_rank_origin <- function(comp_df) {
  comp_df %>%
    dplyr::mutate(
      origin_row = dplyr::case_when(
        hazard_perturbed %in% c("h_HU", "h_HD") ~ "H",
        hazard_perturbed %in% c("h_UH", "h_UD") ~ "U",
        TRUE ~ NA_character_
      ),
      in_origin_row = dplyr::case_when(
        origin_row == "H" & transition %in% c("p_HH","p_HU","p_HD") ~ TRUE,
        origin_row == "U" & transition %in% c("p_UH","p_UU","p_UD") ~ TRUE,
        TRUE ~ FALSE
      )
    ) %>%
    dplyr::group_by(sex, age, hazard_perturbed) %>%
    dplyr::summarise(
      # Spearman correlation of |Δp| within origin row
      spearman_abs_origin = {
        x <- abs(delta_true[in_origin_row])
        y <- abs(delta_taylor[in_origin_row])
        if (sum(is.finite(x) & is.finite(y)) >= 2) {
          suppressWarnings(stats::cor(x, y, method = "spearman", use = "complete.obs"))
        } else {
          NA_real_
        }
      },
      .groups = "drop"
    )
}
# 1) Get your CR-based hazards
haz_rowwise <- probs_to_hazards_cr(probs, dt = 1)

# 2) True expm perturbations
perturb_results <- run_hazard_perturbation_ctmc(
  haz_df = haz_rowwise,
  dt     = 1,
  eps    = 1e-3
)

# 3) Taylor perturbations (choose n_terms)
perturb_results_taylor <- run_hazard_perturbation_ctmc_taylor(
  haz_df  = haz_rowwise,
  dt      = 1,
  eps     = 1e-3,
  n_terms = 5   # or 3, 5, etc.
)

# 4) Compare
comp_true_vs_taylor <- compare_expm_vs_taylor(perturb_results, perturb_results_taylor)

# 5) Global sign & relative error summary
summarise_sign_agreement(comp_true_vs_taylor)

# 6) Rank ordering within origin row
rank_summary <- summarise_rank_origin(comp_true_vs_taylor) %>%
  dplyr::group_by(hazard_perturbed) %>%
  dplyr::summarise(
    median_spearman = median(spearman_abs_origin, na.rm = TRUE),
    min_spearman    = min(spearman_abs_origin, na.rm = TRUE),
    .groups = "drop"
  )
rank_summary

