

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
    L <- expm::logm(P, method = "Eigen")
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



f2dd <- function(probs_vec, df, init = c(H=1,U=0), expectancy = "h",dt=1){
  
  df$p <- probs_vec
  df   <- pivot_wider(df, names_from = transition, values_from = p)
  f2(hd = df$HD,
     hu = df$HU,
     ud = df$UD,
     uh = df$UH,
     init = init,
     expectancy = expectancy,
     interval = dt)
}
# Now we need the same decomp but using probs_cr
f2dd_haz_ctmc <- function(haz_vec,
                          df,
                          init = c(H = 1, U = 0),
                          expectancy = "h",
                          dt = 1) {
  
  # Attach hazards to long template
  df$h <- haz_vec
  
  # Wide hazards per age, with h_ prefix added ONCE
  haz_wide <- df %>%
    tidyr::pivot_wider(names_from = transition, values_from = h) %>%
    dplyr::arrange(age) %>%
    dplyr::mutate(sex = "x") %>%  # dummy
    dplyr::select(sex, age, h_HU, h_HD, h_UH, h_UD)
  
  # Hazards → probabilities (names now p_*)
  probs <- hazards_to_probs_ctmc(haz_wide, dt = dt)
  
  # Feed probabilities to P2 calculator
  f2(
    hd = probs$p_HD,
    hu = probs$p_HU,
    ud = probs$p_UD,
    uh = probs$p_UH,
    init = init,
    expectancy = expectancy,
    interval = dt
  )
}







