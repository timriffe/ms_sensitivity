library(expm)
library(dplyr)
library(purrr)
library(tidyr)
library(readr)
library(expm)

## Build P from scalar probabilities (one row)
build_P_from_scalars <- function(p_HH, p_HU, p_HD,
                                 p_UH, p_UU, p_UD) {
  matrix(
    c(
      p_HH, p_HU, p_HD,
      p_UH, p_UU, p_UD,
      0,    0,    1
    ),
    nrow = 3,
    byrow = TRUE
  )
}

## Matrix log: P -> Q
mlog_markov <- function(P, dt = 1) {
  P <- as.matrix(P)
  # enforce numeric 3x3
  if (length(P) != 9) {
    stop("P does not have 9 elements; got length = ", length(P))
  }
  dim(P) <- c(3, 3)
  
  L <- expm::logm(P)
  
  # Drop imaginary parts (numerical noise)
  if (any(abs(Im(L)) > 1e-8)) {
    warning("Non-negligible imaginary parts in logm(P); check P.")
  }
  Q <- Re(L) / dt
  
  # Enforce row sums 0 by fixing diagonals
  rs <- rowSums(Q)
  diag(Q) <- diag(Q) - rs
  
  Q
}

## Fraction of negative relative change absorbed by p_ii
fraction_absorbed_by_staying <- function(Q, i, j, eps = 1e-4, dt = 1) {
  Q <- as.matrix(Q)
  if (!all(dim(Q) == c(3, 3))) {
    stop("Q is not 3x3; dim(Q) = ", paste(dim(Q), collapse = "x"))
  }
  if (i == j) stop("i and j must differ")
  
  # Baseline P
  P     <- expm::expm(Q * dt)
  p_row <- P[i, ]
  
  h_ij <- Q[i, j]
  if (abs(h_ij) < 1e-12) return(NA_real_)  # essentially zero hazard
  
  # Perturb hazard h_ij and adjust only diagonal to keep row sum 0
  Qp <- Q
  delta_h  <- h_ij * eps
  Qp[i, j] <- h_ij + delta_h
  Qp[i, i] <- Q[i, i] - delta_h
  
  # New P
  Pp     <- expm::expm(Qp * dt)
  p_rowp <- Pp[i, ]
  
  # Relative changes
  rel_change <- (p_rowp - p_row) / p_row
  
  # Negative components are the "counterperturbation"
  neg_idx <- which(rel_change < 0 & is.finite(rel_change))
  if (length(neg_idx) == 0) return(0)    # nothing shrank
  
  # If p_ii didn't shrink, it absorbed none of the counterperturbation
  if (!(i %in% neg_idx)) return(0)
  
  numer <- abs(rel_change[i])                  # |Δrel p_ii|
  denom <- sum(abs(rel_change[neg_idx]))       # total |Δrel| of shrinking probs
  
  if (denom == 0) return(0)
  
  numer / denom    # guaranteed in [0,1]
}
## Compute all four fractions for ONE row of probs_df
compute_fracs_for_row <- function(row_vec, dt = 1, eps = 1e-4) {
  # row_vec: numeric with elements in order:
  # p_HH, p_HU, p_HD, p_UH, p_UU, p_UD
  p_HH <- row_vec[1]
  p_HU <- row_vec[2]
  p_HD <- row_vec[3]
  p_UH <- row_vec[4]
  p_UU <- row_vec[5]
  p_UD <- row_vec[6]
  
  P <- build_P_from_scalars(p_HH, p_HU, p_HD, p_UH, p_UU, p_UD)
  Q <- mlog_markov(P, dt = dt)
  
  frac_p_HH_h_HU <- fraction_absorbed_by_staying(Q, i = 1, j = 2, eps = eps, dt = dt)
  frac_p_HH_h_HD <- fraction_absorbed_by_staying(Q, i = 1, j = 3, eps = eps, dt = dt)
  
  frac_p_UU_h_UH <- fraction_absorbed_by_staying(Q, i = 2, j = 1, eps = eps, dt = dt)
  frac_p_UU_h_UD <- fraction_absorbed_by_staying(Q, i = 2, j = 3, eps = eps, dt = dt)
  
  c(
    frac_p_HH_h_HU = frac_p_HH_h_HU,
    frac_p_HH_h_HD = frac_p_HH_h_HD,
    frac_p_UU_h_UH = frac_p_UU_h_UH,
    frac_p_UU_h_UD = frac_p_UU_h_UD
  )
}

probs_df <- read_csv("transitions_lievre2003_annual.csv") |> 
  pivot_longer(-c(sex,age), names_to = "transition",values_to = "p") |> 
  mutate(transition = paste0("p_",transition)) |> 
  pivot_wider(names_from = transition, values_from = p) 
  
  # Extract just the probs as a plain data.frame
prob_mat <- as.data.frame(probs_df[, c("p_HH", "p_HU", "p_HD", "p_UH", "p_UU", "p_UD")])

# Apply per row; result is 122 x 4 matrix
fracs_mat <- t(apply(prob_mat, 1, compute_fracs_for_row, dt = 1, eps = 1e-4))

# Convert to data.frame with proper names
fracs_df <- as.data.frame(fracs_mat)

# Bind back to original tibble
probs_df_with_fracs <- dplyr::bind_cols(probs_df, fracs_df)

# Peek
head(probs_df_with_fracs)
probs_df_with_fracs %>%
  filter(sex == "f") %>%
  select(age, frac_p_HH_h_HD, frac_p_HH_h_HU,frac_p_UU_h_UH,frac_p_UU_h_UD) %>%
  tidyr::pivot_longer(-age, names_to = "which", values_to = "frac") %>%
  ggplot(aes(age, frac, colour = which)) +
  geom_line() +
  labs(y = "Fraction of counterperturbation absorbed by staying")




library(expm)

# already defined earlier, but for completeness:
build_P_from_scalars <- function(p_HH, p_HU, p_HD,
                                 p_UH, p_UU, p_UD) {
  matrix(
    c(
      p_HH, p_HU, p_HD,
      p_UH, p_UU, p_UD,
      0,    0,    1
    ),
    nrow = 3,
    byrow = TRUE
  )
}

mlog_markov <- function(P, dt = 1) {
  P <- as.matrix(P)
  if (length(P) != 9) stop("P must have 9 elements")
  dim(P) <- c(3, 3)
  
  L <- expm::logm(P)
  if (any(abs(Im(L)) > 1e-8)) {
    warning("Non-negligible imaginary parts in logm(P); check P.")
  }
  Q <- Re(L) / dt
  
  rs <- rowSums(Q)
  diag(Q) <- diag(Q) - rs
  
  Q
}

# NEW: hazards from a single row of probs_df
compute_hazards_for_row <- function(row_vec, dt = 1) {
  # row_vec: numeric c(p_HH, p_HU, p_HD, p_UH, p_UU, p_UD)
  p_HH <- row_vec[1]
  p_HU <- row_vec[2]
  p_HD <- row_vec[3]
  p_UH <- row_vec[4]
  p_UU <- row_vec[5]
  p_UD <- row_vec[6]
  
  P <- build_P_from_scalars(p_HH, p_HU, p_HD, p_UH, p_UU, p_UD)
  Q <- mlog_markov(P, dt = dt)
  
  c(
    h_HU = Q[1, 2],
    h_HD = Q[1, 3],
    h_UH = Q[2, 1],
    h_UD = Q[2, 3]
  )
}
haz_mat <- t(apply(prob_mat, 1, compute_hazards_for_row, dt = 1))
haz_df  <- as.data.frame(haz_mat)

# bind back
probs_df_with_haz <- bind_cols(probs_df, haz_df)


probs_df_with_haz %>%
  select(sex, age, h_HU, h_HD, h_UH, h_UD) %>%
  pivot_longer(-c(sex,age), names_to = "hazard", values_to = "value") %>%
  ggplot(aes(x = age, y = value, colour = hazard)) +
  geom_line() +
  labs(
    x = "Age",
    y = "Hazard (per age)",
    colour = "Transition",
    title = "Implied CTMC hazards by age"
  ) +
  theme_minimal() +
  facet_wrap(~sex)


probs_df |> 
  pivot_longer(-c(sex,age), names_to = "transition", values_to = "p") |> 
  filter(!transition %in% c("p_HH","p_UU")) |> 
  ggplot(aes(x=age,y=p,color=transition)) +
  geom_line() +
  facet_wrap(~sex) +
  theme_minimal() +
  labs(
    x = "Age",
    y = "Probability (per age)",
    colour = "transition",
    title = "Original probabilities by age"
  ) 


P <- matrix(c(.8,.1,.1,0,0,0,0,0,0),3,byrow=TRUE)
P
logm(P)
