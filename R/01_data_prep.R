# This script generates transition probabilities from coefficients given in:
# Lièvre, Agnès, Nicolas Brouard, and Christopher Heathcote.
# "The estimation of health expectancies from cross-longitudinal surveys."
# Mathematical Population Studies 10.4 (2003): 211-248.

library(tidyverse)
library(compositions)
library(expm)

tab2 <- tibble(
  from_to  = c("HU","HD","UH","UD"),
  intercept= c(-11.9,-11.9,-1.4,-5.4),
  age      = c(0.0853,0.0812,-0.0358,0.0207),
  sex      = c(-2.4,3,-.2,-2.5),
  age_sex  = c(0.0316,-0.0472,0.0010,0.0249)
)

# Eq 33 in Lièvre et al. (2003): ALR for two transitions with the third as denominator
eq33 <- function(tab2, .from_to = "HU", ages = 70:90, sex = 1){
  coefs <- tab2 |> filter(from_to == .from_to)
  coefs$intercept + ages * coefs$age + sex * coefs$sex + ages * sex * coefs$age_sex
}

# --- helper: build a 3x3 P matrix from a single row with HH,HU,HD,UH,UU,UD
row_to_P <- function(row) {
  matrix(
    c(
      row[["HH"]], row[["HU"]], row[["HD"]],
      row[["UH"]], row[["UU"]], row[["UD"]],
      0,           0,           1
    ),
    nrow = 3, byrow = TRUE,
    dimnames = list(c("H","U","D"), c("H","U","D"))
  )
}

# --- annualize by product of 12 month matrices inside each integer age
# expects df with columns age (monthly grid), HU,HD,HH,UH,UD,UU
annualize_monthly_probs <- function(df_month, months_per_year = 12) {
  
  stopifnot("age" %in% names(df_month))
  
  df2 <- df_month %>%
    mutate(age_year = floor(age)) %>%
    arrange(age)
  
  # sanity check: do we actually have 12 months per age-year?
  counts <- df2 %>% count(age_year, name = "n")
  if (any(counts$n != months_per_year)) {
    bad <- counts %>% filter(n != months_per_year)
    stop(
      "annualize_monthly_probs(): expected ", months_per_year,
      " rows per floor(age), but got:\n",
      paste0(capture.output(print(bad, n = 50)), collapse = "\n"),
      "\n\nThis can happen if rounding causes some months to duplicate or drop.\n",
      "Consider generating age as 50 + (0:...)/12 exactly (not rounded)."
    )
  }
  
  df2 %>%
    group_by(age_year) %>%
    group_modify(~{
      d <- .x %>% arrange(age)
      
      # Build monthly P matrices (H,U transient; D absorbing)
      P_list <- lapply(seq_len(nrow(d)), function(i) {
        matrix(
          c(
            d$HH[i], d$HU[i], d$HD[i],
            d$UH[i], d$UU[i], d$UD[i],
            0,       0,       1
          ),
          nrow = 3, byrow = TRUE,
          dimnames = list(c("H","U","D"), c("H","U","D"))
        )
      })
      
      # Chronological product across the 12 months
      P_year <- Reduce(`%*%`, P_list)
      
      tibble::tibble(
        age = .y$age_year,   # <-- group key lives here!
        HH  = P_year["H","H"],
        HU  = P_year["H","U"],
        HD  = P_year["H","D"],
        UH  = P_year["U","H"],
        UU  = P_year["U","U"],
        UD  = P_year["U","D"]
      )
    }) %>%
    ungroup()
}


# -----------------------------
# Monthly age grid (not just January-of-age-x)
ages_y <- 50:110
ages_m <- seq(min(ages_y), max(ages_y) + 11/12, by = 1/12)

# --- women monthly
fromHf_m <- cbind(
  HU = eq33(tab2, "HU", ages_m, 1),
  HD = eq33(tab2, "HD", ages_m, 1)
) %>%
  alrInv() %>%
  as_tibble() %>%
  rename(HH = V3) %>%
  mutate(age = ages_m)

fromUf_m <- cbind(
  UH = eq33(tab2, "UH", ages_m, 1),
  UD = eq33(tab2, "UD", ages_m, 1)
) %>%
  alrInv() %>%
  as_tibble() %>%
  rename(UU = V3) %>%
  mutate(age = ages_m)

p_tibble_orig_monthly_f <- left_join(fromHf_m, fromUf_m, by = "age") %>%
  relocate(age, .before = 1)

# --- men monthly
fromHm_m <- cbind(
  HU = eq33(tab2, "HU", ages_m, 0),
  HD = eq33(tab2, "HD", ages_m, 0)
) %>%
  alrInv() %>%
  as_tibble() %>%
  rename(HH = V3) %>%
  mutate(age = ages_m)

fromUm_m <- cbind(
  UH = eq33(tab2, "UH", ages_m, 0),
  UD = eq33(tab2, "UD", ages_m, 0)
) %>%
  alrInv() %>%
  as_tibble() %>%
  rename(UU = V3) %>%
  mutate(age = ages_m)

p_tibble_orig_monthly_m <- left_join(fromHm_m, fromUm_m, by = "age") %>%
  relocate(age, .before = 1)

# -----------------------------
# Annualize: product of 12 month-specific matrices within each age-year
p_tibble_orig_annual_f <- annualize_monthly_probs(p_tibble_orig_monthly_f) %>%
  mutate(sex = "f", .before = age)

p_tibble_orig_annual_m <- annualize_monthly_probs(p_tibble_orig_monthly_m) %>%
  mutate(sex = "m", .before = age)

# save out annual
bind_rows(p_tibble_orig_annual_m, p_tibble_orig_annual_f)|> 
  select(-age_year) |> 
  write_csv("transitions_lievre2003_annual.csv")

# save out monthly
bind_rows(
  p_tibble_orig_monthly_m |> 
    mutate(sex = "m", .before = age),
  
  p_tibble_orig_monthly_f |> 
    mutate(sex = "f", .before = age)) |> 
  write_csv("transitions_lievre2003_monthly.csv")
