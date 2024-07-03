# This script generates transition probabilities from coefficients given in:

# Lièvre, Agnès, Nicolas Brouard, and Christopher Heathcote. 
# "The estimation of health #expectancies from cross-longitudinal surveys." 
# Mathematical population studies 10.4 (2003): #211-248.

# This script generates annual transitions used elsewhere in the code 
# chain
library(tidyverse)
library(compositions)
library(expm)

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

# a custom function to convert month-scale to annualized transition
# probabilities
pmonth_to_pannual <- function(x){
  x %>% 
    pivot_longer(-1,
                 names_to = "from_to", 
                 values_to = "p") %>% 
    mutate(from = substr(from_to,1,1), 
           to = substr(from_to,2,2)) %>% 
    select(-age,-from_to) %>% 
    add_row(from = "D", to = "D", p = 1) %>% 
    pivot_wider(names_from = to, 
                values_from = p, 
                values_fill = 0) %>% 
    select(from, H, U, D) %>% 
    column_to_rownames("from") %>% 
    as.matrix() %>% 
    `%^%`(12) %>% 
    as.data.frame() %>% 
    rownames_to_column("from") %>% 
    pivot_longer(-from, 
                 names_to = "to",
                 values_to = "p") %>% 
    filter(from != "D") %>% 
    unite(from_to, from, to,sep="") %>% 
    pivot_wider(names_from = 
                  from_to, 
                values_from = p)
}

# get monthly transitions, one sex and one origin state at a time.
# eq33, which returns ALR transforms of HU and HD (UH and UD) such that
# we can retrieve all six transitions using the alrInv function.
# these will be monthly transitions.
ages <- 50:110

# for women
fromHf <- 
  cbind(HU = eq33(tab2, "HU", ages, 1),
        HD = eq33(tab2, "HD", ages, 1)) %>% 
  alrInv() %>% 
  as_tibble() %>% 
  rename(HH = V3) %>% 
  mutate(age = ages) 

fromUf <- 
  cbind(UH = eq33(tab2, "UH", ages, 1),
        UD = eq33(tab2, "UD", ages, 1)) %>% 
  alrInv() %>% 
  as_tibble() %>% 
  rename(UU = V3) %>% 
  mutate(age = ages)

fromHm <- 
  cbind(HU = eq33(tab2, "HU", ages, 0),
        HD = eq33(tab2, "HD", ages, 0)) %>% 
  alrInv() %>% 
  as_tibble() %>% 
  rename(HH = V3) %>% 
  mutate(age = ages) 

fromUm <- 
  cbind(UH = eq33(tab2, "UH", ages, 0),
        UD = eq33(tab2, "UD", ages, 0)) %>% 
  alrInv() %>% 
  as_tibble() %>% 
  rename(UU = V3) %>% 
  mutate(age = ages)

left_join(fromHf,
          fromUf,
          by="age") %>% 
  pivot_longer(-age,
               names_to = "from_to", 
               values_to = "p") %>% 
  ggplot(aes(x = age, y = p, color = from_to)) +
  geom_line() +
  theme_minimal() +
  labs(title = "monthly transitions from Lievre 2003")

p_tibble_orig_monthly_f <- left_join(fromHf,
                                     fromUf,
                                     by = "age") %>% 
  relocate(age, .before = 1)
p_tibble_orig_monthly_m <- left_join(fromHm,
                                     fromUm,
                                     by = "age") %>% 
  relocate(age,.before=1)

# convert to annualized transitions. We do this one age at a time, then bind.
p_tibble_orig_annual_f <-
  p_tibble_orig_monthly_f %>% 
  split(p_tibble_orig_monthly_f$age) %>% 
  lapply(pmonth_to_pannual) %>% 
  bind_rows() %>% 
  mutate(age = ages, .before = 1)  %>% 
  mutate(sex = "f", .before = age) 

p_tibble_orig_annual_m <-
  p_tibble_orig_monthly_m %>% 
  split(p_tibble_orig_monthly_m$age) %>% 
  lapply(pmonth_to_pannual) %>% 
  bind_rows() %>% 
  mutate(age = ages, .before = 1) %>% 
  mutate(sex = "m", .before =age) 

# These are the transitions we use in the paper.
trans <- 
  bind_rows(p_tibble_orig_annual_m,
            p_tibble_orig_annual_f) 
trans |> 
  write_csv("transitions_lievre2003_annual.csv")













