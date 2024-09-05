# Script performs all sensitivity and decomposition calculations

library(tidyverse)
library(compositions)
library(expm)
source("R/00_sensitivity_functions.R")
source("R/00_functions_classic.R")

# this csv is created from scratch in the previous script 01_data_prep.R
trans <- read_csv("transitions_lievre2003_annual.csv")

# For Table 1: Expectancies and differences
expectancies <-
  trans |> 
  pivot_longer(-c(sex,age),
               names_to = "transition",
               values_to = "p") |> 
  group_by(sex) %>%
  summarize(hle = f1t(data = pick(everything()), 
                      expectancy = "h"),
            ule = f1t(data = pick(everything()), 
                      expectancy = "u"),
            le = f1t(data = pick(everything()), 
                     expectancy = "t"))

# ------------------------------------------- #
# Check equivalency between parameter cases,
# just eyeball output in the console
expectancies # case 1

# case 2
trans |> 
  pivot_longer(-c(sex,age),
               names_to = "transition",
               values_to = "p") |> 
  group_by(sex) %>%
  summarize(hle = f2t(data=pick(everything()),
                      expectancy = "h"),
            ule = f2t(data=pick(everything()),
                      expectancy = "u"),
            le = f2t(data=pick(everything()),
                     expectancy = "t"))
# case 3
trans |> 
  pivot_longer(-c(sex,age),
               names_to = "transition",
               values_to = "p") |> 
  group_by(sex) %>%
  summarize(hle = f3t(data=pick(everything()),
                      expectancy = "h"),
            ule = f3t(data=pick(everything()),
                      expectancy = "u"),
            le = f3t(data=pick(everything()),
                     expectancy = "t"))
# ---------------------------------------------- #
# calculate the three sensitivities and bind for #
# comparison                                     #
# ---------------------------------------------- #

# parameter case 1
s1all <-
  trans |> 
  pivot_longer(-c(sex,age),
               names_to = "transition",
               values_to = "p") |> 
  group_by(sex) %>%
  group_modify(~s1t(data = .x, expectancy = "all", interval = 1)) |> 
  ungroup() |> 
  mutate(case = 1, 
         .before = 1)
# parameter case 2
s2all <-
  trans |> 
  pivot_longer(-c(sex,age),
               names_to = "transition",
               values_to = "p") |> 
  group_by(sex) %>%
  group_modify(~s2t(data = .x,expectancy = "all")) |> 
  ungroup() |> 
  mutate(case = 2, 
         .before = 1)

# parameter case 3
s3all <-
  trans |> 
  pivot_longer(-c(sex,age),names_to = "transition",values_to = "p") |> 
  group_by(sex) %>%
  group_modify(~s3t(data = .x,expectancy = "all")) |> 
  ungroup() |> 
  mutate(case = 3, 
         .before = 1)

# join
sen_all <-
  bind_rows(s1all,
            s2all,
            s3all)

# we take the arithmetic average of male and female transitions for main
# manuscript results. Here we also have delta
trans_avg <-
  trans |> 
  pivot_longer(HH:UD, names_to = "transition", values_to = "p") |> 
  pivot_wider(names_from = sex, values_from = p) |> 
  mutate(p = (m + f) / 2,
         delta = f - m) |> 
  select(-f, -m)

delta <- 
  trans_avg |> 
  select(-p)

# rather laborious setup for initial conditions. This is laborious because
# we generate the initial conditions using the transitions in the first time
# step, but in the decomposition we would like to treat the initial conditions
# as given, since most often one will use prevalence or something else.
df_init <-
  trans |> 
  filter(age == 50) |> 
  group_by(sex) 

init_m <- df_init |> 
  filter(sex == "m") %>% 
  '['(,c("HH","UH","HU","UU")) |> 
  c() |> 
  init_constant()

init_f <- df_init |> 
  filter(sex == "f") %>% 
  '['(,c("HH","UH","HU","UU")) |> 
  c() |> 
  init_constant()

# To compute delta for init, use difference in init["H"]
init_delta <- c(init_f - init_m)[1] 

# bind init delta to main delta
delta <-
  tibble(age = 50, 
         transition = "init", 
         delta = init_delta) |> 
  bind_rows(delta)

# We now have sensitivities and deltas, and can perform the decompositions

# parameter case 1
d1all <-
  trans_avg |> 
  group_modify(~s1t(data = .x,expectancy = "all")) |> 
  ungroup() |> 
  mutate(case = 1, 
         .before = 1) |> 
  rename(s = effect) |> 
  mutate(age = age + 50) |> 
  left_join(delta, by = join_by(age, transition)) |> 
  mutate(cc = s * delta) |> 
  filter(age < 111)


# parameter case 2
d2all <-
  trans_avg |> 
  group_modify(~s2t(data = .x,expectancy = "all")) |> 
  ungroup() |> 
  mutate(case = 2, 
         .before = 1) |> 
  rename(s = effect) |> 
  mutate(age = age + 50) |> 
  left_join(delta, by = join_by(age, transition)) |> 
  mutate(cc = s * delta) |> 
  filter(age < 111)

# parameter case 3
d3all <-
  trans_avg |> 
  group_modify(~s3t(data = .x,expectancy = "all")) |> 
  ungroup() |> 
  mutate(case = 3, 
         .before = 1) |> 
  rename(s = effect) |> 
  mutate(age = age + 50) |> 
  left_join(delta, by = join_by(age, transition)) |> 
  mutate(cc = s * delta) |> 
  filter(age < 111)

# combine all decomposition results
d_all <-
  bind_rows(d1all,
            d2all,
            d3all)

# save out sensitivities to filename referenced in manuscript
write_csv(sen_all,"all_sensitivities.csv")
# save decomposition results to file named in manuscript
write_csv(d_all, file = "all_decompositions.csv")




