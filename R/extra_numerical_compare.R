
# ------------------------------------------------ #
# This script compares all analytic sensitivities  # 
# with results of a numerical gradient function    #
# we get an exact-as-reasonably-possible match     #
# with all residuals attributable to numerical     #
# precision. This serves as a unit test of our     # 
# implementation of the analytic sensitivities     # 
# ------------------------------------------------ #

source("R/00_functions_classic.R")
source("R/00_sensitivity_functions.R")

trans <- read_csv("transitions_lievre2003_annual.csv")

# parameter case 1, HLE
s1nh <-
  trans |> 
  pivot_longer(-c(sex,age),
               names_to = "transition",
               values_to = "p") |> 
  group_by(sex) %>%
  group_modify(~s1nt(data = .x, expectancy = "h")) |> 
  ungroup() |> 
  mutate(expectancy = "h", 
         .before = 1)|> 
  mutate(case = 1, 
         .before = 1)
# parameter case 1, ULE
s1nu <-
  trans |> 
  pivot_longer(-c(sex,age),
               names_to = "transition",
               values_to = "p") |> 
  group_by(sex) %>%
  group_modify(~s1nt(data = .x, expectancy = "u")) |> 
  ungroup() |> 
  mutate(expectancy = "u", 
         .before = 1)|> 
  mutate(case = 1, 
         .before = 1)
# parameter case 1, LE
s1ntot <-
  trans |> 
  pivot_longer(-c(sex,age),
               names_to = "transition",
               values_to = "p") |> 
  group_by(sex) %>%
  group_modify(~s1nt(data = .x, expectancy = "t")) |> 
  ungroup() |> 
  mutate(expectancy = "t", 
         .before = 1)|> 
  mutate(case = 1, 
         .before = 1)
s1n_all <- bind_rows(s1nh, s1nu, s1ntot)
# --------------------------------------------------#
# parameter case 2, HLE
s2nh <-
  trans |> 
  pivot_longer(-c(sex,age),
               names_to = "transition",
               values_to = "p") |> 
  group_by(sex) %>%
  group_modify(~s2nt(data = .x, expectancy = "h")) |> 
  ungroup() |> 
  mutate(expectancy = "h", 
         .before = 1)|> 
  mutate(case = 2, 
         .before = 1)
# parameter case 2, ULE
s2nu <-
  trans |> 
  pivot_longer(-c(sex,age),
               names_to = "transition",
               values_to = "p") |> 
  group_by(sex) %>%
  group_modify(~s2nt(data = .x, expectancy = "u")) |> 
  ungroup() |> 
  mutate(expectancy = "u", 
         .before = 1)|> 
  mutate(case = 2, 
         .before = 1)
# parameter case 2, LE
s2ntot <-
  trans |> 
  pivot_longer(-c(sex,age),
               names_to = "transition",
               values_to = "p") |> 
  group_by(sex) %>%
  group_modify(~s2nt(data = .x, expectancy = "t")) |> 
  ungroup() |> 
  mutate(expectancy = "t", 
         .before = 1)|> 
  mutate(case = 2, 
         .before = 1)
s2n_all <- bind_rows(s2nh, s2nu, s2ntot)
# --------------------------------------------------#
# parameter case 3, HLE
s3nh <-
  trans |> 
  pivot_longer(-c(sex,age),
               names_to = "transition",
               values_to = "p") |> 
  group_by(sex) %>%
  group_modify(~s3nt(data = .x, expectancy = "h")) |> 
  ungroup() |> 
  mutate(expectancy = "h", 
         .before = 1)|> 
  mutate(case = 3, 
         .before = 1)
# parameter case 2, ULE
s3nu <-
  trans |> 
  pivot_longer(-c(sex,age),
               names_to = "transition",
               values_to = "p") |> 
  group_by(sex) %>%
  group_modify(~s3nt(data = .x, expectancy = "u")) |> 
  ungroup() |> 
  mutate(expectancy = "u", 
         .before = 1)|> 
  mutate(case = 3, 
         .before = 1)
# parameter case 2, LE
s3ntot <-
  trans |> 
  pivot_longer(-c(sex,age),
               names_to = "transition",
               values_to = "p") |> 
  group_by(sex) %>%
  group_modify(~s3nt(data = .x, expectancy = "t")) |> 
  ungroup() |> 
  mutate(expectancy = "t", 
         .before = 1)|> 
  mutate(case = 3, 
         .before = 1)
s3n_all <- bind_rows(s3nh, s3nu, s3ntot)

# ------------------------------------------------ #
sn_all <- bind_rows(s1n_all, s2n_all, s3n_all) |> 
  rename(sn = s)

s_all <- read_csv("all_sensitivities.csv")|> 
  rename(s = effect)

s_compare <- left_join(sn_all, s_all, 
                       by = join_by(sex,expectancy,transition,age,case))
# ------------------------------------------------ #
# all residuals are due to precision of numerical derivatives
# ?numDeriv::grad()
s_compare |> 
  mutate(diff = sn - s) |> 
  ggplot(aes(x = age, y = diff, color = transition, linetype = sex)) +
  geom_point() +
  facet_wrap(expectancy~case)

# ------------------------------------------------ #
# extra check for our shorthand solution for 
# sensitivity to initial conditions. This one is 
# satisfying.
s_compare |> 
  filter(transition == "init")
# ------------------------------------------------ #

# end

