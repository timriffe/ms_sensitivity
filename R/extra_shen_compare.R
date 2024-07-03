source("R/00_sensitivity_functions.R")
source("R/00_functions_classic.R")
# Compare with data and results from:

# Shen, T, Riffe, T, Payne, C, Canudas-Romo, V. 2023, 'Decomposition of Differentials in Health Expectancies From Multistate Life Tables: A Research Note', Demography 60 (6): 1675–1688. https://doi.org/10.1215/00703370-11058373

# Code repository:
# https://github.com/tyaSHEN/HLEdecom

# Transitions we reassess:
dat <- read_csv("https://raw.githubusercontent.com/tyaSHEN/HLEdecom/main/RMLE/PROB.csv",
                show_col_types = FALSE)
head(dat)
dat$ragender |> unique()

# Initial wrangling
p_med <- 
  dat |> 
  rename(H = A,
         U = L,
         D = H,
         state_from = pre_state) |> 
  pivot_longer(H:D, names_to = "state_to", values_to = "p") |> 
  mutate(state_from = if_else(state_from == "A","H","U")) |> 
  rename(sex = ragender) |> 
  mutate(sex = if_else(sex == 1, "m","f")) |> 
  group_by(state_from, state_to, age, sex) |> 
  summarize(p = median(p), .groups= "drop") |> 
  mutate(transition = paste0(state_from, state_to)) |> 
  select(-state_from, -state_to)

# Females have higher disability onset (HU), so we expect this to contribute negatively to the HLE gap
p_med |> 
  filter(transition == "HU") |> 
  ggplot(aes(x = age, y = p, linetype = sex)) +
  geom_line() +
  theme_minimal() +
  labs(title = "transition from good to poor health, Phu",
       y = "probability")

# here are the deltas used in the decomposition, 
# focus on the green line for HU
deltas <- p_med |> 
  pivot_wider(names_from = sex, values_from = p) |> 
  mutate(diff = f - m) 

deltas |> 
  ggplot(aes(x= age, y = diff, color = transition)) +
  geom_line() +
  theme_minimal()

# calculate P1 and P2 sensitivities from the Shen data
p_med_s1 <-
  p_med |> 
  group_by(sex) |> 
  group_modify(~s1t(data = .x, expectancy = "all"))|> 
  mutate(case = 1, .before = 1) |> 
  mutate(age = age + 55)

p_med_s2 <-
  p_med |> 
  group_by(sex) |> 
  group_modify(~s2t(data = .x, expectancy = "all")) |> 
  mutate(case = 2, .before = 1)|> 
  mutate(age = age + 55)

s_compare <- 
  bind_rows(p_med_s1,p_med_s2)

# Here we see 
# HU has a positive sensitivity for P1 and negative for P2 (green lines)
# Meaning HU up = HLE up for P1, whereas HU up = HLE down in P2.

# I propose an axiom:
# "Leaving a state cannot contribute positively to time spent in the 
# state being left."
s_compare |> 
  filter(expectancy == "h",
         transition != "init") |> 
  mutate(case = paste0("P",case)) |> 
  ggplot(aes(x = age, y = effect, color = transition, linetype = sex)) +
  geom_line() +
  facet_wrap(~case) +
  theme_minimal()

# Now the decomposition result we can compare directly,
# since the Shen paper has a figure for this.

# We need to recalculate the sensitivities at the parameter average
p_avg_s1 <-
  p_med |> 
  group_by(transition, age) |> 
  summarize(p = mean(p), .groups = "drop") |> 
  group_modify(~s1t(data = .x, expectancy = "all"))|> 
  mutate(case = 1, .before = 1) |> 
  mutate(age = age + 55)

p_avg_s2 <-
  p_med |> 
  group_by(transition, age) |> 
  summarize(p = mean(p), .groups = "drop") |> 
  group_modify(~s2t(data = .x, expectancy = "all")) |> 
  mutate(case = 2, .before = 1)|> 
  mutate(age = age + 55)

# combine


# final calc for decomposition
d_compare <-
  bind_rows(p_avg_s1,p_avg_s2) |> 
  filter(expectancy == "h",
         transition != "init") |>
  left_join(deltas |> 
              select(age, transition, delta = diff),
            by = join_by(age, transition)) |> 
  mutate(cc = effect * delta,
         case = paste0("P",case)) 
  
# Compare P1 (left side of below figure) to Fig 3 from 
# Shen, T, Riffe, T, Payne, C, Canudas-Romo, V. 2023, 'Decomposition of Differentials in Health Expectancies From Multistate Life Tables: A Research Note', Demography 60 (6): 1675–1688. https://doi.org/10.1215/00703370-11058373

# Here HU is green, and in Shen paper HU is also green (solid line)
d_compare |> 
  ggplot(aes(x = age, y = cc, color = transition)) +
  geom_line() +
  facet_wrap(~case) +
  theme_minimal() +
  ylim(-.05,.075)
# This verifies that our P1 maps to Shen's decomposition.
# Note also that HH of P1 is not the same as HD of P2,
# UH values are different, and P1 UU is different from
# P2 UD. 
# Explanation: Sensitivities for P1 HH and P2 HD are same 
# magnitude,opposite sign, but the deltas are different.

# end

