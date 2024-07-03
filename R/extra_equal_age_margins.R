# This script demonstrates that age margins are all equal, 
# no matter which parameterization we choose

library(tidyverse)
dec <- read_csv("all_decompositions.csv")

dec |> 
  filter(transition != "init") |> 
  group_by(case, expectancy, age) |> 
  summarize(cc = sum(cc)) |> 
  mutate(case = paste0("P",case)) |> 
  ggplot(aes(x = age, y = cc)) +
  geom_line() +
  facet_wrap(expectancy~case) +
  theme_minimal() +
  labs(title = "Each parameterization (row) is equal here",
       subtitle = "Differences appear when we want to know how these margins partition over transitions")
