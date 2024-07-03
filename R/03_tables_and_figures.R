source("R/02_calculate.R")
library(xtable)

# -------------------------------------------------#
# make Figure 2, showing all transitions by sex

fig2 <- 
  trans |> 
  pivot_longer(-c(sex,age),names_to="transition",values_to = "probability") |> 
  ggplot(aes(x=age,y=probability,color=transition, linetype = sex))+
  geom_line() +
  theme_minimal(base_size = 22) +
  scale_color_discrete(labels=c(HH=expression(p[hh]),
                                HU=expression(p[hu]),
                                HD=expression(p[hd]),
                                UU=expression(p[uu]),
                                UH=expression(p[uh]),
                                UD=expression(p[ud]) )) 
fig2
ggsave("fig_transitions.pdf",fig2)

# Table 1: Expectancies and differences, LaTeX code later modified somewhat
expectancies |> 
  pivot_longer(hle:le, names_to = "expectancy", values_to = "value") |> 
  pivot_wider(names_from = sex, values_from = value) |> 
  mutate(diff = f - m) |> 
  pivot_longer(f:diff, names_to = "sex", values_to = "value") |> 
  pivot_wider(names_from = expectancy, values_from = value) |> 
  xtable()

fig3 <-
   sen_all |> 
   filter(transition != "init",
          expectancy != "t") |> 
   mutate(case = paste0("P",case),
          age=age+50) |> 
   rename(sensitivity = effect) |> 
   ggplot(aes(x=age,y=sensitivity,color=transition,linetype=sex)) +
   geom_line() +
   facet_grid(rows = vars(expectancy),
              cols = vars(case),
              scales="free_y") +
   scale_color_discrete(labels=c(HH=expression(p[hh]),
                                 HU=expression(p[hu]),
                                 HD=expression(p[hd]),
                                 UU=expression(p[uu]),
                                 UH=expression(p[uh]),
                                 UD=expression(p[ud]) )) +
   theme_minimal() +
   coord_cartesian(clip = 'off') 
fig3
ggsave("fig_sensitivities.pdf", fig3, width = 5, height = 8)
# Appendix Figure of sensitivities for females, all three cases
fig_sen_appendix <-
  sen_all |> 
  filter(transition != "init") |> 
  mutate(case = paste0("P",case),
         age=age+50,
         expectancy = case_when(
           expectancy == "h" ~ "DFLE",
           expectancy == "u" ~ "DLE",
           expectancy == "t" ~ "LE"
         )) |> 
  rename(sensitivity = effect) |> 
  ggplot(aes(x=age,y=sensitivity,color=transition,linetype=sex)) +
  geom_line() +
  facet_grid(rows = vars(expectancy),
             cols = vars(case),
             scales="free_y") +
  scale_color_discrete(labels=c(HH=expression(p[hh]),
                                HU=expression(p[hu]),
                                HD=expression(p[hd]),
                                UU=expression(p[uu]),
                                UH=expression(p[uh]),
                                UD=expression(p[ud]) )) +
  theme_minimal() +
  coord_cartesian(clip = 'off') 


fig_sen_appendix
ggsave("fig_appendix_sensitivities.pdf", fig_sen_appendix, width = 8, height = 7)



# Appendix table initial condition sensitivities
sen_all |> 
  filter(case == 1,
         transition == "init") |> 
  mutate(expectancy = case_when(
    expectancy == "h" ~ "DFLE",
    expectancy == "u" ~ "DLE",
    TRUE ~ "LE"
  )) |> 
  pivot_wider(names_from = expectancy, values_from = effect) |> 
  select(sex, DFLE, DLE, LE) |> 
  xtable()



# Appendix A Figure showing transition probability differences 
# (deltas used for decomp)
f_deltas <-
  trans_avg |> 
  filter(transition != "init") |> 
  ggplot(aes(x = age, y = delta, color = transition)) +
  geom_line() +
  theme_minimal() +
  labs(y = "difference in probability (female - male)") +
  scale_color_discrete(labels=c(HH = expression(p[hh]),
                                HU = expression(p[hu]),
                                HD = expression(p[hd]),
                                UU = expression(p[uu]),
                                UH = expression(p[uh]),
                                UD = expression(p[ud]))) 
ggsave("fig_deltas_appendix.pdf",f_deltas, width=5,height=4)

# Checking crossover age for HU deltas, which is considered in manuscript
# discussion
# trans_avg |> 
#   filter(transition == "HU") |> 
#   select(age, delta) |> 
#   View()

# Manuscript Figure 4,
# shows decomposition results for DFLE (=HLE) only.
fig4 <- 
  d_all |> 
  filter(transition != "init",
         expectancy == "h") |> 
  mutate(case = paste0("P", case)) |> 
  ggplot(aes(x = age, y = cc, color = transition)) +
  geom_line() +
  facet_wrap(~case) +
  theme_minimal() +
  labs(y = "contribution to difference in HLE") + 
  scale_color_discrete(labels=c(HH = expression(p[hh]),
                                HU = expression(p[hu]),
                                HD = expression(p[hd]),
                                UU = expression(p[uu]),
                                UH = expression(p[uh]),
                                UD = expression(p[ud]))) +
  theme_minimal() +
  coord_cartesian(clip = 'off') 
fig4
ggsave("fig_decomp.pdf",fig4, width=5,height=4)

# Appendix C Figure decomp: all expectancies considered
f_d_appendix <-
  d_all |> 
  filter(transition != "init") |> 
  mutate(case = paste0("P",case),
         age=age,
         expectancy = case_when(
           expectancy == "h" ~ "DFLE",
           expectancy == "u" ~ "DLE",
           expectancy == "t" ~ "LE"
         )) |> 
    rename(effect = cc) |> 
  ggplot(aes(x=age,y=effect,color=transition)) +
  geom_line() +
  facet_grid(rows = vars(expectancy),
             cols = vars(case),
             scales="free_y") +
  scale_color_discrete(labels=c(HH=expression(p[hh]),
                                HU=expression(p[hu]),
                                HD=expression(p[hd]),
                                UU=expression(p[uu]),
                                UH=expression(p[uh]),
                                UD=expression(p[ud]) )) +
  theme_minimal() +
  coord_cartesian(clip = 'off') 
f_d_appendix
ggsave("fig_decomp_appendix.pdf",f_d_appendix, width = 8, height = 7)

margins <- 
  d_all |> 
  group_by(case, transition, expectancy) |> 
  summarize(cc = sum(cc), .groups = "drop") |> 
  pivot_wider(names_from = case, values_from = cc) |> 
  arrange(expectancy,transition) |> 
  relocate(expectancy, .before = 1) 


margin_totals <-
  margins |> 
  pivot_longer(`1`:`3`, names_to = "case", values_to = "margin") |> 
  group_by(expectancy, case) |> 
  summarize(Total = sum(margin, na.rm = TRUE), .groups="drop")

margin_long <-
  margins |> 
  pivot_longer(`1`:`3`, names_to = "case", values_to = "margin") 

margin_table <- 
expectancies |> 
  rename(h = hle, u = ule,t = le) |> 
  pivot_longer(h:t, names_to = "expectancy") |> 
  pivot_wider(values_from = value, names_from = sex) |> 
  mutate(diff = f - m) |> 
  select(-f,-m) |> 
  left_join(margin_totals, by = "expectancy") |> 
  mutate(Residual = diff - Total) |> 
  select(-diff)  |> 
  pivot_longer(3:4, names_to = "transition", values_to = "margin") |> 
  bind_rows(margin_long) |> 
  arrange(expectancy, transition) |> 
  mutate(expectancy = case_when(expectancy == "h" ~ "DFLE",
                                expectancy == "u" ~ "DLE",
                                expectancy == "t" ~ "LE",
                                TRUE ~ expectancy)) |> 
  pivot_wider(names_from = case, values_from = margin) 
margin_table |>   
  xtable()

margin_table |> 
  filter(expectancy == "DFLE") |> 
  xtable()



# figwebinar <- 
#   trans |> 
#   rename(edad = age, sexo = sex) |> 
#   pivot_longer(-c(sexo,edad),names_to="transición",values_to = "probabilidad") |> 
#   ggplot(aes(x=edad,y=probabilidad,color=transición, linetype = sexo))+
#   geom_line() +
#   theme_minimal(base_size = 22) +
#   scale_color_discrete(labels=c(HH=expression(p[hh]),
#                                 HU=expression(p[hu]),
#                                 HD=expression(p[hd]),
#                                 UU=expression(p[uu]),
#                                 UH=expression(p[uh]),
#                                 UD=expression(p[ud]) )) 
# figwebinar
expectancies |> 
  rename(h = hle, u = ule,t = le) |> 
  pivot_longer(h:t, names_to = "expectancy") |> 
  pivot_wider(values_from = value, names_from = sex) |> 
  mutate(diff = f - m) |> 
  select(-f,-m) |> 
  left_join(margin_totals, by = "expectancy") |> 
  mutate(Residual = diff - Total)
