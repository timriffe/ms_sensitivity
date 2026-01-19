
# From R1 we have the comment
#"Discuss how other parameterizations (e.g., those incorporating combinations 
# like hh, uh, hd, ud) might affect decomposition results."

# we should first check that additivity is maintained:
library(tidyverse)

dec <- read_csv("all_decompositions.csv")
dec |> 
  filter(expectancy == "h",
         ((case == 3 & transition %in% c("HH","HD","init"))|
         (case == 2 & transition %in% c("UH","UD")) )) |> 
  summarize(cc=sum(cc))
dec |> 
  filter(expectancy == "h") |> 
  summarize(cc=sum(cc), .by=case)

# yes, additivity is maintained, which we could surmise from
# the equal age margins.

dec |> 
  filter(expectancy == "h",
         ((case == 3 & transition %in% c("HH","HD","init"))|
            (case == 2 & transition %in% c("UH","UD")) )) |> 
  ggplot(aes(x=age,y=cc,color=transition)) +
  geom_line()
