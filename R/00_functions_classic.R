
library(data.table)
library(tidyverse)
library(collapse) # for fsubset(), fast subsetting
library(tictoc)
library(tidyfast)

# From Results section:
# We compute the initial conditions at age 50 (ℓ_i(0), where 0 here means age 50)
# assuming that the probabilities we have for age 50 were constant in younger ages.
init_constant <- function(x){
  u <- matrix(x[c("HH","UH","HU","UU")] %>% unlist(),2,byrow=TRUE)
  v <- eigen(u)$vectors[,1]
  init <- v / sum(v)
  names(init) <- c("H","U")
  init
}



pi_block <- function(p, state_from, state_to, age) {
  
  state_fromi  <- state_from[1]
  state_toi    <- state_to[1]
  age          <- c(age, max(age) + 1)
  P            <- diag(p)
  P            <- cbind(rbind(0, P), 0)
  from_names   <- paste(state_fromi, age, sep = "_")
  to_names     <- paste(state_toi,   age, sep = "_")
  dimnames(P)  <- list(to_names, from_names)
  return(P)
}

pi_block_outer <- function(chunk) {
  
  pi_block(chunk[["p"]] |>
             as.double(),
           chunk[["from"]],
           chunk[["to"]],
           chunk[["age"]]) |>
    as_tibble()
}

Ptibble2U <- function(Ptibble) {
  
  age <- Ptibble[["age"]] |>
    unique() |>
    sort()
  
  age <- c(age, age[length(age)] + 1)
  pre <- Ptibble |>
    dt_pivot_longer(-age, 
                    names_to  = "fromto", 
                    values_to = "p") |>
    fmutate(from = substr(fromto, 0, 1),
            to   = substr(fromto, 2, 2)) |>
    fselect(-fromto) |>
    fsubset(to != "D")
  
  states <- pre$from |>
    unique()
  
  pre |>
    group_by(from, to) |>
    nest() |>
    fmutate(data = map(data, ~ .x |>
                         pi_block_outer())) |> 
    dt_pivot_wider(names_from  = from, 
                   values_from = data) |> 
    unnest(cols = all_of(states),
           names_sep = "") |> 
    ungroup() |> 
    fmutate(to = paste(rep(states, each = length(age)), 
                       rep(age,    length(states)), # each
                       sep = "_")) |> 
    column_to_rownames("to") |> 
    as.matrix()
}

# 24-07-2024, updated to add interval arg
Ptibble2N <- function(Ptibble, discount = FALSE, interval = 1) {
  
  U <- Ptibble2U(Ptibble)
  I <- diag(rep(1, nrow(U)))
  N <- solve(I - U) * interval
  if (discount) {
    # untested
    N <- N - (I * interval) / 2
  }
  return(N)
}







