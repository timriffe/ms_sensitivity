
#' function modified from s1()
s1_constrained <- function(hh, hu, uu, uh,
               init = c(H = .97, U = .03), 
               expectancy = c("h", "u", "t", "all"),
               interval = 1){
  # time steps
  n     = length(hh)
  # transient states
  s     = 2

  # survivor stock by state
  lh    = rep(0, n + 1)
  lu    = rep(0, n + 1)
  lh[1] = init["H"]
  lu[1] = init["U"]
  
  # P1 in generic terms
  for (i in 1:n){
    lh[i+1] = lh[i] * hh[i] + lu[i] * uh[i]
    lu[i+1] = lu[i] * uu[i] + lh[i] * hu[i]
    
    # P <- matrix(c(hh[i], uh[i], hu[i], uu[i]), ncol = 2) # same as delta_x_l
    # x_next_mat <- c(lh[i], lu[i]) %*% P      # row-vector update
    # stopifnot(max(abs(c(lh[i+1], lu[i+1]) - as.numeric(x_next_mat))) < 1e-12)
  }
  
  # ------------------------------
 
  delta_x_l = list()
  for (i in 1:n){
    # eq 13 from RevisionAttempt.pdf "direct effect of states on themselves"
    #delta_x_l[[i]] <- matrix(c(hh[i], uh[i], hu[i], uu[i]), ncol = s)
    # RevisionAttempt eq 36
    delta_x_l[[i]] <- matrix(
      c((1 - lh[i]) * hh[i] - lu[i] * uh[i], 
        -lu[i] * uu[i] + (1 - lh[i]) * hu[i], 
       -lh[i] * hh[i] + (1 - lu[i]) * uh[i], 
       (1 - lu[i]) * uu[i] - lh[i] * hu[i]),
                             ncol = s, byrow = TRUE)
  }
  # ------------------------------
  # eq 12 in R0 (implied but not yet included in R1)
  delta_x <- 
    Matrix::bdiag(delta_x_l) |> 
    mpad(side = "l", n = s, value = 0) |> 
    mpad(side = "b", n = s, value = 0)
  # comment out test
  diag(delta_x) <- 1

  # eq 12 R0
  # We use 2 instead of 1 because delta_x has 1s
  # in the diagonal; slightly different from standard 
  # transient matrix.
  d_x = solve(diag(nrow(delta_x)) * 2  - delta_x)
  # d_x <- solve(diag(nrow(delta_x)) - delta_x)
  # we use a shorthand solution for sensitivity to initial conditions.
  # eq 14
  s_in <- s_init(hh = hh, hu = hu, uu = uu, uh = uh,
                 interval = interval)

  # replaces eq 31 R0
  # eq 17 R1
  delta_u_l = list()
  for (i in 1:n){
    # eq 18 R1
    # pi <- matrix(c(1-hh[i], -hh[i],0,0,
    #                -hu[i],1-hu[i],0,0,
    #                0,0,1-uh[i],-uh[i],
    #                0,0,-uu[i],1-uu[i]),4)
    # pi <- matrix(
    #   c(
    #     1 - hh[i], -hu[i], 0, 0,
    #     -hh[i],  1 - hu[i], 0, 0,
    #     0, 0, 1 - uh[i], -uu[i],
    #     0, 0, -uh[i],  1 - uu[i]
    #   ),
    #   nrow = 4, byrow = TRUE
    # )
    # # eq 19 R1
    # li <- matrix(c(lh[i],0,lu[i],0,0,lh[i],0,lu[i]),
    #              ncol = s)
    # eq 16 from RevisionAttempt.pdf (identical to above)
    # derivative wrt transitions
    uli <- matrix(
      c(lh[i] * (1 - hh[i]),   -lh[i] * hu[i] ,          # HH
        lu[i] * (1 - uh[i]),   -lu[i] * uu[i] ,          # UH
       -lu[i] * uh[i],          lu[i] * (1 - uu[i]),     # UU
       -lh[i] * hh[i],          lh[i] * (1 - hu[i])),    # HU
      ncol = s,
      byrow = TRUE
    )
    # eq 17 R1
    # delta_u_l[[i]] = pi %*% li
    delta_u_l[[i]] = uli
  }

  # R1 eq 25
  delta_u =
    Matrix::bdiag(delta_u_l) |> 
    mpad(side = "l", n = s, value = 0) |> 
    mpad(side = "b", n = s^2, value = 0)

  # full sensitivity,
  # eq 27 R1
  sen = delta_u %*% d_x
  
  # ----------------------------------------------- #
  # The rest of this is just tidy management of the 
  # computed sensitivity
  
  # first 4 rows = first age,
  age_from      = rep(0:n,each=s^2)
  
  # TR: this ordering might be fragile
  state_from_to = rep(c("HH","UH","UU","HU"),n+1)
  rownames(sen) = paste(state_from_to,age_from,sep="_")
  # u1,u2,u3,u4
  
  age_to        = rep(0:n,each=s)
  effect_on     = rep(c("H","U"),n+1)
  colnames(sen) = paste(effect_on,age_to,sep="_")
  

  
  
  # return output reformatted
  senl = 
    sen |> 
    as.matrix() |> 
    as.data.frame() |> 
    rownames_to_column("trans_age") |> 
    pivot_longer(-trans_age,
                 names_to="state_age",
                 values_to ="effect") |> 
    separate_wider_delim(trans_age,
                         names=c("transition","agefrom"),
                         delim="_") |> 
    separate_wider_delim(state_age,
                         names=c("state","age"),
                         delim="_") |> 
    mutate(age = as.integer(age),
           agefrom = as.integer(agefrom)) |> 
    filter(age >= agefrom) |> 
    mutate(effect = effect * interval)
  
  # With respect to which expectancy?
  if (expectancy == "h"){
    out =
      senl |> 
      filter(state == "H") |> 
      group_by(transition, agefrom) |> 
      summarize(effect = sum(effect),
                .groups = "drop")|> 
      rename(age = agefrom)
    
    # add on initial conditions
    # comp_effect = sum(init_effects[init_H_ind])
    out = tibble(age = 0,
                 transition = "init",
                 effect = s_in["h"]) |> 
      bind_rows(out)
  }
  if (expectancy == "u"){
    out =
      senl |> 
      filter(state == "U") |> 
      group_by(transition, agefrom) |> 
      summarize(effect = sum(effect),
                .groups = "drop")|> 
      rename(age = agefrom)
    
    # add on initial conditions
    out = tibble(age = 0,
                 transition = "init",
                 effect = s_in["u"]) |> 
      bind_rows(out)
  }
  if (expectancy == "t"){
    out =
      senl |> 
      group_by(transition, agefrom) |> 
      summarize(effect = sum(effect),
                .groups = "drop") |> 
      rename(age = agefrom)
    
    # add on initial conditions
    out = tibble(age = 0,
                 transition = "init",
                 effect = s_in["t"]) |> 
      bind_rows(out)
  }
  if (expectancy == "all"){
    H =
      senl |> 
      filter(state == "H") |> 
      group_by(transition, agefrom) |> 
      summarize(effect = sum(effect),
                .groups = "drop")|> 
      rename(age = agefrom)|> 
      mutate(expectancy = "h", .before = 1)
    U =
      senl |> 
      filter(state == "U") |> 
      group_by(transition, agefrom) |> 
      summarize(effect = sum(effect),
                .groups = "drop")|> 
      rename(age = agefrom) |> 
      mutate(expectancy = "u", .before = 1)
    Tot =
      out =
      senl |> 
      group_by(transition, agefrom) |> 
      summarize(effect = sum(effect),
                .groups = "drop") |> 
      rename(age = agefrom) |> 
      mutate(expectancy = "t", .before = 1)
    out = bind_rows(H,U,Tot)
    # add on initial effects
    # cH = sum(init_effects[init_H_ind])
    # cU = sum(init_effects[init_U_ind])
    # cT = cH + cU
    out = tibble(expectancy = c("h","u","t"),
                 age = c(0,0,0),
                 transition = rep("init",3),
                 effect = s_in[c("h","u","t")]) |> 
      bind_rows(out)
    
  }
  out
}

s2_constrained <- function(hd, hu, ud, uh, 
               init = c(H = .97, U = .03), 
               expectancy = c("h","u","t","all"),
               interval = 1){
  # time steps
  n     = length(hd)
  # transient states
  s     = 2
  
  # survivor stock by state
  lh    = rep(0,n+1)
  lu    = rep(0,n+1)
  lh[1] = init["H"]
  lu[1] = init["U"]
  
  # P2 in generic terms
  for (i in 1:n){
    lh[i+1] = lh[i] * (1 - hd[i] - hu[i]) + lu[i] * uh[i]
    lu[i+1] = lu[i] * (1 - ud[i] - uh[i]) + lh[i] * hu[i]
  }
  
  # eq 9
  delta_x_l = list()
  for (i in 1:n){
    # eq 33 sitting in here
    # delta_x_l[[i]] = matrix(c(1-hd[i]-hu[i],
    #                           uh[i],
    #                           hu[i],
    #                           1-ud[i]-uh[i]),
    #                         ncol=s)
    # RevisionAttempt eq 36
    delta_x_l[[i]] =matrix(
      c((1 - lh[i]) * (1 - hu[i] - hd[i]) - lu[i] * uh[i], 
        -lu[i] * (1 - ud[i] - uh[i]) + (1 - lh[i]) * hu[i], 
        -lh[i] * (1 - hu[i] - hd[i]) + (1 - lu[i]) * uh[i], 
        (1 - lu[i]) * (1 - ud[i] - uh[i]) - lh[i] * hu[i]),
      ncol = s, byrow = TRUE)
  }
  
  # eq 9 (cont)
  delta_x <- 
    Matrix::bdiag(delta_x_l) |> 
    mpad(side = "l", n = s, value = 0) |> 
    mpad(side = "b", n = s, value = 0)
  diag(delta_x) = 1
  
  
  # eq 12
  d_x = solve(diag(nrow(delta_x)) * 2  - delta_x) * interval
  
  # eq 14
  # initial effects sorting code is ad hoc, 
  # it handles only this state space. General
  # code solution needed
  # init_effects <- d_x[1:s, ]
  # init_H_ind   <- col(init_effects) %% 2 == 1
  # init_U_ind   <- col(init_effects) %% 2 == 0
  s_in <- s_init(hh = (1 - hd - hu),
                 hu = hu,
                 uu = (1 - ud - uh),
                 uh = uh,
                 interval = interval)
  
  # eq 16
  delta_u_l = list()
  # 4x2 matrices
  for (i in 1:n){
    # RA eq 26
    delta_u_l[[i]] = matrix(c(-lh[i] * (1 - hu[i] - hd[i]), -lh[i] * hu[i],              #HD
                               lu[i] * (1 - uh[i]),         -lu[i] * (1 - uh[i] - ud[i]),#UH
                              -lu[i] * uh[i],               -lu[i] * (1 - uh[i] - ud[i]),#UD
                              -lh[i] * (1 - hu[i] - hd[i]),  lh[i] * (1 - hu[i])),       #HU
                            ncol = s,
                            byrow = TRUE)
  }
  
  # eq 16 (cont)
  delta_u =
    Matrix::bdiag(delta_u_l) |> 
    mpad(side="l",n=s,value=0) |> 
    mpad(side="b",n=s^2,value=0)
  
  # full sensitivity,
  # eq 18
  sen = delta_u %*% d_x
  
  # first 4 rows = first age,
  age_from      = rep(0:n,each=s^2)
  state_from_to = rep(c("HD","UH","UD","HU"),n+1)
  rownames(sen) = paste(state_from_to, age_from, sep = "_")
  
  age_to        = rep(0:n,each=s) 
  effect_on     = rep(c("H","U"),n+1)
  colnames(sen) = paste(effect_on, age_to, sep="_")
  
  # return output reformatted
  senl = 
    sen |> 
    as.matrix() |> 
    as.data.frame() |> 
    rownames_to_column("trans_age") |> 
    pivot_longer(-trans_age,
                 names_to="state_age",
                 values_to ="effect") |> 
    separate_wider_delim(trans_age,
                         names=c("transition","agefrom"),
                         delim="_") |> 
    separate_wider_delim(state_age,
                         names=c("state","age"),
                         delim="_") |> 
    mutate(age = as.numeric(age),
           agefrom = as.numeric(agefrom)) |> 
    filter(age >= agefrom) 
  
  # TR: here 2 decimals allows for quarters, 
  # but not a super great solution
  # if (interval == 1){
  #   senl <- senl |> 
  #     mutate(age = as.integer(age), 
  #            agefrom = as.integer(agefrom))
  # } else {
  #   senl <- senl |> 
  #     mutate(age =  round(age,2), 
  #            agefrom = round(agefrom,2))
  # }
  
  if (expectancy == "h"){
    out =
      senl |> 
      filter(state == "H") |> 
      group_by(transition, agefrom) |> 
      summarize(effect = sum(effect),
                .groups = "drop")|> 
      rename(age = agefrom)
    
    # add on initial conditions
    # comp_effect = sum(init_effects[init_H_ind])
    out = tibble(age = 0,
                 transition = "init",
                 effect = s_in["h"]) |> 
      bind_rows(out) |> 
      mutate(age = age * interval)
  }
  if (expectancy == "u"){
    out =
      senl |> 
      filter(state == "U") |> 
      group_by(transition, agefrom) |> 
      summarize(effect = sum(effect),
                .groups = "drop")|> 
      rename(age = agefrom)
    
    # add on initial conditions
    # comp_effect = sum(init_effects[init_U_ind])
    out = tibble(age = 0,
                 transition = "init",
                 effect = s_in["u"]) |> 
      bind_rows(out)|> 
      mutate(age = age * interval)
  }
  if (expectancy == "t"){
    out =
      senl |> 
      group_by(transition, agefrom) |> 
      summarize(effect = sum(effect),
                .groups = "drop") |> 
      rename(age = agefrom)
    
    # add on initial conditions
    # comp_effect = sum(init_effects)
    out = tibble(age = 0,
                 transition = "init",
                 effect = s_in["t"]) |> 
      bind_rows(out)|> 
      mutate(age = age * interval)
  }
  if (expectancy == "all"){
    H =
      senl |> 
      filter(state == "H") |> 
      group_by(transition, agefrom) |> 
      summarize(effect = sum(effect),
                .groups = "drop")|> 
      rename(age = agefrom)|> 
      mutate(expectancy = "h", .before = 1)
    U =
      senl |> 
      filter(state == "U") |> 
      group_by(transition, agefrom) |> 
      summarize(effect = sum(effect),
                .groups = "drop")|> 
      rename(age = agefrom) |> 
      mutate(expectancy = "u", .before = 1)
    Tot =
      out =
      senl |> 
      group_by(transition, agefrom) |> 
      summarize(effect = sum(effect),
                .groups = "drop") |> 
      rename(age = agefrom) |> 
      mutate(expectancy = "t", .before = 1)
    out = bind_rows(H,U,Tot)
    # add on initial effects
    # cH = sum(init_effects[init_H_ind])
    # cU = sum(init_effects[init_U_ind])
    # cT = cH + cU
    out = tibble(expectancy = c("h","u","t"),
                 age = c(0,0,0),
                 transition = rep("init",3),
                 effect = s_in[c("h","u","t")]) |> 
      bind_rows(out)|> 
      mutate(age = age * interval)
    
  }
  out
  
}


s3_constrained <- function(hh, uu, ud, hd, 
               init = c(H = .97, U = .03), 
               expectancy = c("h","u","t","all"),
               interval = 1){
  # time steps
  n     = length(hd)
  # transient states
  s     = 2
  
  
  # survivor stock by state
  lh    = rep(0,n+1)
  lu    = rep(0,n+1)
  lh[1] = init["H"]
  lu[1] = init["U"]
  
  for (i in 1:n){
    lh[i+1] = lh[i] * hh[i] + lu[i] * (1 - ud[i] - uu[i])
    lu[i+1] = lu[i] * uu[i] + lh[i] * (1 - hd[i] - hh[i])
  }
  
  # eq 9
  delta_x_l = list()
  for (i in 1:n){
    # R0 eq 39 sitting in here
    # eq 29 of RevisionAttempt
    # delta_x_l[[i]] = matrix(c(hh[i],
    #                           1-ud[i]-uu[i],
    #                           1-hd[i]-hh[i],
    #                           uu[i]),
    #                         ncol=s)
    # Eq 36 RevisionAttempt
    delta_x_l[[i]] =matrix(
      c((1 - lh[i]) * hh[i] - lu[i] * (1 - ud[i] - uu[i]), 
        -lu[i] * uu[i] + (1 - lh[i]) * (1 - hh[i] - hd[i]), 
        -lh[i] * hh[i] + (1 - lu[i]) * (1 - ud[i] - uu[i]), 
        (1 - lu[i]) * uu[i] - lh[i] * (1 - hh[i] - hd[i])),
      ncol = s, byrow = TRUE)
  }
  
  # eq 9 (cont)
  delta_x <- 
    Matrix::bdiag(delta_x_l) |> 
    mpad(side = "l", n = s, value = 0) |> 
    mpad(side = "b", n = s, value = 0)
  diag(delta_x) = 1
  
  
  # eq 12
  d_x = solve(diag(nrow(delta_x)) * 2  - delta_x)
  
  # eq 14 (see equivalent expression commented out in s1())
  s_in <- s_init(hh = hh,
                 hu = (1 - hh - hd),
                 uu = uu,
                 uh = (1 - ud - uu))
  
  # eq 16
  delta_u_l = list()
  # 4x2 matrices
  for (i in 1:n){
    # RevisionAttempt eq 32
    delta_u_l[[i]] = matrix(c(lh[i] * (1 - hh[i]),          -lh[i] * (1 - hd[i] - hh[i]), # HH
                              -lu[i] * (1 - uu[i] - ud[i]),  lu[i] * (1 - uu[i]),      # UU
                              -lu[i] * (1 - uu[i] - ud[i]), -lu[i] * uu[i],            # UD
                              -lh[i] * hh[i],               -lh[i] * (1 - hd[i] - hh[i])),# HD
                            ncol = s,
                            byrow = TRUE)
  }
  # eq 16 (cont)
  delta_u =
    Matrix::bdiag(delta_u_l) |> 
    mpad(side="l",n=s,value=0) |> 
    mpad(side="b",n=s^2,value=0)
  
  # full sensitivity,
  # eq 18
  sen = delta_u %*% d_x
  
  # first 4 rows = first age,
  age_from      = rep(0:n,each=s^2)
  state_from_to = rep(c("HH","UU","UD","HD"),n+1)
  rownames(sen) = paste(state_from_to, age_from, sep = "_")
  # u1,u2,u3,u4
  
  age_to        = rep(0:n,each=s)
  effect_on     = rep(c("H","U"),n+1)
  colnames(sen) = paste(effect_on, age_to, sep = "_")
  
  # return output reformatted
  senl = 
    sen |> 
    as.matrix() |> 
    as.data.frame() |> 
    rownames_to_column("trans_age") |> 
    pivot_longer(-trans_age,
                 names_to="state_age",
                 values_to ="effect") |> 
    separate_wider_delim(trans_age,
                         names=c("transition","agefrom"),
                         delim="_") |> 
    separate_wider_delim(state_age,
                         names=c("state","age"),
                         delim="_") |> 
    mutate(age = as.integer(age),
           agefrom = as.integer(agefrom)) |> 
    filter(age >= agefrom) |> 
    mutate(effect = effect * interval)
  
  if (expectancy == "h"){
    out =
      senl |> 
      filter(state == "H") |> 
      group_by(transition, agefrom) |> 
      summarize(effect = sum(effect),
                .groups = "drop")|> 
      rename(age = agefrom)
    
    # add on initial conditions
    # comp_effect = sum(init_effects[init_H_ind])
    out = tibble(age = 0,
                 transition = "init",
                 effect = s_in["h"]) |> 
      bind_rows(out)
  }
  if (expectancy == "u"){
    out =
      senl |> 
      filter(state == "U") |> 
      group_by(transition, agefrom) |> 
      summarize(effect = sum(effect),
                .groups = "drop")|> 
      rename(age = agefrom)
    
    # add on initial conditions
    # comp_effect = sum(init_effects[init_U_ind])
    out = tibble(age = 0,
                 transition = "init",
                 effect = s_in["u"]) |> 
      bind_rows(out)
  }
  if (expectancy == "t"){
    out =
      senl |> 
      group_by(transition, agefrom) |> 
      summarize(effect = sum(effect),
                .groups = "drop") |> 
      rename(age = agefrom)
    
    # add on initial conditions
    # comp_effect = sum(init_effects)
    out = tibble(age = 0,
                 transition = "init",
                 effect = s_in["t"]) |> 
      bind_rows(out)
  }
  if (expectancy == "all"){
    H =
      senl |> 
      filter(state == "H") |> 
      group_by(transition, agefrom) |> 
      summarize(effect = sum(effect),
                .groups = "drop")|> 
      rename(age = agefrom)|> 
      mutate(expectancy = "h", .before = 1)
    U =
      senl |> 
      filter(state == "U") |> 
      group_by(transition, agefrom) |> 
      summarize(effect = sum(effect),
                .groups = "drop")|> 
      rename(age = agefrom) |> 
      mutate(expectancy = "u", .before = 1)
    Tot =
      out =
      senl |> 
      group_by(transition, agefrom) |> 
      summarize(effect = sum(effect),
                .groups = "drop") |> 
      rename(age = agefrom) |> 
      mutate(expectancy = "t", .before = 1)
    out = bind_rows(H,U,Tot)
    # add on initial effects
    # cH = sum(init_effects[init_H_ind])
    # cU = sum(init_effects[init_U_ind])
    # cT = cH + cU
    out = tibble(expectancy = c("h","u","t"),
                 age = c(0,0,0),
                 transition = rep("init",3),
                 effect = s_in[c("h","u","t")]) |> 
      bind_rows(out)
    
  }
  out
}


s1t_constrained <- function(data,init,expectancy, interval = 1){
  pt <-
    data |> 
    select(age, transition, p) |> 
    pivot_wider(names_from = transition, values_from = p)
  if (missing(init))  init = init_constant(pt[1,])
  s1_constrained(
     hh = pt$HH, 
     hu = pt$HU, 
     uu = pt$UU, 
     uh = pt$UH, 
     init = init, 
     expectancy = expectancy,
     interval = interval)
}
s2t_constrained <- function(data,init,expectancy, interval = 1){
  pt <-
    data |> 
    select(age, transition, p) |> 
    pivot_wider(names_from = transition, values_from = p)
  if (missing(init))  init = init_constant(pt[1,])
  s2_constrained(
    hd = pt$HD, 
    hu = pt$HU, 
    ud = pt$UD, 
    uh = pt$UH, 
    init = init, 
    expectancy = expectancy,
    interval = interval)
}
s3t_constrained <- function(data,init,expectancy, interval = 1){
  pt <-
    data |> 
    select(age, transition, p) |> 
    pivot_wider(names_from = transition, values_from = p)
  if (missing(init))  init = init_constant(pt[1,])
  s3_constrained(
    hd = pt$HD, 
    hh = pt$HH, 
    ud = pt$UD, 
    uu = pt$UU, 
    init = init, 
    expectancy = expectancy,
    interval = interval)
}

data <- trans |> 
  pivot_longer(-c(sex,age),
               names_to = "transition",
               values_to = "p") |> 
  filter(sex=="m")

s1all_constrained <-
  trans |> 
  pivot_longer(-c(sex,age),
               names_to = "transition",
               values_to = "p") |> 
  group_by(sex) %>%
  group_modify(~s1t_constrained(data = .x, expectancy = "all", interval = 1)) |> 
  ungroup() |> 
  mutate(case = 1, 
         .before = 1)

s1all_constrained |> 
  filter(transition != "init",
         expectancy == "h") |> 
  ggplot(aes(x=age,y=effect, color = transition)) +
  geom_line() +
  facet_wrap(~sex)

s2all_constrained <-
  trans |> 
  pivot_longer(-c(sex,age),
               names_to = "transition",
               values_to = "p") |> 
  group_by(sex) %>%
  group_modify(~s2t_constrained(data = .x, expectancy = "all", interval = 1)) |> 
  ungroup() |> 
  mutate(case = 2, 
         .before = 1)
s3all_constrained <-
  trans |> 
  pivot_longer(-c(sex,age),
               names_to = "transition",
               values_to = "p") |> 
  group_by(sex) %>%
  group_modify(~s3t_constrained(data = .x, expectancy = "all", interval = 1)) |> 
  ungroup() |> 
  mutate(case = 3, 
         .before = 1)


source("R/00_sensitivity_functions.R")
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
s2all <-
  trans |> 
  pivot_longer(-c(sex,age),
               names_to = "transition",
               values_to = "p") |> 
  group_by(sex) %>%
  group_modify(~s2t(data = .x, expectancy = "all", interval = 1)) |> 
  ungroup() |> 
  mutate(case = 2, 
         .before = 1)
s3all <-
  trans |> 
  pivot_longer(-c(sex,age),
               names_to = "transition",
               values_to = "p") |> 
  group_by(sex) %>%
  group_modify(~s3t(data = .x, expectancy = "all", interval = 1)) |> 
  ungroup() |> 
  mutate(case = 3, 
         .before = 1)

unconstrained_all <-
  bind_rows(s1all,s2all,s3all) |> 
  mutate(version = "unconstrained",.before=1)

  


constrained_all <-
  bind_rows(s1all_constrained,s2all_constrained,s3all_constrained) |> 
  mutate(version = "constrained",.before=1)

sen_all <- bind_rows(unconstrained_all,constrained_all)

sen_all |> 
  filter(transition != "init",
         expectancy =="h",
         sex == "f",
         case == 2) |> 
  ggplot(aes(x=age,y=effect,color = transition)) +
  geom_line(linewidth=1) +
  theme_minimal() +
  labs(y="sensitivity") +
  facet_grid(vars(version),vars(case)) +
  theme(strip.text = element_text(size=14),
        axis.title = element_text(size=14)) +
  ylim(-1,1)

##############################

# Derive decompositions
init_f <- trans |> filter(sex== "f") |> slice(1) |> 
  init_constant()
init_m <- trans |> filter(sex== "m") |> slice(1) |> 
  init_constant()
init_row <- tibble(age = 50, 
                   transition = "init", 
                   m = init_m[["H"]],
                   f = init_f[["H"]])
trans_dec <-
  trans |> 
  pivot_longer(-c(sex,age),names_to = "transition",values_to="p") |> 
  pivot_wider(names_from = sex, values_from = p) |> 
  bind_rows(init_row) |> 
  arrange(age, transition) |> 
  mutate(delta = f - m,
         p = (m + f) / 2)

# what are the empirical gaps we want to explain?
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
                     expectancy = "t")) |> 
  pivot_longer(-sex,names_to = "expectancy", values_to = "e50") |> 
  pivot_wider(names_from = sex, values_from = e50) |> 
  mutate(Delta = f - m)

expectancies

# can we recover these gaps using constrained sensitivities?
trans_dec |> 
  filter(transition != "init") |> 
  s1t_constrained(expectancy = "all") |> 
  mutate(age = age + 50) |> 
  left_join(trans_dec) |> 
  mutate(cc = effect * delta) |> 
  group_by(expectancy) |> 
  summarize(Delta = sum(cc, na.rm = TRUE)) # no

trans_dec |> 
  filter(transition != "init") |> 
  s2t_constrained(expectancy = "all") |> 
  mutate(age = age + 50) |> 
  left_join(trans_dec) |> 
  mutate(cc = effect * delta) |> 
  group_by(expectancy) |> 
  summarize(Delta = sum(cc, na.rm = TRUE)) # no

trans_dec |> 
  filter(transition != "init") |> 
  s3t_constrained(expectancy = "all") |> 
  mutate(age = age + 50) |> 
  left_join(trans_dec) |> 
  mutate(cc = effect * delta) |> 
  group_by(expectancy) |> 
  summarize(Delta = sum(cc, na.rm = TRUE)) # no
