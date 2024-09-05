
# these are the three basic parameterizations of HLE, ULE, LE 

# P1
f1 <- function(hh, hu, uu, uh, 
               init = c(H=.9,U=.1), 
               expectancy = c("h","u","t"), 
               interval = 1){
  expectany = match.arg(expectancy)
  n     = length(hh)
  lh    = rep(0,n+1)
  lu    = rep(0,n+1)
  lh[1] = init["H"]
  lu[1] = init["U"]
  
  for (i in 1:n){
    # eq 1a
    lh[i+1] <- lh[i] * hh[i] + lu[i] * uh[i]
    # eq 1b
    lu[i+1] <- lu[i] * uu[i] + lh[i] * hu[i]
  }
  
  if (expectancy == "h"){
    # 3a
    out <- sum(lh) 
  }
  if (expectancy == "u"){
    # 3b
    out <- sum(lu) 
  }
  if (expectancy == "t"){
    out <- sum(lh) + sum(lu)
  }
  out * interval
}

# P2
f2 <- function(hd, hu, ud, uh, 
               init = c(H=.9,U=.1), 
               expectancy = c("h","u","t"), 
               interval = 1){
  expectany = match.arg(expectancy)
  n     = length(hd)
  lh    = rep(0,n+1)
  lu    = rep(0,n+1)
  lh[1] = init["H"]
  lu[1] = init["U"]
  
  # for (i in 1:n){       
  #   # 4a
  #   lh[i+1] <- lh[i] * (1-hd[i]-hu[i]) + lu[i] * uh[i]
  #   # 4b
  #   lu[i+1] <- lu[i] * (1-ud[i]-uh[i]) + lh[i] * hu[i]
  # }
  for (i in 1:n){ 
    # 4a
    lh[i+1] <- lh[i] - lh[i] * (hd[i]+hu[i]) + lu[i] * uh[i]
    # 4b
    lu[i+1] <- lu[i] - lu[i] * (ud[i]+uh[i]) + lh[i] * hu[i]
  }
  if (expectancy == "h"){
    # 3a
    out <- sum(lh)
  }
  if (expectancy == "u"){
    # 3b
    out <- sum(lu)
  }
  if (expectancy == "t"){
    out <- sum(lh) + sum(lu)
  }
  out * interval
}

# P3
f3 <- function(hd, hh, ud, uu, 
               init = c(H=.9,U=.1), 
               expectancy = c("h","u","t"),
               interval = 1){
  expectany = match.arg(expectancy)
  n     = length(hd)
  lh    = rep(0,n+1)
  lu    = rep(0,n+1)
  lh[1] = init["H"]
  lu[1] = init["U"]
  
  for (i in 1:n){
    # 5a
    lh[i+1] <- lh[i] * hh[i] + lu[i] * (1 - ud[i] - uu[i])
    # 5b
    lu[i+1] <- lu[i] * uu[i] + lh[i] * (1 - hd[i] - hh[i])
  }
  
  if (expectancy == "h"){
    # 3a
    out <- sum(lh)
  }
  if (expectancy == "u"){
    # 3b
    out <- sum(lu)
  }
  if (expectancy == "t"){
    out <- sum(lh) + sum(lu)
  }
  out * interval
}

# ------------------------- #
# Wrap for numeric gradients or generalized decomposition.
# just need to be careful putting pars in the same order as the 
# arguments specified in the corresponding function.
# w stands for wrapper. Here, first element of pars should be initH.
f1w <- function(pars, 
                expectancy = "h", 
                interval = 1){
  init        <- c(pars[1],1-pars[1])
  names(init) <- c("H","U")
  pars        <- pars[-1]
  n           <- length(pars)
  dim(pars)   <- c(n/4,4)

  colnames(pars) <- c("HH","HU","UU","UH")
  #if (missing(init)) init <- init_constant(pars[1,])
  f1(hh = pars[,"HH"], 
     hu =pars[,"HU"],
     uu = pars[,"UU"],
     uh = pars[,"UH"],
     init = init,
     expectancy = expectancy,
     interval = interval)
}

f2w <- function(pars,
                expectancy = "h",
                interval = 1){
  init        <- c(pars[1],1-pars[1])
  names(init) <- c("H","U")
  pars        <- pars[-1]
  n           <- length(pars)
  dim(pars)   <- c(n/4,4)
  
  colnames(pars) <- c("HD", "HU", "UD", "UH")
  # if (missing(init)) {
  #   init_pars <- c(HH = 1- pars[1,"HD"] - pars[1,"HU"],
  #                  HU = pars[1,"HU"],
  #                  UU = 1 - pars[1,"UD"] - pars[1,"UH"],
  #                  UH = pars[1,"UH"])
  #   init <- init_constant(init_pars)
  # }
  f2(hd = pars[,"HD"],
     hu = pars[,"HU"],
     ud = pars[,"UD"],
     uh = pars[,"UH"],
     init = init,
     expectancy = expectancy,
     interval = interval)
}

f3w <- function(pars,
                expectancy = "h",
                interval = 1){
  init        <- c(pars[1],1-pars[1])
  names(init) <- c("H","U")
  pars        <- pars[-1]
  n           <- length(pars)
  dim(pars)   <- c(n/4,4)

  colnames(pars) <- c("HD", "HH", "UD", "UU")
  # if (missing(init)) {
  #   init_pars <- c(HH = pars[1,"HH"],
  #                  HU = 1 -  pars[1,"HD"] - pars[1,"HU"],
  #                  UU = pars[1,"UU"],
  #                  UH = 1 - pars[1,"UD"] - pars[1,"UU"])
  #   init <- init_constant(init_pars)
  # }
  f3(hd = pars[,"HD"],
     hh = pars[,"HH"],
     ud = pars[,"UD"],
     uu = pars[,"UU"],
     init = init,
     expectancy = expectancy,
     interval = interval)
}

# and again wrapping parameterizations for tidy calculations,
# t standing for tidy
f1t <- function(data,
                init = c(H = .99, U = .01),
                expectancy = "h",
                interval = 1){
  pt =
    data |> 
    select(age, transition, p) |> 
    pivot_wider(names_from = transition, values_from = p)
  if (missing(init))  init = init_constant(pt[1,])
  f1(hh = pt$HH, 
     hu = pt$HU, 
     uu = pt$UU, 
     uh = pt$UH, 
     init = init, 
     expectancy = expectancy,
     interval = interval)
}

f2t <- function(data,
                init = c(H = .99, U = .01),
                expectancy = "h",
                interval = 1){
  pt <-
    data |> 
    select(age, transition, p) |> 
    pivot_wider(names_from = transition, values_from = p)
  if (missing(init))  init = init_constant(pt[1, ])
  f2(hd = pt$HD, 
     hu = pt$HU, 
     ud = pt$UD, 
     uh = pt$UH, 
     init = init, 
     expectancy = expectancy,
     interval = interval)
}

f3t <- function(data,
                init = c(H = .99, U = .01),
                expectancy = "h",
                interval = 1){
  pt <-
    data |> 
    select(age, transition, p) |> 
    pivot_wider(names_from = transition, values_from = p)
  if (missing(init))  init = init_constant(pt[1,])
  f3(hh = pt$HH, 
     hd = pt$HD, 
     uu = pt$UU, 
     ud = pt$UD, 
     init = init, 
     expectancy = expectancy,
     interval = interval)
}

# helper function for padding block matrices
mpad <- function(X,
                 side = c("l","t","r","b"),
                 n=1,
                 value = 0){
  side = match.arg(side)
  
  if (side %in% c("l","r")){
    block = matrix(value,ncol=n,nrow=nrow(X))
    if (side == "l"){
      out <- cbind(block,X)
    }
    if (side == "r"){
      out <- cbind(X,block)
    }
  }
  if (side %in% c("t","b")){
    block = matrix(value,nrow = n, ncol = ncol(X))
    if (side == "t"){
      out <- rbind(block,X)
    }
    if (side == "b"){
      out <- rbind(X,block)
    }
  }
  out
}

# sensitivity functions using vector args for transitions
# TR: interval arg added, but need to double-check where it needs to be used
s1 <- function(hh, hu, uu, uh, 
               init = c(H = .97, U = .03), 
               expectancy = c("h", "u", "t", "all"),
               interval = 1){
  # time steps
  n     = length(hh)
  # transient states
  s     = 2
  
  u1    = hh[1:n]
  u2    = uh[1:n]
  u3    = uu[1:n]
  u4    = hu[1:n]
  
  # survivor stock by state
  x1    = rep(0, n + 1)
  x2    = rep(0, n + 1)
  x1[1] = init["H"]
  x2[1] = init["U"]
  
  # P1 in generic terms
  for (i in 1:n){
    x1[i+1] = x1[i] * u1[i] + x2[i] * u2[i]
    x2[i+1] = x2[i] * u3[i] + x1[i] * u4[i]
  }
  
  # eq 9
  delta_x_l = list()
  for (i in 1:n){
    # eq 28 sitting in here
    delta_x_l[[i]] = matrix(c(u1[i],
                              u2[i],
                              u4[i],
                              u3[i]),
                            ncol=s)
  }
  
  # eq 9 (cont)
  delta_x <- 
    Matrix::bdiag(delta_x_l) |> 
    mpad(side = "l", n = s, value = 0) |> 
    mpad(side = "b", n = s, value = 0)
  diag(delta_x) = 1
  
  
  # eq 12
  # We use 2 instead of 1 because delta_x has 1s
  # in the diagonal; slightly different from standard 
  # transient matrix.
  d_x = solve(diag(nrow(delta_x)) * 2  - delta_x)
  
  # we use a shorthand solution for sensitivity to initial conditions.
  # eq 14
  s_in <- s_init(hh,hu,uu,uh,
                 interval = interval)
  
  # This gives identical values for h and u, which would sum to t,
  # but I think s_init() is somehow slicker
  # evens <- 1:ncol(init_effects) %% 2 == 0
  # odds  <- 1:ncol(init_effects) %% 2 == 1
  # s_in2 <- c(h = init_effects[1, odds] |> sum() - 
  # init_effects[2, odds] |> sum(),
  # u = init_effects[1, evens] |> sum() - 
  #     init_effects[2, evens] |> sum())
  # To do the same thing robustly, we would need costly 
  # name-based selection. So we stick with the s_init() solution.
  
  delta_u_l = list()
  # 4x2 matrices
  # eq 16
  for (i in 1:n){
    # eq 31 sitting in here
    delta_u_l[[i]] = matrix(c(x1[i], 
                              x2[i],
                              0,0,0,0,
                              x2[i],x1[i]),
                            ncol = s)
  }
  # eq 16  (cont)
  delta_u =
    Matrix::bdiag(delta_u_l) |> 
    mpad(side = "l", n = s, value = 0) |> 
    mpad(side = "b", n = s^2, value = 0)
  
  # full sensitivity,
  # eq 18
  sen = delta_u %*% d_x
  
  # ----------------------------------------------- #
  # The rest of this is just tidy management of the 
  # computed sensitivity
  
  # first 4 rows = first age,
  age_from      = rep(0:n,each=s^2)
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


s2 <- function(hd, hu, ud, uh, 
               init = c(H = .97, U = .03), 
               expectancy = c("h","u","t","all"),
               interval = 1){
  # time steps
  n     = length(hd)
  # transient states
  s     = 2
  
  u1    = hd[1:n]
  u2    = uh[1:n]
  u3    = ud[1:n]
  u4    = hu[1:n]
  
  # survivor stock by state
  x1    = rep(0,n+1)
  x2    = rep(0,n+1)
  x1[1] = init["H"]
  x2[1] = init["U"]
  
  # P2 in generic terms
  for (i in 1:n){
    x1[i+1] = x1[i] * (1 - u1[i] - u4[i]) + x2[i] * u2[i]
    x2[i+1] = x2[i] * (1 - u3[i] - u2[i]) + x1[i] * u4[i]
  }
  
  # eq 9
  delta_x_l = list()
  for (i in 1:n){
    # eq 33 sitting in here
    delta_x_l[[i]] = matrix(c(1-u1[i]-u4[i],
                              u2[i],
                              u4[i],
                              1-u3[i]-u2[i]),
                            ncol=s)
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
                 hu,
                 uu = (1 - ud - uh),
                 uh,
                 interval = interval)
  
  # eq 16
  delta_u_l = list()
  # 4x2 matrices
  for (i in 1:n){
    # eq 36 sitting in here
    delta_u_l[[i]] = matrix(c(-x1[i], 
                              x2[i],
                              0,
                              -x1[i],
                              0,
                              -x2[i],
                              -x2[i],
                              x1[i]),
                            ncol = s)
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


s3 <- function(hh, uu, ud, hd, 
               init = c(H = .97, U = .03), 
               expectancy = c("h","u","t","all"),
               interval = 1){
  # time steps
  n     = length(hd)
  # transient states
  s     = 2
  
  u1    = hh[1:n]
  u2    = uu[1:n]
  u3    = ud[1:n]
  u4    = hd[1:n]
  
  # survivor stock by state
  x1    = rep(0,n+1)
  x2    = rep(0,n+1)
  x1[1] = init["H"]
  x2[1] = init["U"]
  
  for (i in 1:n){
    x1[i+1] = x1[i] * u1[i] + x2[i] * (1 - u3[i] - u2[i])
    x2[i+1] = x2[i] * u2[i] + x1[i] * (1 - u4[i] - u1[i])
  }
  
  # eq 9
  delta_x_l = list()
  for (i in 1:n){
    # eq 39 sitting in here
    delta_x_l[[i]] = matrix(c(u1[i],
                              1-u3[i]-u2[i],
                              1-u4[i]-u1[i],
                              u2[i]),
                            ncol=s)
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
    # eq 42
    delta_u_l[[i]] = matrix(c(x1[i], 
                              -x2[i],
                              -x2[i],
                              0,
                              -x1[i],
                              x2[i],
                              0,
                              -x1[i]),
                            ncol = s)
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

# again, but for use in tidy pipeline
s1t <- function(data,init,expectancy, interval = 1){
  pt <-
    data |> 
    select(age, transition, p) |> 
    pivot_wider(names_from = transition, values_from = p)
  if (missing(init))  init = init_constant(pt[1,])
  s1(hh = pt$HH, 
     hu = pt$HU, 
     uu = pt$UU, 
     uh = pt$UH, 
     init = init, 
     expectancy = expectancy,
     interval = interval)
}

s2t <- function(data,init,expectancy, interval = 1){
  pt <-
    data |> 
    select(age, transition, p) |> 
    pivot_wider(names_from = transition, values_from = p)
  if (missing(init))  init = init_constant(pt[1,])
  s2(hd = pt$HD, 
     hu = pt$HU, 
     ud = pt$UD, 
     uh = pt$UH, 
     init = init, 
     expectancy = expectancy,
     interval = interval)
}

s3t <- function(data,init,expectancy,interval=1){
  pt <-
    data |> 
    select(age, transition, p) |> 
    pivot_wider(names_from = transition, values_from = p)
  if (missing(init))  init = init_constant(pt[1,])
  s3(hh = pt$HH, 
     hd = pt$HD, 
     uu = pt$UU, 
     ud = pt$UD, 
     init = init, 
     expectancy = expectancy,
     interval=interval)
}

# numerical sensitivities (gradients, etc), for comparison
# These offer a check that code is doing the right thing. You
# need to call s_init() beforehand or else specify the fraction healthy
# at start
s1n <- function(hh, hu, uu, uh, 
                initH = .99, 
                expectancy = "h",
                interval=1){

  pars1    <- c(initH,hh,hu,uu,uh)
  pd1      <- numDeriv::grad(f1w, 
                             pars1, 
                             expectancy = expectancy,
                             interval = interval)
  
  s_init.   <- pd1[1]
  pd1      <- pd1[-1]
  n        <- length(pd1)/4
  dim(pd1) <- c(n,4)
  colnames(pd1) <- c("HH","HU","UU","UH")
  
  out <-
    pd1 |> 
    as.data.frame() |> 
    mutate(age = 1:n - 1,
           case = "1") |> 
    pivot_longer(1:4, names_to = "transition", values_to = "s")
  
  bind_rows(tibble(age=0,case="1",transition="init",s=s_init.),
            out)
}

s2n <- function(hd, hu, ud, uh, 
                initH = .99, 
                expectancy = "h",
                interval = 1){
  pars2    <- c(initH, hd, hu, ud, uh)
  pd2      <- numDeriv::grad(f2w, 
                             pars2, 
                             expectancy = expectancy,
                             interval = interval)
 
  s_init.   <- pd2[1]
  pd2      <- pd2[-1]
  
  n        <- length(pd2)/4
  dim(pd2) <- c(n,4)
  colnames(pd2) <- c("HD","HU","UD","UH")
  
  out <-
  pd2 |> 
    as.data.frame() |> 
    mutate(age = 1:n - 1,
           case = "2") |> 
    pivot_longer(1:4, names_to = "transition", values_to = "s")
  
  bind_rows(tibble(age=0,case="2",transition="init",s=s_init.),
            out)
}

s3n <- function(hd, hh, ud, uu, 
                initH = .99, 
                expectancy = "h",
                interval = 1){
  pars3    <- c(initH, hd,hh,ud,uu)
  pd3      <- numDeriv::grad(f3w, 
                             pars3, 
                             expectancy = expectancy,
                             interval = interval)
  
  s_init.   <- pd3[1]
  pd3      <- pd3[-1]
 
   n        <- length(pd3)/4
  dim(pd3) <- c(n,4)
  colnames(pd3) <- c("HD","HH","UD","UU")
  
  out <- pd3 |> 
    as.data.frame() |> 
    mutate(age = 1:n - 1,
           case = "3") |> 
    pivot_longer(1:4, names_to = "transition", values_to = "s")
  
  bind_rows(tibble(age=0,case="3",transition="init",s=s_init.),
            out)
}

# wrappers to do numeric derivatives in tidy workflow
s1nt <- function(data,init,expectancy="h", interval = 1){
  pt <-
    data |> 
    select(age, transition, p) |> 
    pivot_wider(names_from = transition, values_from = p)
  if (missing(init))  init = init_constant(pt[1,])
  s1n(hh = pt$HH, 
     hu = pt$HU, 
     uu = pt$UU, 
     uh = pt$UH, 
     initH = init[1], 
     expectancy = expectancy,
     interval = interval)
}

s2nt <- function(data,init,expectancy="h",interval=1){
  pt <-
    data |> 
    select(age, transition, p) |> 
    pivot_wider(names_from = transition, values_from = p)
  if (missing(init))  init = init_constant(pt[1,])
  s2n(hd = pt$HD, 
      hu = pt$HU, 
      ud = pt$UD, 
      uh = pt$UH, 
      initH = init[1], 
      expectancy = expectancy,
      interval = interval)
}

s3nt <- function(data,init,expectancy="h",interval=1){
  pt <-
    data |> 
    select(age, transition, p) |> 
    pivot_wider(names_from = transition, values_from = p)
  if (missing(init))  init = init_constant(pt[1,])
  s3n(hd = pt$HD, 
      hh = pt$HH, 
      ud = pt$UD, 
      uu = pt$UU, 
      initH = init[1], 
      expectancy = expectancy,
      interval=interval)
}
# vectorized wrapper for s1(), etc; these can also be passed to 
# DemoDecomp::ltre(), can help reduce residual to 0. 
# first element of pars should be initH, followed by hh,hu,uu,uh in that order
s1w <- function(func, pars, expectancy="h",interval=1){
  initH <- pars[1]
  pars <- pars[-1]
  init <- c(initH, 1 - initH)
  names(init) <- c("H","U")
  nn = length(pars)
  dim(pars) = c(nn/4,4)
  sen = s1(hh=pars[,1],
           hu=pars[,2],
           uu=pars[,3],
           uh=pars[,4], 
           init = init, 
           expectancy = expectancy,
           interval = interval) |> 
    filter(age < max(age))
  c(filter(sen, transition == "init")$effect,
    filter(sen, transition == "HH")$effect,
    filter(sen, transition == "HU")$effect,
    filter(sen, transition == "UU")$effect,
    filter(sen, transition == "UH")$effect)
}

# first element of pars should be initH, followed by hd,hu,ud,uh in that order
s2w <- function(func, pars, expectancy = "h", interval = 1){
  initH <- pars[1]
  pars <- pars[-1]
  init <- c(H = initH, U = 1 - initH)
  names(init) <- c("H","U")
  nn = length(pars)
  dim(pars) = c(nn/4,4)
  sen = s2(hd=pars[,1],
           hu=pars[,2],
           ud=pars[,3],
           uh=pars[,4], 
           init = init, 
           expectancy = expectancy,
           interval = interval) |> 
    filter(age < max(age))
  c(filter(sen, transition == "init")$effect,
    filter(sen, transition == "HD")$effect,
    filter(sen, transition == "HU")$effect,
    filter(sen, transition == "UD")$effect,
    filter(sen, transition == "UH")$effect)
}

# first element of pars should be initH, followed by hd,hh,ud,uu in that order
s3w <- function(func, pars, expectancy = "h", interval = 1){
  initH <- pars[1]
  pars <- pars[-1]
  init <- c(H = initH, U = 1 - initH)
  names(init) <- c("H","U")
  nn = length(pars)
  dim(pars) = c(nn/4,4)
  sen = s3(hd=pars[,1],
           hh=pars[,2],
           ud=pars[,3],
           uu=pars[,4], 
           init = init, 
           expectancy = expectancy,
           interval = interval) |> 
    filter(age < max(age))
  c(filter(sen, transition == "init")$effect,
    filter(sen, transition == "HD")$effect,
    filter(sen, transition == "HH")$effect,
    filter(sen, transition == "UD")$effect,
    filter(sen, transition == "UU")$effect)
}

# A wrapper that just calculates sensitivity,
# to initial conditions. this uses P1 args,
# just transform args as needed when calling this function.
# compare w eq 14. The intuitive explanation: The gain in
# HLE for a person that switches from u to h in the initial state
# is (HLE | i=h,x=0) - (HLE | i=u,x=0)
# -If mixture is fast over age then initial conditions are less important.
# -If mortality is equal between states then initial conditions have no leverage
# -If health deterioration is irreversible then initial conditions are 
# very important
s_init <- function(hh,hu,uu,uh,interval = 1){
  
  sh <- f1(hh, hu, uu, uh, 
           init = c(H = 1, U = 0), 
           expectancy = "h",
           interval=interval) - 
    f1(hh, hu, uu, uh, 
       init = c(H = 0, U = 1), 
       expectancy = "h",
       interval=interval) 
  
  su <- f1(hh, hu, uu, uh, 
           init = c(H = 1, U = 0), 
           expectancy = "u",
           interval=interval) - 
    f1(hh, hu, uu, uh, 
       init = c(H = 0, U = 1), 
       expectancy = "u",
       interval=interval)
  
  st <- sh + su
  c(h = sh, u = su, t = st)
}

# end


