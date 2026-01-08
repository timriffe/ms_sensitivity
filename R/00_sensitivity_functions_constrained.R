
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
  }
  
  # ------------------------------
  # REVIEW HERE
  # replaces eq 28 R0
  # eq 14 R1
  delta_x_l = list()
  for (i in 1:n){
    pt_eq3 <- matrix(c(hh[i],hu[i],uh[i],uu[i]),ncol=s)
    eq13   <- matrix(c(1-lh[i],-lh[i],-lu[i],1-lu[i]), ncol=s)
    # eq 14 R1
    delta_x_l[[i]] <- eq13 %*% t(pt_eq3)
    
    # R0 eq 28
    # delta_x_l[[i]] = matrix(c(hh[i],
    #                           hu[i],
    #                           uh[i],
    #                           uu[i]),
    #                         ncol=s,
    #                         byrow = TRUE)
  }
  # ------------------------------
  # eq__ # not sure yet
  delta_x <- 
    Matrix::bdiag(delta_x_l) |> 
    mpad(side = "l", n = s, value = 0) |> 
    mpad(side = "b", n = s, value = 0)
  diag(delta_x) <- 1

  # eq__
  # We use 2 instead of 1 because delta_x has 1s
  # in the diagonal; slightly different from standard 
  # transient matrix.
  d_x = solve(diag(nrow(delta_x)) * 2  - delta_x)
  
  # d_x[1,seq(1,124) %% 2 == 1] |> sum() 
  # d_x[2,seq(1,124) %% 2 == 1] |> sum()
  # d_x[1,seq(1,124) %% 2 == 0] |> sum()
  # d_x[2,seq(1,124) %% 2 == 0] |> sum()
  
  
  
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
  
  # replaces eq 31 R0
  # eq 17 R1
  delta_u_l = list()
  for (i in 1:n){
    # eq 18 R1
    pi <- matrix(c(1-hh[i], -hh[i],0,0,
                   -hu[i],1-hu[i],0,0,
                   0,0,1-uh[i],-uh[i],
                   0,0,-uu[i],1-uu[i]),4)

    # eq 19 R1
    li <- matrix(c(lh[i],0,lu[i],0,0,lh[i],0,lu[i]),
                 ncol = s)
    # eq 17 R1
    delta_u_l[[i]] = pi %*% li
  }

  # eq 13 R1  (cont)
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
  state_from_to = rep(c("HH","HU","UH","UU"),n+1)
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
                           expectancy = c("h", "u", "t", "all"),
                           interval = 1){
  
  expectancy <- match.arg(expectancy)
  n <- length(hd)
  s <- 2
  
  # survivor stock by state
  lh <- rep(0, n + 1)
  lu <- rep(0, n + 1)
  lh[1] <- init["H"]
  lu[1] <- init["U"]
  
  # P2 recurrence (column-vector convention)
  for (i in 1:n){
    lh[i+1] <- lh[i] * (1 - hd[i] - hu[i]) + lu[i] * uh[i]
    lu[i+1] <- lu[i] * (1 - uh[i] - ud[i]) + lh[i] * hu[i]
  }
  
  # ------------------------------
  # constrained direct effects of states on themselves
  delta_x_l <- vector("list", n)
  for (i in 1:n){
    # Pa for P2 (column-vector convention)
    p11 <- 1 - hd[i] - hu[i]
    p12 <- uh[i]
    p21 <- hu[i]
    p22 <- 1 - ud[i] - uh[i]
    pt_eq3 <- matrix(c(p11, p21, p12, p22), ncol = s)
    
    eq13 <- matrix(c(1 - lh[i], -lh[i],
                     -lu[i],  1 - lu[i]), ncol = s)
    
    delta_x_l[[i]] <- eq13 %*% t(pt_eq3)
  }
  
  delta_x <- Matrix::bdiag(delta_x_l) |>
    mpad(side = "l", n = s, value = 0) |>
    mpad(side = "b", n = s, value = 0)
  diag(delta_x) <- 1
  
  d_x <- solve(diag(nrow(delta_x)) * 2 - delta_x)
  
  # sensitivity to initial conditions (reuse existing helper)
  # (derive full P1-style hh/uu from constraints)
  hh <- 1 - hd - hu
  uu <- 1 - ud - uh
  s_in <- s_init(hh = hh, hu = hu, uu = uu, uh = uh, interval = interval)
  
  # ------------------------------
  # constrained direct effects of transitions on states (RevisionAttempt Eq. 26)
  delta_u_l <- vector("list", n)
  for (i in 1:n){
    phh <- 1 - hd[i] - hu[i]
    phu <- hu[i]
    puh <- uh[i]
    puu <- 1 - ud[i] - uh[i]
    
    delta_u_l[[i]] <- matrix(
      c(
        -lh[i] * phh,          -lh[i] * phu,          # d/d(phd)
         lu[i] * (1 - puh),    -lu[i] * puu,          # d/d(puh)
        -lu[i] * puh,          -lu[i] * puu,          # d/d(pud)
        -lh[i] * phh,           lh[i] * (1 - phu)     # d/d(phu)
      ),
      ncol = s, byrow = TRUE
    )
  }
  
  delta_u <- Matrix::bdiag(delta_u_l) |>
    mpad(side = "l", n = s, value = 0) |>
    mpad(side = "b", n = s^2, value = 0)
  
  sen <- delta_u %*% d_x
  
  # -----------------------------------------------
  # tidy output (same pattern as s1_constrained)
  age_from <- rep(0:n, each = s^2)
  state_from_to <- rep(c("HD","UH","UD","HU"), n + 1)
  rownames(sen) <- paste(state_from_to, age_from, sep = "_")
  
  age_to <- rep(0:n, each = s)
  effect_on <- rep(c("H","U"), n + 1)
  colnames(sen) <- paste(effect_on, age_to, sep = "_")
  
  senl <- sen |>
    as.matrix() |>
    as.data.frame() |>
    rownames_to_column("trans_age") |>
    pivot_longer(-trans_age, names_to = "state_age", values_to = "effect") |>
    separate_wider_delim(trans_age, names = c("transition","agefrom"), delim = "_") |>
    separate_wider_delim(state_age,  names = c("state","age"),     delim = "_") |>
    mutate(age = as.integer(age), agefrom = as.integer(agefrom)) |>
    filter(age >= agefrom) |>
    mutate(effect = effect * interval)
  
  if (expectancy == "h"){
    out <- senl |>
      filter(state == "H") |>
      group_by(transition, agefrom) |>
      summarize(effect = sum(effect), .groups = "drop") |>
      rename(age = agefrom)
    out <- tibble(age = 0, transition = "init", effect = s_in["h"]) |> bind_rows(out)
  }
  if (expectancy == "u"){
    out <- senl |>
      filter(state == "U") |>
      group_by(transition, agefrom) |>
      summarize(effect = sum(effect), .groups = "drop") |>
      rename(age = agefrom)
    out <- tibble(age = 0, transition = "init", effect = s_in["u"]) |> bind_rows(out)
  }
  if (expectancy == "t"){
    out <- senl |>
      group_by(transition, agefrom) |>
      summarize(effect = sum(effect), .groups = "drop") |>
      rename(age = agefrom)
    out <- tibble(age = 0, transition = "init", effect = s_in["t"]) |> bind_rows(out)
  }
  if (expectancy == "all"){
    H <- senl |>
      filter(state == "H") |>
      group_by(transition, agefrom) |>
      summarize(effect = sum(effect), .groups = "drop") |>
      rename(age = agefrom) |>
      mutate(expectancy = "h", .before = 1)
    
    U <- senl |>
      filter(state == "U") |>
      group_by(transition, agefrom) |>
      summarize(effect = sum(effect), .groups = "drop") |>
      rename(age = agefrom) |>
      mutate(expectancy = "u", .before = 1)
    
    Tot <- senl |>
      group_by(transition, agefrom) |>
      summarize(effect = sum(effect), .groups = "drop") |>
      rename(age = agefrom) |>
      mutate(expectancy = "t", .before = 1)
    
    out <- bind_rows(H, U, Tot)
    out <- tibble(
      expectancy = c("h","u","t"),
      age = c(0,0,0),
      transition = rep("init", 3),
      effect = s_in[c("h","u","t")]
    ) |> bind_rows(out)
  }
  
  out
}


s3_constrained <- function(hh, uu, ud, hd,
                           init = c(H = .97, U = .03),
                           expectancy = c("h", "u", "t", "all"),
                           interval = 1){
  
  expectancy <- match.arg(expectancy)
  n <- length(hh)
  s <- 2
  
  # survivor stock by state
  lh <- rep(0, n + 1)
  lu <- rep(0, n + 1)
  lh[1] <- init["H"]
  lu[1] <- init["U"]
  
  # P3 recurrence (column-vector convention)
  for (i in 1:n){
    lh[i+1] <- lh[i] * hh[i] + lu[i] * (1 - ud[i] - uu[i])
    lu[i+1] <- lu[i] * uu[i] + lh[i] * (1 - hd[i] - hh[i])
  }
  
  # ------------------------------
  # constrained direct effects of states on themselves
  delta_x_l <- vector("list", n)
  for (i in 1:n){
    # Pa for P3 (column-vector convention)
    p11 <- hh[i]
    p12 <- 1 - ud[i] - uu[i]
    p21 <- 1 - hd[i] - hh[i]
    p22 <- uu[i]
    pt_eq3 <- matrix(c(p11, p21, p12, p22), ncol = s)
    
    eq13 <- matrix(c(1 - lh[i], -lh[i],
                     -lu[i],  1 - lu[i]), ncol = s)
    
    delta_x_l[[i]] <- eq13 %*% t(pt_eq3)
  }
  
  delta_x <- Matrix::bdiag(delta_x_l) |>
    mpad(side = "l", n = s, value = 0) |>
    mpad(side = "b", n = s, value = 0)
  diag(delta_x) <- 1
  
  d_x <- solve(diag(nrow(delta_x)) * 2 - delta_x)
  
  # sensitivity to initial conditions (reuse existing helper)
  # (derive hu/uh from constraints)
  hu <- 1 - hh - hd
  uh <- 1 - uu - ud
  s_in <- s_init(hh = hh, hu = hu, uu = uu, uh = uh, interval = interval)
  
  # ------------------------------
  # constrained direct effects of transitions on states (RevisionAttempt Eq. 32)
  delta_u_l <- vector("list", n)
  for (i in 1:n){
    phh <- hh[i]
    puu <- uu[i]
    pud <- ud[i]
    phd <- hd[i]
    phu <- 1 - phh - phd
    puh <- 1 - puu - pud
    
    delta_u_l[[i]] <- matrix(
      c(
        lh[i] * (1 - phh),   -lh[i] * phu,          # d/d(phh)
        -lu[i] * puh,         lu[i] * (1 - puu),    # d/d(puu)
        -lu[i] * puh,        -lu[i] * puu,          # d/d(pud)
        -lh[i] * phh,        -lh[i] * phu           # d/d(phd)
      ),
      ncol = s, byrow = TRUE
    )
  }
  
  delta_u <- Matrix::bdiag(delta_u_l) |>
    mpad(side = "l", n = s, value = 0) |>
    mpad(side = "b", n = s^2, value = 0)
  
  sen <- delta_u %*% d_x
  
  # -----------------------------------------------
  # tidy output
  age_from <- rep(0:n, each = s^2)
  state_from_to <- rep(c("HH","UU","UD","HD"), n + 1)
  rownames(sen) <- paste(state_from_to, age_from, sep = "_")
  
  age_to <- rep(0:n, each = s)
  effect_on <- rep(c("H","U"), n + 1)
  colnames(sen) <- paste(effect_on, age_to, sep = "_")
  
  senl <- sen |>
    as.matrix() |>
    as.data.frame() |>
    rownames_to_column("trans_age") |>
    pivot_longer(-trans_age, names_to = "state_age", values_to = "effect") |>
    separate_wider_delim(trans_age, names = c("transition","agefrom"), delim = "_") |>
    separate_wider_delim(state_age,  names = c("state","age"),     delim = "_") |>
    mutate(age = as.integer(age), agefrom = as.integer(agefrom)) |>
    filter(age >= agefrom) |>
    mutate(effect = effect * interval)
  
  if (expectancy == "h"){
    out <- senl |>
      filter(state == "H") |>
      group_by(transition, agefrom) |>
      summarize(effect = sum(effect), .groups = "drop") |>
      rename(age = agefrom)
    out <- tibble(age = 0, transition = "init", effect = s_in["h"]) |> bind_rows(out)
  }
  if (expectancy == "u"){
    out <- senl |>
      filter(state == "U") |>
      group_by(transition, agefrom) |>
      summarize(effect = sum(effect), .groups = "drop") |>
      rename(age = agefrom)
    out <- tibble(age = 0, transition = "init", effect = s_in["u"]) |> bind_rows(out)
  }
  if (expectancy == "t"){
    out <- senl |>
      group_by(transition, agefrom) |>
      summarize(effect = sum(effect), .groups = "drop") |>
      rename(age = agefrom)
    out <- tibble(age = 0, transition = "init", effect = s_in["t"]) |> bind_rows(out)
  }
  if (expectancy == "all"){
    H <- senl |>
      filter(state == "H") |>
      group_by(transition, agefrom) |>
      summarize(effect = sum(effect), .groups = "drop") |>
      rename(age = agefrom) |>
      mutate(expectancy = "h", .before = 1)
    
    U <- senl |>
      filter(state == "U") |>
      group_by(transition, agefrom) |>
      summarize(effect = sum(effect), .groups = "drop") |>
      rename(age = agefrom) |>
      mutate(expectancy = "u", .before = 1)
    
    Tot <- senl |>
      group_by(transition, agefrom) |>
      summarize(effect = sum(effect), .groups = "drop") |>
      rename(age = agefrom) |>
      mutate(expectancy = "t", .before = 1)
    
    out <- bind_rows(H, U, Tot)
    out <- tibble(
      expectancy = c("h","u","t"),
      age = c(0,0,0),
      transition = rep("init", 3),
      effect = s_in[c("h","u","t")]
    ) |> bind_rows(out)
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
         sex == "f") |> 
  ggplot(aes(x=age,y=effect,color = transition)) +
  geom_line(linewidth=1) +
  theme_minimal() +
  labs(y="sensitivity") +
  facet_grid(vars(version),vars(case))

s3all_constrained |> 
  filter(transition != "init",
         expectancy =="h",
         sex == "f") |> 
  ggplot(aes(x=age,y=effect,color = transition)) +
  geom_line(linewidth=1) +
  theme_minimal() +
  labs(y="sensitivity")
