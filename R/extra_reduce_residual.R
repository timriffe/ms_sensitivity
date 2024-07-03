# This script shows two strategies to reduce the decomposition residual,
# as mentioned in the manuscript. 

# (1) We repeat the exercise in fine interpolation
# intervals, then sum decomposition results over the full interpolation space.
# This is inspired by the Horiuchi et al (2008) method.

# (2) alternatively, one could search for the one unique interpolation
# point between male and female parameters where the sensitivity evaluation 
# implies a decomposition result summing exactly to the gap. This turns out
# to be more efficient than (1), somewhat to our surprise. If a shorthand 
# rule for this point were known, this would be far and away the most 
# efficient and precise approach.

source("R/00_functions_classic.R")
source("R/00_sensitivity_functions.R")
trans <- read_csv("transitions_lievre2003_annual.csv")

# observed sex differences that we want the decomposition results
# to match in the margin.
gaps <- 
  trans |> 
  pivot_longer(-c(sex,age),
               names_to = "transition",
               values_to = "p") |> 
  group_by(sex) %>%
  summarize(h = f1t(data = pick(everything()), 
                      expectancy = "h"),
            u = f1t(data = pick(everything()), 
                      expectancy = "u"),
            t = f1t(data = pick(everything()), 
                     expectancy = "t")) |> 
  pivot_longer(h:t, names_to = "expectancy", values_to = "value") |> 
  pivot_wider(names_from = sex, values_from = value) |> 
  mutate(gap = f-m) |> 
  pull(gap)
names(gaps) <- c("h","u","t")

pt_f <- trans |> filter(sex == "f")
pt_m <- trans |> filter(sex == "m")

init_f <- init_constant(pt_f[1,])
init_m <- init_constant(pt_m[1,])

hhf <- pt_f$HH
huf <- pt_f$HU
hdf <- pt_f$HD
uuf <- pt_f$UU
uhf <- pt_f$UH
udf <- pt_f$UD

hhm <- pt_m$HH
hum <- pt_m$HU
hdm <- pt_m$HD
uum <- pt_m$UU
uhm <- pt_m$UH
udm <- pt_m$UD

# get parameter vectors
pars_p1_f <- c(init_f[1], hhf, huf, uuf, uhf)
pars_p1_m <- c(init_m[1], hhm, hum, uum, uhm)
 
pars_p2_f <- c(init_f[1], hdf, huf, udf, uhf)
pars_p2_m <- c(init_m[1], hdm, hum, udm, uhm)
 
pars_p3_f <- c(init_f[1], hdf, hhf, udf, uuf)
pars_p3_m <- c(init_m[1], hdm, hhm, udm, uum)

# ------------------------------------------ #
# Experiment
# You can repeat this for other expectancies and
# parameterizations. Just be sure to change values
# where you see annotations indicating to do so
Ns <- c(1,5,20,50,100,200)
ccs <- Ns * 0
for (j in 1:length(Ns)){
  # to change parameterization, switch these vectors!
  pars_f <- pars_p1_f
  pars_m <- pars_p1_m
  # -----------------
  
  delta <- (pars_f - pars_m) / Ns[j]
  p_mat <- t(apply(cbind(pars_f,pars_m), 1, function(y, N) {
    xout <- seq(0, 1, length = N + 2)
    xout <- xout[-c(1, length(xout))]
    c(approx(x = c(0, 1), y = y, xout = xout)$y)}, N = Ns[j]))
  if (j == 1){
    p_mat <- t(p_mat)
  }
  smat <- p_mat
  for (i in 1:Ns[j]){
    # change functions s1w and f1w to match parameterization;
    # e.g. s2w;f2w or s3w;f3w
    # change expectancy to "u" or "t" to check those
    smat[,i] <- s1w(pars = p_mat[,i],expectancy = "h", func = f1w)
  }
  ccs[j] <- sum(smat * delta) 
}

plot(Ns, 
     ccs, 
     log = 'x', 
     type = 'o', 
     main = "increasing precision,\n and decreasing returns\nto finer interpolations in LTRE",
     xlab = "N (axis logged)")
# if using a different expectancy, pick out either "u" or "t" from gaps.
abline(a=gaps["h"],b=0, col = "red")
text(3,gaps["h"],"observed gap in HLE",col="red",pos=1)

# End linear-integral-w-analytic-sensitivity exercise

# ------------------------------------------------------------#
# Start: optimized parameter mean where we average parameters 
# according to some weight that is constant over all ages and transitions.
# in this case, 'pf' means proportion female, such that the avg parameter
# turns out to be f * pf + m * (1 - pf). We optimize pf such that the 
# gap is exactly reproduced. pf turns out to be same over all cases
# for a given expectancy.

transfm <- 
  trans |> 
  pivot_longer(-c(age,sex),names_to = "transition", values_to = "p") |> 
  pivot_wider(names_from = sex, values_from = p) 

gap_min <- function(pf, 
                    transfm, 
                    gap, 
                    case = 1, 
                    expectancy = "h", 
                    return_dec = FALSE){
  
  s_fun <- switch(case,
                  s1t,
                  s2t,
                  s3t)
  dat <- 
  transfm |> 
    mutate(p_avg = f*pf + m * (1-pf),
           delta = f-m)
  
  # a bit knarly to handle initial conditions
  # first, for the given parameter weighting, we get the init
  p0 <- dat |> 
    filter(age == min(age)) |> 
    select(transition,p_avg) |> 
    pivot_wider(names_from = transition, values_from = p_avg) |> 
    unlist()
  
  p_i <- list(HH = p0["HH"],
             HU = p0["HU"],
             UH = p0["UH"],
             UU = p0["UU"])
  # second, to calculate the delta for initH
  init <- init_constant(p_i)
  
  initm <- transfm |> 
    select(-f) |> 
    pivot_wider(names_from = transition, values_from = m) |> 
    slice(1) |> 
    init_constant()
  
  initf <- transfm |> 
    select(-m) |> 
    pivot_wider(names_from = transition, values_from = f) |> 
    slice(1) |> 
    init_constant()

  init_delta <- initf["H"] - initm["H"]
  
  init_d_t <- tibble(age = 50, transition = "init", delta = init_delta)
  
  # all deltas, including 
  datb <- dat |> 
    select(age,transition, delta) |> 
    bind_rows(init_d_t)
  
  dec <- dat |> 
    select(age, transition, p = p_avg) |> 
    s_fun(init = init, expectancy = expectancy) |> 
    mutate(age= age + 50) |> 
    filter(age < 111) |> 
    left_join(datb, by = join_by(age, transition)) |> 
    mutate(cc = effect * delta) 
  
  if (return_dec){
    dec <-
      dec |> 
      mutate(case = case, 
             expectancy = expectancy)
    return(dec)
  }
  gapi <-
    dec |> 
    pull(cc) |> 
    sum()
 
 abs(gap - gapi)
}

# what is the optimal parameter weighting; i.e. close to, but not exactly 0.5
pf <- optimise(gap_min,lower=.45,upper=.55, transfm=transfm, gap = gaps["h"], expectancy = "h", case = 3, tol = 1e-8)$minimum
pf # meaning if we give female parameters this weight, and male parameters the complement, then the resulting decomposition is almost exact.

d_exact <- gap_min(pf = pf, 
        transfm = transfm,
        expectancy = "h",
        case = 1,
        return_dec = TRUE)

# There, you happy now?
d_exact$cc |> sum() - gaps["h"]

# Outstanding questions that we do not seek to answer in the present analysis:
# 1) The optimized mean between parameters decomposition result equal to the linear-
# integral-using-analytic-sensitivity decomposition result under large N (i.e. 200)?
# 2) if so, is the optimize result to be preferred when exactitude is required, 
# simply on grounds of a speedier evaluation?
# 3) is there a shorthand calculation to arrive at the pf we get from optimizing?







