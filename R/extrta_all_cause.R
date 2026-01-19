
library(DemoDecomp) # on CRAN
library(LEdecomp)   # on CRAN
library(tidyverse)
library(janitor)

grab_data <- FALSE
if (grab_data){
library(HMDHFDplus)
dat <- readHMDweb(CNTRY = "ESP",item = "mltper_1x1",
            username = Sys.getenv("us"), 
            password = Sys.getenv("pw"))  |> 
  filter(Year %in% c(2000,2019), Age <= 100) |> 
  clean_names() |> 
  select(year, age, qx)  |> 
  pivot_wider(names_from = year, values_from = qx)
  dput(dat)
}
# pasted from the above
dat <- structure(list(age = 0:100, 
                      `2000` = c(0.0047, 0.00045, 0.00029, 
0.00022, 0.00017, 0.00019, 0.00019, 0.00016, 0.00015, 0.00015, 
0.00015, 0.00013, 0.00015, 0.00019, 0.00032, 0.00035, 0.00052, 
7e-04, 7e-04, 0.00073, 0.00085, 0.00085, 8e-04, 0.00083, 0.00091, 
9e-04, 0.00098, 9e-04, 0.00106, 0.00101, 0.00105, 0.00122, 0.00127, 
0.00135, 0.00138, 0.00164, 0.00167, 0.00174, 0.0019, 0.0019, 
0.00221, 0.00217, 0.00244, 0.00266, 0.00291, 0.00309, 0.00314, 
0.00377, 0.00382, 0.00418, 0.00489, 0.00491, 0.00543, 0.0059, 
0.00635, 0.00676, 0.00771, 0.00811, 0.00875, 0.00964, 0.01039, 
0.01183, 0.01266, 0.01399, 0.01542, 0.01701, 0.01781, 0.02029, 
0.02139, 0.02406, 0.02643, 0.02928, 0.03227, 0.03562, 0.03884, 
0.04349, 0.04812, 0.05344, 0.0596, 0.06607, 0.0734, 0.08027, 
0.08861, 0.09588, 0.10642, 0.11546, 0.13233, 0.14266, 0.15691, 
0.16783, 0.18371, 0.20978, 0.22078, 0.24228, 0.25855, 0.27478, 
0.29392, 0.31333, 0.33287, 0.35242, 0.37183), 
`2019` = c(0.00284, 0.00023, 0.00015, 9e-05, 0.00011, 7e-05, 6e-05, 7e-05, 6e-05, 
6e-05, 8e-05, 0.00011, 5e-05, 0.00012, 0.00011, 0.00015, 0.00017, 
0.00019, 0.00027, 0.00026, 0.00028, 0.00031, 3e-04, 0.00037, 
0.00035, 0.00042, 0.00036, 0.00043, 0.00039, 0.00039, 0.00047, 
5e-04, 0.00051, 0.00047, 0.00052, 0.00054, 6e-04, 0.00066, 0.00072, 
0.00085, 0.00075, 0.00089, 0.00101, 0.0011, 0.00109, 0.00141, 
0.00165, 0.00184, 0.00212, 0.00251, 0.00271, 0.00314, 0.00351, 
0.00391, 0.00462, 0.00489, 0.00554, 0.00579, 0.00607, 0.00698, 
0.00772, 0.00845, 0.00921, 0.00983, 0.01144, 0.01184, 0.01293, 
0.01419, 0.01458, 0.01629, 0.01733, 0.01966, 0.0219, 0.0232, 
0.02579, 0.02815, 0.03018, 0.03336, 0.03783, 0.04089, 0.04888, 
0.05426, 0.0607, 0.06859, 0.07613, 0.08703, 0.09805, 0.11026, 
0.12374, 0.13751, 0.15522, 0.17022, 0.18978, 0.21214, 0.22629, 
0.24824, 0.26987, 0.29207, 0.31465, 0.3374, 0.36011)), row.names = c(NA, 
-101L), class = c("tbl_df", "tbl", "data.frame"))

# some helpers
a0ak_from_q0 <- function (q0, Sex) {
  # coefs coming from Andreev Kingkade (2015)
  #Andreev, Evgeny M., and W. Ward Kingkade. 
  #"Average age at death in infancy and infant mortality level: 
  # Reconsidering the Coale-Demeny formulas at current levels of low mortality." 
  # Demographic Research 33 (2015): 363-390.
  
  Sex <- rep(Sex, length(q0))
  ifelse(Sex == "m", ifelse(q0 < 0.0226, {
    0.1493 - 2.0367 * q0
  }, ifelse(q0 < 0.0785, {
    0.0244 + 3.4994 * q0
  }, 0.2991)), ifelse(q0 < 0.017, {
    0.149 - 2.0867 * q0
  }, ifelse(q0 < 0.0658, {
    0.0438 + 4.1075 * q0
  }, 0.3141)))
}
qx_to_ax <- function (qx, 
          nx = rep(1, length(qx)), 
          age = 0:(length(qx) - 1), 
          sex = "t", 
          closeout = TRUE) {

  ax <- nx/2
  sex <- substr(tolower(sex), 1, 1)
  sex <- ifelse(sex == "b", "t", sex)
  stopifnot(sex %in% c("m", "f", "t"))
  if (min(age) == 0) {
    if (sex != "t") {
      a0 <- a0ak_from_q0(q0 = qx[1], Sex = sex)
    }
    else {
      a0m <- a0ak_from_q0(q0 = qx[1], Sex = "m")
      a0f <- a0ak_from_q0(q0 = qx[1], Sex = "f")
      a0 <- (a0m + a0f)/2
    }
    ax[1] <- a0
  }
  # if (closeout) {
  #   q_close qx[length(qx)]
  #   ax[length(ax)] <- 1/qx[length(qx)]
  # }
  ax
}

qx_to_e0 <- function(qx, 
                     nx = rep(1, length(qx)), 
                     age = 0:(length(qx) - 1), 
                     sex = "t", 
                     closeout = TRUE){
  ax <- qx_to_ax(qx = qx, 
                 nx = nx, 
                 age = age, 
                 sex = sex, 
                 closeout = FALSE)
  # can be solved from HMD Methods Protocol eq 74
  # https://mortality.org/File/GetDocument/Public/Docs/MethodsProtocolV6.pdf
  mx <- qx / (nx - (nx - ax) * qx)
  mx_to_e0(mx = mx, 
           nx = nx, 
           age = age, 
           sex = sex, 
           closeout = closeout)
}

px_to_e0 <- function(px, 
                     nx = rep(1, length(qx)), 
                     age = 0:(length(qx) - 1), 
                     sex = "t", 
                     closeout = TRUE){
  qx = 1 - px
  qx_to_e0(qx = qx, 
           nx = nx, 
           age = age, 
           sex = sex, 
           closeout = closeout)
}

cc_qx <- horiuchi(qx_to_e0,
         pars1 = dat[["2000"]],
         pars2 = dat[["2019"]],
         age = 0:100,
         sex = "m",
         N = 20)

cc_px <- horiuchi(px_to_e0,
                  pars1 = 1-dat[["2000"]],
                  pars2 = 1-dat[["2019"]],
                  age = 0:100,
                  sex = "m",
                  N = 20)

# show these are equal
all(abs(cc_qx - cc_px) < 1e-12)

# This demonstrates that all-cause decompositions can be either-or. Not sure what "both" would look like. This is
# used for an analogy only in the manuscript.