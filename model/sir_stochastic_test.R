library(tidyverse)
library(odin.dust)

# I update odin.dust by force
# remotes::install_github("mrc-ide/odin.dust")

# gen_sir <- odin.dust::odin_dust("model/sir_basic_trial.R")
gen_sir <- odin.dust::odin_dust("model/sir_stochastic.R")

# Create contact_matrix 3 demographic groups:
# > 2
# 2-64
# 65+
age.limits = c(0, 2, 65)
N_age <- length(age.limits)

contact_demographic <- socialmixr::contact_matrix(socialmixr::survey(
  socialmixr::polymod$participants, 
  socialmixr::polymod$contacts
  ),
  countries = "United Kingdom",
  age.limits = age.limits,
  symmetric = TRUE
)

transmission <- contact_demographic$matrix /
  rep(contact_demographic$demography$population, each = ncol(contact_demographic$matrix))
transmission

# init_A_ini <- 0.01479864

# Running the SIR model with dust
pars <- list(m = transmission,
             N_ini = contact_demographic$demography$population,
             log_A_ini = c(-4, -4, -4),
             D_ini = c(0, 0, 0),
             R_ini = c(0, 0, 0),
             time_shift_1 = 0.366346711348848,
             time_shift_2 = 0.366346711348848,
             beta_0 = 0.063134635077278,
             beta_1 = 0.161472506104886,
             beta_2 = 0.261472506104886,
             gamma = 0.1,
             scaled_wane = (0.9),
             log_delta = (-4.03893492453891), # will be fitted to logN(-10, 0.7)
             psi = (0.5),
             sigma_2 = (1)
)

n_times <- 4745*4 # 4745 or similar to the number of date range (of the provided data), or try 500 for trial
n_particles <- 15L

sir_model <- gen_sir$new(pars = pars,
                         time = 1,
                         n_particles = n_particles,
                         n_threads = 4L,
                         seed = 1L)

model <- array(NA, dim = c(sir_model$info()$len, n_particles, n_times))

# compartment position check
sir_model$info()
sir_model$info()$index$D

# Beta check
time <- seq(1, n_times, 1)
time_shift <- 70
beta <- pars$beta_0*(1+pars$beta_1*sin(2*pi*(time_shift+time)/365))
max(beta)
min(beta)

# R0 estimation (R0 changes due to seasonality)
R0 <- (beta/(pars$log_delta+pars$sigma_1)) +  ((pars$log_delta)*(beta)) / ((pars$log_delta + 192/(4064*4745))*(pars$sigma_2 + 192/(4064*4745))) # print R0
max(R0)
min(R0)
# plot(time, R0)
# pars$beta_1/(pars$delta) + (pars$qu*pars$delta)/(pars$delta*pars$sigma) # print R0


for (t in seq_len(n_times)) {
  model[ , , t] <- sir_model$run(t)
}
time <- model[1, 1, ] # because in the position of [1, 1, ] is time
# model <- model[-1, , ] # compile all matrix into 1 huge df, delete time (position [-1, , ])
library(tidyverse)
glimpse(model)

## 1. Data Load ################################################################
# Plotting the trajectories
# See gen_sir$new(pars = pars, time = 0, n_particles = 1L)$info()
incidence <- read.csv("inputs/incidence_weekly.csv")

par(mfrow = c(1,3), oma=c(2,3,0,0))
for (i in 1:N_age) {
  par(mar = c(3, 4, 2, 0.5))
  cols <- c(S = "#8c8cd9", A = "darkred", D = "orange", R = "#999966", n_AD_daily = "#cc0099", n_AD_cumul = "green")
  matplot(time, t(model[i + 7 + 5*N_age, , ]), type = "l", # Offset to access numbers in age compartment
          xlab = "", ylab = "", yaxt="none", main = paste0("Age ", contact_demographic$demography$age.group[i]),
          col = cols[["n_AD_daily"]], lty = 1)#, ylim=range(model[-1:-3,,]))
  matlines(time, )
  legend("right", lwd = 1, col = cols, legend = names(cols), bty = "n")
  axis(2, las =2)
}
mtext("Number of individuals", side=2,line=1, outer=T)
mtext("Time", side = 1, line = 0, outer =T)
