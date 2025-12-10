library(tidyverse)
library(odin.dust)

# I update odin.dust by force
# remotes::install_github("mrc-ide/odin.dust")

# gen_sir <- odin.dust::odin_dust("model/sir_basic_trial.R")
gen_sir <- odin.dust::odin_dust("model/sir_stochastic_ageGroup3.R")

# Running the SIR model with dust
pars <- list(N_ini = c(0.12*6.7e7, 0.88*6.7e7),
             log_A_ini = c(-1, -1),
             D_ini = c(0,0),
             R_ini = c(0,0),
             # trial m matrix failed
             # m = c(0.5, 0.5),
             time_shift_1 = 0.265185074455071,
             time_shift_2 = 0.2688027206357,
             beta_0 = 0.80904100678898,
             beta_1 = 0.193999573638097,
             beta_2 = 0.184928540835887,
             scaled_wane = (0.486156008428636),
             log_delta_kids = (-5.79347000840983), # will be fitted to logN(-10, 0.7)
             log_delta_adults = (-5.79347000840983), # will be fitted to logN(-10, 0.7)
             psi = (0.5),
             sigma_2 = (0.10738841030217)
)

n_times <- round(seq(1, by = 365/52, length.out = 52*3)) # per-week, 22 years
sir_model <- gen_sir$new(pars = pars,
                         time = 1,
                         n_particles = 15L,
                         n_threads = 4L,
                         seed = 1L)

# compartment position check
# sir_model$info()
# sir_model$info()$index$n_AD_weekly
# update_state is required "every single time" to run & produce matrix output (don't know why)
sir_model$update_state(pars = pars,
                       time = 0) # make sure time is 0

# all_date <- incidence$day
# all_date <- data.frame(col = integer(4745))
# incidence <- read.csv("inputs/incidence_week_12F_allAge.csv") %>% 
#   dplyr::mutate(day = week*7) 

n_times <- 200 # 500 for trial
n_particles <- 15L
x <- array(NA, dim = c(sir_model$info()$len, n_particles, n_times))

# R0 estimation
R0 <- pars$beta_0/pars$sigma_2
R0

for (t in seq_len(n_times)) {
  x[ , , t] <- sir_model$run(t)
}
# time <- x[1, 1, ] # because in the position of [1, 1, ] is time
# x <- x[-1, , ] # compile all matrix into 1 huge df, delete time (position [-1, , ])

# Some viz
par(mfrow = c(1, 1), mar = c(5.1, 5.1, 0.5, 0.5), mgp = c(3.5, 1, 0), las = 1)
cols <- c(S = "#8c8cd9", A = "darkred", D = "#999966", R = "green", n_AD_weekly = "maroon")
matplot(t(x[1, , ]), t(x[2, , ]), type = "l",
        xlab = "t(x[1, , ])", ylab = "Number of individuals",
        col = cols[["S"]], lty = 1, ylim = c(0,pars$N_ini),
        main = "All model")
matlines(t(x[1, , ]), t(x[3, , ]), col = cols[["A"]], lty = 1)
matlines(t(x[1, , ]), t(x[4, , ]), col = cols[["D"]], lty = 1)
# matlines(incidence)
matlines(t(x[1, , ]), t(x[5, , ]), col = cols[["R"]], lty = 1)
matlines(t(x[1, , ]), t(x[5, , ]), col = cols[["n_AD_weekly"]], lty = 1)

legend("right", lwd = 1, col = cols, legend = names(cols), bty = "n")
