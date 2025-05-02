library(tidyverse)
library(odin.dust)

# I update odin.dust by force
# remotes::install_github("mrc-ide/odin.dust")

# gen_sir <- odin.dust::odin_dust("model/sir_basic_trial.R")
gen_sir <- odin.dust::odin_dust("model/sir_stochastic_allAge.R")

# Running the SIR model with dust
pars <- list(N_ini = 6.7e7,
             log_A_ini = -1,
             D_ini = 0,
             R_ini = 0,
             time_shift_1 = 0.265185074455071,
             time_shift_2 = 0.2688027206357,
             beta_0 = 0.00904100678898,
             beta_1 = 0.193999573638097,
             beta_2 = 0.184928540835887,
             scaled_wane = (0.486156008428636),
             psi = (0.5),
             sigma_2 = (0.10738841030217),
             log_delta = (-2),
             alpha = 100,
             gamma_annual = 0.01,
             nu_annual = 0.01
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
# sir_model$update_state(pars = pars,
#                        time = 0) # make sure time is 0

# all_date <- incidence$day
# all_date <- data.frame(col = integer(4745))
# incidence <- read.csv("inputs/incidence_week_12F_allAge.csv") %>% 
#   dplyr::mutate(day = week*7) 

n_times <- 200 # 500 for trial
n_particles <- 15L
model <- array(NA, dim = c(sir_model$info()$len, n_particles, n_times))

# R0 estimation
R0 <- pars$beta_0/pars$sigma_2
R0

for (t in seq_len(n_times)) {
  model[ , , t] <- sir_model$run(t)
}
# time <- x[1, 1, ] # because in the position of [1, 1, ] is time
# x <- x[-1, , ] # compile all matrix into 1 huge df, delete time (position [-1, , ])

daily_incidence_modelled <- 
  reshape2::melt(model) %>% 
  dplyr::rename(index = Var1,     # Var1 = dimension that stored SADR values
                replicate = Var2, # Var2 = particles
                time = Var3       # Var3 = time
  ) %>% 
  # dplyr::filter(index < 5) %>% 
  dplyr::mutate(compartment = 
                  dplyr::case_when(index == 1 ~ "Time",
                                   index == 2 ~ "A",
                                   index == 3 ~ "D",
                                   index == 4 ~ "S",
                                   index == 5 ~ "R",
                                   index == 6 ~ "n_AD_weekly",
                                   index == 7 ~ "Ne",
                                   index == 8 ~ "cases_55",
                                   index == 9 ~ "cases_non55",
                                   index == 10 ~ "cases_12F"
                  )) %>% 
  dplyr::rename(value_model = value) %>% 
  dplyr::select(-index) %>% 
  glimpse()

# ggplot(daily_incidence_modelled,
#        aes(x = time, y = value_model,
#            group = interaction(compartment,replicate),
#            colour = compartment)) +
#   geom_line() +
#   scale_y_continuous(trans = "log1p") +
#   ggtitle("Cases (Aggregated by Days)") +
#   xlab("Time") +
#   ylab("Number of People") +
#   theme_bw()

ggplot(daily_incidence_modelled %>% 
         dplyr::filter(grepl("cases|D|R", compartment)),
       aes(x = time, y = value_model,
           group = interaction(compartment,replicate),
           colour = compartment)) +
  geom_line() +
  scale_y_continuous(trans = "log1p") +
  ggtitle("Cases (Aggregated by Days)") +
  xlab("Time") +
  ylab("Number of People") +
  theme_bw()




















# Some viz
par(mfrow = c(1, 1), mar = c(5.1, 5.1, 0.5, 0.5), mgp = c(3.5, 1, 0), las = 1)
cols <- c(S = "#8c8cd9", A = "darkred", D = "#999966", R = "green", n_AD_weekly = "maroon")
matplot(t(x[1, , ]), t(x[4, , ]), type = "l",
        xlab = "t(x[1, , ])", ylab = "Number of individuals",
        col = cols[["S"]], lty = 1, ylim = c(0,pars$N_ini),
        main = "All model")
matlines(t(x[2, , ]), t(x[3, , ]), col = cols[["A"]], lty = 1)
matlines(t(x[3, , ]), t(x[4, , ]), col = cols[["D"]], lty = 1)
# matlines(incidence)
matlines(t(x[5, , ]), t(x[5, , ]), col = cols[["R"]], lty = 1)
matlines(t(x[6, , ]), t(x[5, , ]), col = cols[["n_AD_weekly"]], lty = 1)

legend("right", lwd = 1, col = cols, legend = names(cols), bty = "n")
