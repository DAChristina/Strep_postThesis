library(tidyverse)
library(odin.dust)

# I update odin.dust by force
# remotes::install_github("mrc-ide/odin.dust")
source("global/all_function_allAge.R")

# load chains
n_sts <- 1000
dir_name <- paste0("outputs/genomics/trial_", n_sts, "/")
dir.create(paste0(dir_name, "/figs"), FALSE, TRUE)
results <- read.csv(paste0(dir_name, "tune_initial_with_CI.csv"))


# gen_sir <- odin.dust::odin_dust("model/sir_basic_trial.R")
gen_sir <- odin.dust::odin_dust("model/sir_stochastic_allAge.R")

# Running the SIR model with dust
pars <- list(N_ini = 6.7e7,
             log_A_ini = results[1,2],
             D_ini = 0,
             R_ini = 0,
             time_shift_1 = results[2,2],
             time_shift_2 = results[3,2],
             beta_0 = results[4,2],
             beta_1 = results[5,2],
             beta_2 = results[6,2],
             scaled_wane = results[7,2],
             # psi = (0.5),
             # sigma_2 = (0.10738841030217),
             log_delta = results[8,2],
             alpha = results[9,2],
             gamma_annual = results[10,2],
             nu_annual = results[11,2]
)

time_points <- round(seq(0, by = (365/52), length.out = 52*3)) # per-week, 22 years
# n_times <- length(time_points)
n_times <- 13800 # roughly from 1987-2025 in days
sir_model <- gen_sir$new(pars = pars,
                         time = 1,
                         n_particles = 1L,
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
n_particles <- 1L
model <- array(NA, dim = c(sir_model$info()$len, n_particles, n_times))

# R0 estimation
R0 <- pars$beta_0/pars$sigma_2
R0

for (t in seq_len(n_times)) {
  model[ , , t] <- sir_model$run(t)
}
# time <- x[1, 1, ] # because in the position of [1, 1, ] is time
# x <- x[-1, , ] # compile all matrix into 1 huge df, delete time (position [-1, , ])

sir_data <- readRDS("inputs/pmcmc_data_week_allAge.rds") %>% 
  tidyr::pivot_longer(
    cols = c(contains("count_"), "Ne"),
    names_to = "compartment",
    values_to = "value"
  ) %>% 
  dplyr::select(yearWeek, compartment, value) %>% 
  glimpse()

allDays <- data.frame(
  date = seq.Date(from = min(sir_data$yearWeek), to = max(sir_data$yearWeek), by = "day")
) %>% 
  dplyr::mutate(
    time_day = seq_along(date),
    iso_week = paste0(year(date), "-W", sprintf("%02d", week(date)), "-1"),
    yearWeek =ISOweek::ISOweek2date(iso_week),
    time_week = dplyr::dense_rank(yearWeek) - 1
  ) %>% 
  glimpse()
  


daily_incidence_modelled <- 
  reshape2::melt(model) %>% 
  dplyr::rename(index = Var1,     # Var1 = dimension that stored SADR values
                replicate = Var2, # Var2 = particles
                time = Var3       # Var3 = time
  ) %>% 
  # dplyr::filter(index < 5) %>% 
  dplyr::mutate(compartment = 
                  dplyr::case_when(index == 1 ~ "model_time",
                                   index == 2 ~ "model_A",
                                   index == 3 ~ "model_D",
                                   index == 4 ~ "model_S",
                                   index == 5 ~ "model_R",
                                   index == 6 ~ "model_n_AD_weekly",
                                   index == 7 ~ "model_Ne",
                                   index == 8 ~ "model_cases_55",
                                   index == 9 ~ "model_cases_non55",
                                   index == 10 ~ "model_cases_12F"
                  )) %>% 
  dplyr::mutate(
    week = floor(time / 7)
  ) %>% 
  dplyr::group_by(
    week, compartment
  ) %>% 
  dplyr::summarise(
    value = sum(value, na.rm = TRUE), .groups = "drop"
  ) %>% 
  # dplyr::select(-index) %>% 
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
         dplyr::filter(grepl("cases|D", compartment)),
       aes(x = week, y = value,
           group = compartment, #interaction(compartment,replicate),
           colour = compartment)) +
  geom_line() +
  scale_y_continuous(trans = "log1p") +
  ggtitle("Cases (Aggregated by Week)") +
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
