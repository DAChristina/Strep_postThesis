dir_name <- paste0("outputs/genomics/trial_", "100250_final_best_vcvmodif_NOreruns_corrected_seasonality", "/")
dir.create(paste0(dir_name, "/figs"), FALSE, TRUE)
# run 4_post_pmcmc_pics.R first
results <- read.csv(paste0(dir_name, "tune_initial_with_CI.csv")) %>% 
  glimpse()


# gen_sir <- odin.dust::odin_dust("model/sir_basic_trial.R")
gen_sir <- odin.dust::odin_dust("model/sir_stochastic_allAge.R")

# Running the SIR model with dust
pars <- list(N_ini = 6.7e7,
             log_A_ini = results[1,2],
             D_ini = 0,
             R_ini = 0,
             time_shift_1 = results[2,2],
             time_shift_2 = results[3,2],
             beta_0 = results[3,2],
             beta_1 = results[4,2],
             beta_2 = results[6,2],
             # scaled_wane = results[7,2],
             # psi = (0.5),
             hypo_sigma_2 = (1),
             log_delta = results[5,2]
             # alpha = results[9,2],
             # gamma_annual = results[10,2],
             # nu_annual = results[11,2]
)

pars_lo_CI <- list(N_ini = 6.7e7,
                   log_A_ini = results[1,3],
                   D_ini = 0,
                   R_ini = 0,
                   time_shift_1 = results[2,3],
                   time_shift_2 = results[3,2],
                   beta_0 = results[3,3],
                   beta_1 = results[4,3],
                   beta_2 = results[6,2],
                   # scaled_wane = results[7,2],
                   # psi = (0.5),
                   hypo_sigma_2 = (1),
                   log_delta = results[5,3]
                   # alpha = results[9,2],
                   # gamma_annual = results[10,2],
                   # nu_annual = results[11,2]
)


pars_hi_CI <- list(N_ini = 6.7e7,
                   log_A_ini = results[1,4],
                   D_ini = 0,
                   R_ini = 0,
                   time_shift_1 = results[2,4],
                   time_shift_2 = results[3,2],
                   beta_0 = results[3,4],
                   beta_1 = results[4,4],
                   beta_2 = results[6,2],
                   # scaled_wane = results[7,2],
                   # psi = (0.5),
                   hypo_sigma_2 = (1),
                   log_delta = results[5,4]
                   # alpha = results[9,2],
                   # gamma_annual = results[10,2],
                   # nu_annual = results[11,2]
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

# betas & R0 estimation
time <- seq(0, 365*22, 1)
beta <- pars$beta_0*((1+pars$beta_1*cos(2*pi*((pars$time_shift_1*365)+time)/365)))

beta_lo_CI <- pars_lo_CI$beta_0*((1+pars_lo_CI$beta_1*cos(2*pi*((pars_lo_CI$time_shift_1*365)+time)/365)))
beta_hi_CI <- pars_hi_CI$beta_0*((1+pars_hi_CI$beta_1*cos(2*pi*((pars_hi_CI$time_shift_1*365)+time)/365)))

print(c(max(beta), min(beta))) # beta simulated with no infant vaccination

# R0 estimation (R0 changes due to seasonality)
mu_0 <- 1/(80.70*365)
mu_1 <- 0
pars$sigma_1 <- 1/28
pars$sigma_2 <- 1/1
R0_vacc <- beta/((mu_0+(10^(pars$log_delta))+pars$sigma_1)*((pars$sigma_2)+(mu_0+mu_1)))
R0_vacc_lo_CI <- beta_lo_CI/((mu_0+(10^(pars_lo_CI$log_delta))+pars$sigma_1)*((pars$sigma_2)+(mu_0+mu_1)))
R0_vacc_hi_CI <- beta_hi_CI/((mu_0+(10^(pars_hi_CI$log_delta))+pars$sigma_1)*((pars$sigma_2)+(mu_0+mu_1)))

# save beta & R0 simulation
beta_R0_df <- data.frame(
  time = time,
  beta = beta,
  beta_hi_CI = beta_hi_CI,
  beta_lo_CI = beta_lo_CI,
  #
  R0_vacc = R0_vacc,
  R0_vacc_hi_CI = R0_vacc_hi_CI,
  R0_vacc_lo_CI = R0_vacc_lo_CI
) %>% 
  dplyr::mutate(
    date = as.Date("2010-01-01") + days(time - 1),
    weekly = ceiling(time/7))

# png("pictures/R0_weekly_simulated4745times4.png", width = 24, height = 12, unit = "cm", res = 1200)
ggplot(beta_R0_df, aes(x = date, y = R0_vacc,
                                                            # group = variable,
                                                            # colour = variable
)) +
  geom_ribbon(aes(ymin = R0_vacc_lo_CI, ymax = R0_vacc_hi_CI), fill = "steelblue", alpha = 0.2) +
  geom_line() + 
  # geom_vline(data = vaccine_UK, aes(xintercept = as.Date(vaccine_UK$date),
  #                                   colour = vaccine),
  #            linetype = "dashed") +
  # scale_color_manual(values = c(col_compartment),
  #                    name = "States",
  #                    breaks = c("Susceptible", "Asymptomatic", "Diseased", "Recovered", "Data_dis", "PCV7", "PCV13"),
  #                    labels = c("Susceptible", "Asymptomatic", "Diseased", "Recovered", "Diseased Data", "PCV7", "PCV13")
  # ) +
  scale_x_date(date_breaks = "1 month",
               date_labels = "%m",
               limits = as.Date(c("2010-01-01", "2012-12-31"))) +
  # xlim(as.Date('1/1/2010'), as.Date('1/1/2012'), format="%m%Y")) +
  ggtitle("Simulation Model for Serotype 1 Cases (Aggregated by Days)") +
  ggtitle("Serotype 1 Cases (Aggregated by Weeks)") +
  xlab("Time (in Week)") +
  ylab("R0") +
  theme_bw()
# dev.off()




