dir_1 <- paste0("outputs/genomics/trial_5e+05_GPSC55_FINAL/")
dir_name <- paste0("outputs/genomics/trial_5e+05_GPSC55_FINAL/")

dir.create(paste0(dir_name, "/figs"), FALSE, TRUE)
# run 4_post_pmcmc_pics.R first
results <- read.csv(paste0(dir_name, "tune_initial_with_CI.csv")) %>% 
  glimpse()

# Running the SIR model with dust
pars <- list(N_ini = 6.7e7,
             log_A_ini = results[1,2],
             time_shift_1 = results[2,2],
             beta_0 = results[3,2],
             beta_1 = results[4,2],
             log_delta1 = results[5,2],
             log_delta2 = results[6,2],
             
             # fixed parameters
             sigma_1 = 1/28,
             sigma_2 = 1,
             mu_0 = 1/(80.70*365)
)

time <- seq(0, 365*15, 1) # from 2010 to 2025, 15 years
beta <- pars$beta_0*((1+pars$beta_1*cos(2*pi*((pars$time_shift_1*365)+time)/365)))

plot(time, beta)

print(c(max(beta), min(beta))) # beta simulated
((max(beta)-min(beta))/2)+min(beta) # average beta ~ pars$beta_0
pars$beta_0

# R0 estimation (R0 changes due to seasonality)
R0 <- beta*(10^(pars$log_delta1)+10^(pars$log_delta2)+pars$sigma_2+pars$mu_0)/
  ((10^(pars$log_delta1)+10^(pars$log_delta2)+pars$sigma_1+pars$mu_0)*
     (pars$sigma_2+pars$mu_0))

R0_ave <- pars$beta_0*(10^(pars$log_delta1)+10^(pars$log_delta2)+pars$sigma_2+pars$mu_0)/
  ((10^(pars$log_delta1)+10^(pars$log_delta2)+pars$sigma_1+pars$mu_0)*
     (pars$sigma_2+pars$mu_0))

print(c(max(R0), min(R0))) # R0 simulated
((max(R0)-min(R0))/2)+min(R0) # average R0 ~ R0_ave

R0_ave

# save beta & R0 simulation
beta_R0_1 <- data.frame(
  time = seq(0, 365*15, 1), # from 2010 to 2025, 15 years,
  beta_t = pars$beta_0*((1+pars$beta_1*
                           cos(2*pi*((pars$time_shift_1*365)+time)/365))),
  #
  R0 = beta*(10^(pars$log_delta1)+10^(pars$log_delta2)+pars$sigma_2+pars$mu_0)/
    ((10^(pars$log_delta1)+10^(pars$log_delta2)+pars$sigma_1+pars$mu_0)*
       (pars$sigma_2+pars$mu_0)),
  R_ave = pars$beta_0*(10^(pars$log_delta1)+10^(pars$log_delta2)+pars$sigma_2+pars$mu_0)/
    ((10^(pars$log_delta1)+10^(pars$log_delta2)+pars$sigma_1+pars$mu_0)*
       (pars$sigma_2+pars$mu_0))
) %>% 
  dplyr::mutate(
    date = as.Date("2010-01-01") + days(time - 1),
    weekly = ceiling(time/7)) %>% 
  glimpse()

# png("pictures/R0_weekly_simulated4745times4.png", width = 24, height = 12, unit = "cm", res = 1200)
ggplot(beta_R0_df, aes(x = date, y = R0,
)) +
  geom_line() + 
  scale_x_date(date_breaks = "1 month",
               date_labels = "%m",
               limits = as.Date(c("2010-01-01", "2011-12-31"))) +
  # xlim(as.Date('1/1/2010'), as.Date('1/1/2012'), format="%m%Y")) +
  ggtitle("GPSC55 Cases (Aggregated by Weeks)") +
  xlab("Time (in Week)") +
  ylab("R0") +
  theme_bw()
# dev.off()




