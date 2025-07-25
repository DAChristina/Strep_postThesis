library(tidyverse)
library(odin.dust)

# I update odin.dust by force
# remotes::install_github("mrc-ide/odin.dust")

# gen_sir <- odin.dust::odin_dust("model/sir_basic_trial.R")
gen_sir <- odin.dust::odin_dust("model/sir_stochastic_allAge.R")

# Running the SIR model with dust
pars <- list(N_ini = 6.7e7,
             log_A_ini = -4.0763973707201,
             time_shift_1 = 0.0639227346367733, #0.302114578070083, # 0.100043419341372, # 
             # time_shift_2 = 0.3766235,
             beta_0 = 0.0365789436634438, # 0.0381562615720545, #  # 
             beta_1 = 0.2647984930712, # 0.464821184134391, #  # 
             # beta_2 = 0.58190970,
             # scaled_wane = 0.0682579543,
             # psi = (0.5),
             hypo_sigma_2 = (1),
             log_delta = (-4.65219600756188) # (-4.65135010884371) #  # 
             # alpha = 0.01,
             # gamma_annual = 0.01,
             # nu_annual = 0.01
)

# time_points <- round(seq(0, by = (365/52), length.out = 52*3)) # per-week, 22 years
# n_times <- length(time_points)
n_times <- 5000
n_pars <- 1L
sir_model <- gen_sir$new(pars = pars,
                         time = 1,
                         n_particles = n_pars,
                         n_threads = n_pars,
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
model <- array(NA, dim = c(sir_model$info()$len, n_pars, n_times))

# R0 estimation
R0 <- pars$beta_0/pars$sigma_2
R0

for (t in seq_len(n_times)) {
  model[ , , t] <- sir_model$run(t)
}
# time <- x[1, 1, ] # because in the position of [1, 1, ] is time
# x <- x[-1, , ] # compile all matrix into 1 huge df, delete time (position [-1, , ])
data <- readRDS("inputs/pmcmc_data_week_allAge.rds") %>% 
  glimpse()

sir_data <- data %>% 
  dplyr::transmute(
    replicate = 1,
    # steps = time_start+1,
    weekly = seq_along(replicate),
    value = count_WGS_GPSC55,
    compartment = "data_count_WGS_GPSC55"
  ) %>%
  glimpse()

# all_dates <- data.frame(date = seq(min(data$yearWeek), max(data$yearWeek), by = "day")) %>%
#   dplyr::mutate(
#     steps = seq_along(date)
#   ) %>%
#   glimpse()
all_dates <- data %>%
  dplyr::select(yearWeek) %>% 
  dplyr::mutate(
    weekly = seq_along(yearWeek)
  ) %>%
  glimpse()

# focused on n_AD_weekly (already in weeks)
incidence_modelled <- 
  reshape2::melt(model) %>% 
  dplyr::rename(index = Var1,     # Var1 = dimension that stored SADR values
                replicate = Var2, # Var2 = particles
                steps = Var3       # Var3 = steps are in days, but n_AD_weekly is aggregated in weeks
  ) %>% 
  # dplyr::filter(index < 5) %>%
  dplyr::mutate(compartment = 
                  dplyr::case_when(index == 1 ~ "Time",
                                   index == 2 ~ "A",
                                   index == 3 ~ "D",
                                   index == 4 ~ "S",
                                   index == 5 ~ "R",
                                   index == 6 ~ "model_n_AD_weekly",
                                   index == 7 ~ "Ne",
                                   index == 8 ~ "cases_55",
                                   index == 9 ~ "cases_non55",
                                   index == 10 ~ "cases_12F"
                  )) %>% 
  dplyr::select(-index) %>%
  dplyr::mutate(weekly = ceiling(steps/7)) %>% 
  dplyr::group_by(replicate, weekly, compartment) %>% 
  dplyr::summarise(value = sum(value, na.rm = T),
                   # date = max(date),
                   .groups = "drop") %>% 
  dplyr::ungroup() %>% 
  dplyr::bind_rows(sir_data) %>%
  dplyr::full_join(
    all_dates
    ,
    by = "weekly"
  ) %>%
  # dplyr::filter(date %in% data$yearWeek) %>%
  glimpse()


ggplot(incidence_modelled %>% 
         dplyr::filter(
           # grepl("cases|D|data", compartment),
           # compartment %in% c("D", "model_n_AD_weekly", "data_count_WGS_GPSC55"), # redesign the model, would rather fit to D
           # compartment %in% c("model_n_AD_weekly", "data_count_WGS_GPSC55"),
           # compartment %in% c("D", "n_AD_weekly"),
           compartment %in% c("D", "data_count_WGS_GPSC55"),
           # compartment %in% c("D"),
           compartment != "Time",
           # compartment %in% c("S")
         )
       ,
       aes(x = yearWeek, y = value,
           group = interaction(compartment,replicate),
           colour = compartment)) +
  geom_line() +
  # scale_y_continuous(trans = "log1p") +
  # scale_y_continuous(limits = c(0, 50)) +
  # scale_x_continuous(limits = c(0, 700)) +
  scale_x_date(limits = c(as.Date(min(all_dates$yearWeek)), as.Date(max(all_dates$yearWeek))),
               date_breaks = "year",
               date_labels = "%Y") +
  ggtitle("Cases (Aggregated by Week)") +
  xlab("Time") +
  ylab("Number of People") +
  theme_bw() +
  theme(legend.position = c(0.15, 0.85),
        legend.title = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent", color = "transparent"))


transformations <- data.frame(
  log_beta_0 = seq(-5, 0, 0.5)
) %>% 
  dplyr::mutate(beta_0 = as.numeric(10^log_beta_0),
                scaled_A_ini = seq(0, 1, 0.1),
                max_A_ini = 0,
                min_A_ini = -20,
                log_A_ini = scaled_A_ini*(max_A_ini-min_A_ini)+min_A_ini,
                A_ini = 10^(log_A_ini)*6.7e7) %>% 
  glimpse()

# test A_ini
scaled_A_ini = 0.7484698
# scaled to log
(scaled_A_ini*(0-(-15)+(-15))) # scaled_x * (max - min) + min

log_A_ini = -3.77295312521097
# log to scaled
(log_A_ini-(-15)) / (0-(-15)) # (x - min) / (max - min)









