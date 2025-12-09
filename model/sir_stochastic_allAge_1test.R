library(tidyverse)
library(odin.dust)

# I update odin.dust by force
# remotes::install_github("mrc-ide/odin.dust")

# gen_sir <- odin.dust::odin_dust("model/sir_basic_trial.R")
gen_sir <- odin.dust::odin_dust("model/sir_stochastic_allAge.R")

# Running the SIR model with dust
pars <- list(N_ini = 6.7e7,
             log_A_ini1 = 0.7094744,
             log_A_ini2 = 0.7229443,
             time_shift_1 = 0.0639227346367733, #0.302114578070083, # 0.100043419341372, # 
             # time_shift_2 = 0.3766235,
             beta_0 = 0.018, # 0.0381562615720545, #  # 
             beta_1 = 0.2647984930712, # 0.464821184134391, #  # 
             # beta_2 = 0.58190970,
             # scaled_wane = 0.0682579543,
             # psi = (0.5),
             hypo_sigma_2 = (1),
             log_delta1 = (-4.65219600756188)*1.8,
             log_delta2 = (-4.65219600756188)*1.2,
             sigma_1 = 0.1
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

for (t in seq_len(n_times)) {
  model[ , , t] <- sir_model$run(t)
}
# time <- x[1, 1, ] # because in the position of [1, 1, ] is time
# x <- x[-1, , ] # compile all matrix into 1 huge df, delete time (position [-1, , ])
data <- readRDS("inputs/pmcmc_data_week_allAge_ser1_test_2agegroups.rds") %>% 
  glimpse()

sir_data <- dplyr::bind_rows(
  data %>% 
    dplyr::transmute(
      replicate = 1,
      # steps = time_start+1,
      weekly = seq_along(replicate),
      value = count_s1_1,
      compartment = "data_count_s1_1"
    )
  ,
  data %>% 
    dplyr::transmute(
      replicate = 1,
      # steps = time_start+1,
      weekly = seq_along(replicate),
      value = count_s1_2,
      compartment = "data_count_s1_2"
    )
) %>%
  tidyr::complete(weekly, compartment,
                  fill = list(value = 0)) %>% 
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
                                   index == 2 ~ "A1",
                                   index == 3 ~ "A2",
                                   index == 4 ~ "model_D1",
                                   index == 5 ~ "model_D2",
                                   index == 6 ~ "S",
                                   index == 7 ~ "R",
                                   index == 8 ~ "n_AD1_weekly",
                                   index == 9 ~ "n_AD2_weekly"
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
           # compartment %in% c("S",
           #                    "A1", "A2",
           #                    "model_D1", "model_D2",
           #                    "R"),
           compartment %in% c("model_D1", "model_D2",
                              "data_count_s1_1", "data_count_s1_2"),
           compartment != "Time",
           # compartment %in% c("S")
         )
       ,
       aes(x = yearWeek, y = value,
           group = interaction(compartment,replicate),
           colour = compartment)) +
  geom_line() +
  scale_y_continuous(trans = "log1p") +
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

p1 <- ggplot(incidence_modelled %>% 
               dplyr::filter(
                 compartment %in% c("model_D1", "data_count_s1_1"),
                 compartment != "Time",
               )
             ,
             aes(x = yearWeek, y = value,
                 group = interaction(compartment,replicate),
                 colour = compartment)) +
  geom_line() +
  geom_vline(aes(xintercept = as.Date("2010-04-01"),
                 colour = "PCV13 (April 2010)"),
             linetype = "dashed") +
  scale_x_date(limits = c(as.Date(min(all_dates$yearWeek)), as.Date(max(all_dates$yearWeek))),
               date_breaks = "year",
               date_labels = "%Y") +
  ggtitle("Cases (Aggregated by Week) for age 0-9") +
  xlab("Time") +
  ylab("Number of People") +
  theme_bw() +
  theme(legend.position = c(0.15, 0.85),
        legend.title = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent", color = "transparent"))

p2 <- ggplot(incidence_modelled %>% 
               dplyr::filter(
                 compartment %in% c("model_D2", "data_count_s1_2"),
                 compartment != "Time",
               )
             ,
             aes(x = yearWeek, y = value,
                 group = interaction(compartment,replicate),
                 colour = compartment)) +
  geom_line() +
  scale_x_date(limits = c(as.Date(min(all_dates$yearWeek)), as.Date(max(all_dates$yearWeek))),
               date_breaks = "year",
               date_labels = "%Y") +
  ggtitle("Cases (Aggregated by Week) for age 10+") +
  xlab("Time") +
  ylab("Number of People") +
  theme_bw() +
  theme(legend.position = c(0.15, 0.85),
        legend.title = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent", color = "transparent"))

p_combined <- cowplot::plot_grid(p1, p2,
                                 nrow =2,
                                 labels = c("A", "B"))


print(p_combined)

