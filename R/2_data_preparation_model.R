rm(list = ls())
library(tidyverse)

if (!dir.exists("inputs")) {
  dir.create("inputs")
}

# Data preparation for serotype 1
# mcstate data preparation #####################################################
# non-heterogeneity (allAges), weekly
allAges_weekly_ser1 <- read.csv("raw_data/serotype1_UKHSA_imperial_date_age_region_MOLIS_sequenced_postThesis_cleaned.csv") %>% 
  dplyr::mutate(Earliest.specimen.date = as.Date(Earliest.specimen.date),
                iso_week = paste0(year(Earliest.specimen.date), "-W", sprintf("%02d", week(Earliest.specimen.date)), "-1"),
                yearWeek =ISOweek::ISOweek2date(iso_week)
  ) %>% 
  dplyr::group_by(yearWeek) %>% 
  dplyr::summarise(count_serotype = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(
    count_serotype = as.numeric(count_serotype),
  ) %>% 
  dplyr::arrange(yearWeek) %>% 
  dplyr::mutate(
    yearWeek = as.Date(yearWeek),
    day = as.numeric(round((yearWeek - (as.Date("2003-01-01")-2)))), # min(dat_g$Earliest.specimen.date)-2 to make it 7
    # day = seq_len(n())
  ) %>% 
  dplyr::filter(day > 0) %>% 
  mcstate::particle_filter_data(.,
                                time = "day", # I use steps instead of day
                                rate = 1, # I change the model to weekly, therefore weekly rate is required
                                initial_time = 0
  ) %>%
  glimpse()

saveRDS(allAges_weekly_ser1, "inputs/pmcmc_data_week_allAge_ser1.rds")

# test plot
ggplot(allAges_weekly_ser1
       , aes(x = yearWeek, y = count_serotype)) +
  geom_line(size = 1) +
  geom_vline(xintercept = as.Date("2017-08-01"), color = "steelblue", linetype = "dashed") +
  scale_x_date(limits = c(as.Date("2002-12-31"), as.Date("2023-12-31")),
               date_breaks = "1 year",
               date_labels = "%Y") +
  # scale_y_log10() +
  theme_bw() +
  labs(
    title = "Serotype 1 Counts",
    y = "Serotype 1 counts"
  ) +
  theme(legend.position = c(0.15, 0.85),
        legend.title = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent", colour = "transparent"))










