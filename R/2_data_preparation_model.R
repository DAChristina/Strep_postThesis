rm(list = ls())

library(tidyverse)

if (!dir.exists("inputs")) {
  dir.create("inputs")
}

# 12F
dat_combined_week <- read.csv("raw_data/12F_Jan_2025_combined_cleaned.csv") %>% 
  dplyr::group_by(week_date) %>% 
  dplyr::summarise(counts = sum(counts)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(week_date = as.Date(week_date))

all_date <- data.frame(week_date = seq.Date(from = min(dat_combined_week$week_date),
                                          to = max(dat_combined_week$week_date), 
                                          by = 1)) %>% 
  dplyr::mutate(week_step = 1:nrow(.))

incidence_week_allAge <- dplyr::left_join(all_date, 
                                          dat_combined_week, 
                                          by = "week_date"
                                          ) %>% 
  dplyr::select(week_step, counts)

write.csv(incidence_week_allAge, "inputs/incidence_week_12F_allAge.csv", row.names = FALSE)

