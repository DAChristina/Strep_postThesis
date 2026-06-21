rm(list = ls())
library(tidyverse)

if (!dir.exists("inputs")) {
  dir.create("inputs")
}

source("global/all_function_allAge.R")
# global/all_function_allAge.R also incorporated:
# burnin_days

# Data preparation for serotype 12F (GPSC55)
# mcstate data preparation #####################################################
# load epidata
dat_c <- read.csv("inputs/12F_Jan_2025_combined_cleaned.csv") %>% 
  dplyr::filter(ageGroup3 != "Unknown") %>% 
  dplyr::mutate(
    ageGroup12F = case_when(
      ageGroup6 == "45-64" | ageGroup6 == "65+" ~ "45+",
      TRUE ~ "0-44"
    )
  ) %>% 
  glimpse()

# load genomic data
gen <- read.csv("raw_data/genomic_data_cleaned.csv") %>% 
  dplyr::filter(!is.na(strain),
                collection_date >= as.Date("2017-08-01")) %>%  # after 2017-08-01
  glimpse()

earlier_ne_df_1 <- read.csv("raw_data/GPSC55_mlesky_cleaned_interpolated_predictedModel_binom_1.csv") %>% 
  glimpse()

earlier_ne_df_2 <- read.csv("raw_data/GPSC55_mlesky_cleaned_interpolated_predictedModel_binom_2.csv") %>% 
  glimpse()

# load ne
# ne_55 <- read.csv("raw_data/GPSC55_mlesky_cleaned.csv") %>% 
#   dplyr::filter(date >= min(dat_c$week_date))
interpolated_ne <- read.csv("raw_data/GPSC55_mlesky_cleaned_interpolated.csv") %>% 
  dplyr::select(itr_Ne, yearWeek) %>% 
  dplyr::rename(Ne = itr_Ne) %>% 
  dplyr::mutate(yearWeek = as.Date(yearWeek)) %>% 
  glimpse()

# non-heterogeneity (allAges), weekly
ageGroup12F_weekly <- dat_c %>% 
  dplyr::mutate(week_date = as.Date(week_date),
                iso_week = paste0(year(week_date), "-W", sprintf("%02d", week(week_date)), "-1"),
                yearWeek =ISOweek::ISOweek2date(iso_week)
  ) %>% 
  dplyr::group_by(yearWeek, ageGroup12F) %>% 
  dplyr::summarise(count_12F = sum(counts)) %>% 
  dplyr::ungroup() %>% 
  tidyr::pivot_wider(
    .,
    names_from = contains("ageGroup"),
    names_prefix = "count_",
    values_from = "count_12F"
  ) %>% 
  dplyr::rename(
    count_12F_1 = "count_0-44",
    count_12F_2 = "count_45+"
  ) %>% 
  dplyr::mutate(yearWeek = as.Date(yearWeek)
  ) %>% 
  dplyr::full_join(
    dplyr::bind_rows(
      gen %>% 
        dplyr::filter(strain == "GPSC55") %>% 
        dplyr::mutate(week_date = as.Date(week_date),
                      iso_week = paste0(year(week_date), "-W", sprintf("%02d", week(week_date)), "-1"),
                      yearWeek =ISOweek::ISOweek2date(iso_week)
        ) %>% 
        dplyr::group_by(yearWeek, ageGroup12F) %>% 
        dplyr::summarise(count_WGS_GPSC55 = n()) %>% 
        dplyr::ungroup() %>% 
        tidyr::pivot_wider(
          .,
          names_from = contains("ageGroup"),
          names_prefix = "count_",
          values_from = "count_WGS_GPSC55"
        ) %>% 
        dplyr::rename(
          count_55_1 = "count_0-44",
          count_55_2 = "count_45+"
        ) %>% 
        dplyr::mutate(yearWeek = as.Date(yearWeek))
      ,
      dplyr::full_join(
        earlier_ne_df_1 %>% 
          dplyr::transmute(
            yearWeek = as.Date(yearWeek),
            count_55_1 = predicted_count_55_1
          ) %>% 
          dplyr::filter(yearWeek <= as.Date("2017-08-01"))
        ,
        earlier_ne_df_2 %>% 
          dplyr::transmute(
            yearWeek = as.Date(yearWeek),
            count_55_2 = predicted_count_55_2
          ) %>% 
          dplyr::filter(yearWeek <= as.Date("2017-08-01"))
        ,
        by = "yearWeek"
        ,
        relationship = "many-to-many"
      )
    ) # %>% 
    # distinct(yearWeek, .keep_all = T)
    , by = "yearWeek"
  ) %>% 
  dplyr::full_join(
    interpolated_ne
    ,
    by = "yearWeek"
    ,
    relationship = "many-to-many"
  ) %>%
  dplyr::mutate(
    count_12F_1 = as.numeric(count_12F_1),
    count_12F_2 = as.numeric(count_12F_2),
    count_55_1 = as.numeric(round(count_55_1), 0),
    count_55_2 = as.numeric(round(count_55_2), 0),
    count_55_all = as.numeric(count_55_1 + count_55_2),
    count_12F_all = as.numeric(count_12F_1 + count_12F_2),
    Ne = as.numeric(Ne)
  ) %>% 
  dplyr::arrange(yearWeek) %>% 
  dplyr::filter(
    yearWeek >= as.Date("2010-01-01") # & yearWeek <= as.Date("2020-01-01") # filter out data not based on initial Ne but the first time GPSC55 was predicted 
  ) %>%
  distinct(yearWeek, .keep_all = T) %>% 
  dplyr::mutate(
    yearWeek = as.Date(yearWeek),
    day = as.numeric(round((yearWeek - as.Date("2010-01-04")))),
    # day = seq_len(n())
  ) %>%
  dplyr::filter(day > 0) %>% 
  mcstate::particle_filter_data(.,
                                time = "day", # I use steps instead of day
                                rate = 1, # I change the model to weekly, therefore weekly rate is required
                                initial_time = 0
  ) %>%
  glimpse()

saveRDS(ageGroup12F_weekly, "inputs/pmcmc_data_week_ageGroup12F.rds")

# test viz combined GPSC55
ageGroup12F_weekly_long <- ageGroup12F_weekly %>% 
  tidyr::pivot_longer(cols = c(count_12F_1,
                               count_12F_2,
                               count_55_1,
                               count_55_2,
                               Ne),
                      names_to = "group",
                      values_to = "count") %>% 
  glimpse()

# plot with Ne
ggplot(ageGroup12F_weekly_long
       , aes(x = yearWeek, y = count, colour = group)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("count_12F_1" = "lightcoral",
                                "count_12F_2" = "maroon",
                                "count_55_1" = "grey40",
                                "count_55_2" = "grey10",
                                "Ne" = "gold2")) +
  geom_vline(xintercept = as.Date("2017-08-01"), color = "steelblue", linetype = "dashed") +
  scale_x_date(limits = c(as.Date("2010-01-01"), as.Date("2022-06-01")), 
               date_breaks = "1 year",
               date_labels = "%Y") +
  scale_y_log10() +
  theme_bw() +
  labs(
    title = "GPSC55 Counts Prediction + Real Data",
    y = "GPSC55 counts (in log10-scaled)"
  ) +
  theme(legend.position = c(0.15, 0.85),
        legend.title = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent", colour = "transparent"))

# plot without Ne
ggplot(ageGroup12F_weekly_long %>% 
         dplyr::filter(group != "Ne")
       , aes(x = yearWeek, y = count, colour = group)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("count_12F_1" = "lightcoral",
                                "count_12F_2" = "maroon",
                                "count_55_1" = "grey40",
                                "count_55_2" = "grey10",
                                "Ne" = "gold2")) +
  geom_vline(xintercept = as.Date("2017-08-01"), color = "steelblue", linetype = "dashed") +
  scale_x_date(limits = c(as.Date("2010-01-01"), as.Date("2022-06-01")), 
               date_breaks = "1 year",
               date_labels = "%Y") +
  # scale_y_log10() +
  theme_bw() +
  labs(
    title = "GPSC55 Counts Prediction + Real Data",
    y = "GPSC55 counts"
  ) +
  theme(legend.position = c(0.15, 0.85),
        legend.title = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent", colour = "transparent"))


# test only available GPSC55 data for non-heterogeneity
gen <- read.csv("raw_data/genomic_data_cleaned.csv") %>% 
  dplyr::filter(!is.na(strain),
                # collection_date >= as.Date("2017-08-01")
  ) %>%  # after 2017-08-01
  glimpse()

edited_data <- dat_c %>% 
  dplyr::mutate(week_date = as.Date(week_date),
                iso_week = paste0(year(week_date), "-W", sprintf("%02d", week(week_date)), "-1"),
                yearWeek =ISOweek::ISOweek2date(iso_week)
  ) %>% 
  dplyr::group_by(yearWeek) %>% 
  dplyr::summarise(count_serotype = sum(counts)) %>% 
  dplyr::ungroup() %>% 
  dplyr::full_join(
    dplyr::bind_rows(
      gen %>% 
        dplyr::filter(strain == "GPSC55") %>% 
        dplyr::mutate(week_date = as.Date(week_date),
                      iso_week = paste0(year(week_date), "-W", sprintf("%02d", week(week_date)), "-1"),
                      yearWeek =ISOweek::ISOweek2date(iso_week)
        ) %>% 
        dplyr::group_by(yearWeek) %>% 
        dplyr::summarise(count_WGS_GPSC55 = n()) %>% 
        dplyr::ungroup() %>% 
        dplyr::mutate(yearWeek = as.Date(yearWeek))
      # ,
      # earlier_ne_df %>% 
      #   dplyr::select(yearWeek, predicted_count_GPSC55) %>% 
      #   dplyr::mutate(yearWeek = as.Date(yearWeek)) %>% 
      #   dplyr::filter(yearWeek <= as.Date("2017-08-01")) %>% 
      #   dplyr::rename(count_WGS_GPSC55 = predicted_count_GPSC55)
    )
    ,
    by = "yearWeek"
  ) %>% 
  dplyr::full_join(
    gen %>% 
      dplyr::filter(strain == "non55") %>% 
      dplyr::mutate(week_date = as.Date(week_date),
                    iso_week = paste0(year(week_date), "-W", sprintf("%02d", week(week_date)), "-1"),
                    yearWeek =ISOweek::ISOweek2date(iso_week)
      ) %>% 
      dplyr::group_by(yearWeek) %>% 
      dplyr::summarise(count_WGS_non55 = n()) %>% 
      dplyr::ungroup()
    ,
    by = "yearWeek"
  ) %>% 
  dplyr::full_join(
    interpolated_ne
    ,
    by = "yearWeek"
    ,
    relationship = "many-to-many"
  ) %>%
  dplyr::mutate(
    count_serotype = as.numeric(count_serotype),
    count_WGS_GPSC55 = as.numeric(count_WGS_GPSC55),
    count_WGS_non55 = as.numeric(count_WGS_non55),
    Ne = as.numeric(Ne)
  ) %>% 
  # tidyr::pivot_longer(
  #   cols = starts_with(c("count_")), # ignore Ne at the moment
  #   names_to = "type",
  #   values_to = "count"
  # ) %>% 
  dplyr::arrange(yearWeek) %>% 
  dplyr::filter(
    yearWeek >= as.Date("2010-01-01") # filter out data not based on initial Ne but the first time GPSC55 was predicted 
  ) %>%
  dplyr::mutate(# count_WGS_GPSC55 = round(count_WGS_GPSC55), # rounded cases
    # count_WGS_GPSC55 = ifelse(count_WGS_GPSC55 == 0, NA_real_, count_WGS_GPSC55),
    yearWeek = as.Date(yearWeek),
    day = as.numeric(round((yearWeek - as.Date("2010-01-04")))),
    # day = seq_len(n())
  ) %>%
  dplyr::filter(day > 0) %>% 
  mcstate::particle_filter_data(.,
                                time = "day", # I use steps instead of day
                                rate = 1, # I change the model to weekly, therefore weekly rate is required
                                initial_time = 0
  ) %>%
  glimpse()

saveRDS(edited_data, "inputs/pmcmc_data_week_allAge_nonGAM.rds")
