# Total population data by age, year for each region
# SOURCE: https://www.nomisweb.co.uk/

library(tidyverse)

# load model result
# (based on sir_stochastic_allAge_2post_pmcmc_picts.R, saved in raw_data)
incidence_modelled <- read.csv("raw_data/incidence_modelled_serotype1.csv") %>% 
  dplyr::filter(compartment == "D") %>% 
  dplyr::mutate(yearWeek = as.Date(yearWeek),
                year = year(yearWeek)) %>% 
  dplyr::left_join(
    # age proportion in England
    dplyr::left_join(
      read.csv("raw_data/nomis_population_long.csv") %>% 
        # region stratification is not needed
        dplyr::group_by(year) %>% 
        dplyr::summarise(PopSize_year = sum(PopSize)) %>% 
        dplyr::ungroup()
      ,
      read.csv("raw_data/nomis_population_long.csv") %>% 
        # region stratification is not needed
        dplyr::group_by(ageGroup6, year) %>% 
        dplyr::summarise(PopSize6 = sum(PopSize)) %>% 
        dplyr::ungroup()
      ,
      by = "year"
    ) %>% 
      dplyr::mutate(
        ageGroup6 = factor(ageGroup6,
                           levels = c("<2", "2-4", "5-14", "15-44", "45-64", "65+")),
        PopProp = round(PopSize6/PopSize_year, 1)
      ) %>% 
      glimpse()
    ,
    by = "year",
    relationship = "many-to-many"
  ) %>% 
  dplyr::mutate(case_modelled = round(value*PopProp, 1)) %>% 
  dplyr::filter(!is.na(ageGroup6)) %>% 
  dplyr::arrange(yearWeek) %>% 
  glimpse()

ggplot(incidence_modelled,
       aes(x = yearWeek, y = case_modelled,
         group = ageGroup6,
         colour = ageGroup6)) +
  geom_line() +
  theme_bw()


# comparison with the real case counts data
serotype1_data <- read.csv("raw_data/serotype1_UKHSA_imperial_date_age_region_MOLIS_sequenced_postThesis_cleaned.csv") %>% 
  glimpse()

all_week <- data.frame(week_date = seq.Date(from = min(as.Date(serotype1_data$Earliest.specimen.date)),
                                            to = max(as.Date(serotype1_data$Earliest.specimen.date)), 
                                            by = 1)) %>% 
  dplyr::mutate(week_step = 1:nrow(.),
                iso_week = paste0(year(week_date), "-W", sprintf("%02d", week(week_date)), "-1"),
                yearWeek =ISOweek::ISOweek2date(iso_week)
  ) %>% 
  distinct(yearWeek) %>% 
  glimpse()

serotype1_data_ageGroup6 <- dplyr::left_join(
  all_week
  ,
  serotype1_data %>% 
    dplyr::mutate(Earliest.specimen.date = as.Date(Earliest.specimen.date),
                  iso_week = paste0(year(Earliest.specimen.date), "-W", sprintf("%02d", week(Earliest.specimen.date)), "-1"),
                  yearWeek =ISOweek::ISOweek2date(iso_week),
                  ageGroup6 = case_when(
                    ageGroup7 == "15-30" | ageGroup7 == "31-44" ~ "15-44",
                    TRUE ~ ageGroup7
                  )
    ) %>% 
    dplyr::group_by(yearWeek, ageGroup6) %>% 
    dplyr::summarise(case_data = n()) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(
      case_data = as.numeric(case_data),
      ageGroup6 = case_when(
        is.na(ageGroup6) ~ "Unknown",
        TRUE ~ ageGroup6
      )
    ) %>% 
    dplyr::arrange(yearWeek)
  ,
  by = "yearWeek"
) %>% 
  dplyr::filter(ageGroup6 != "Unknown") %>% 
  glimpse()

ggplot(serotype1_data_ageGroup6,
       aes(x = yearWeek, y = case_data,
           group = ageGroup6,
           colour = ageGroup6)) +
  geom_line() +
  theme_bw()

# combine df
png("outputs/genomics/trials_serotype1/trial_100250/figs/model_vs_data.png",
    width = 24, height = 17, unit = "cm", res = 600)
dplyr::bind_rows(
  incidence_modelled %>% 
    dplyr::transmute(
      yearWeek = as.Date(yearWeek),
      ageGroup6 = ageGroup6,
      version = "model",
      counts = case_modelled
    )
  ,
  serotype1_data_ageGroup6 %>% 
    dplyr::transmute(
      yearWeek = as.Date(yearWeek),
      ageGroup6 = ageGroup6,
      version = "data",
      counts = case_data
    )
) %>% 
  dplyr::mutate(
    ageGroup6 = factor(ageGroup6,
                       levels = c("<2", "2-4", "5-14",
                                  "15-44", "45-64", "65+")),
    version = factor(version,
                     levels = c("model", "data"))
  ) %>% 
  ggplot(aes(x = yearWeek, y = counts,
             color = version, group = version)) +
  geom_line() +
  labs(
    title = "Serotype 1: model vs data"
  ) +
  facet_wrap(~ ageGroup6, 
             # scales = "free_y"
             ) +
  theme_bw()
dev.off()
  
