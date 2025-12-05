# Total population data by age, year for each region
# SOURCE: https://www.nomisweb.co.uk/

library(tidyverse)

# compartment separation for 1 #################################################
age_proportion <- dplyr::left_join(
  read.csv("inputs/nomis_population_long.csv") %>% 
    # region stratification is not needed
    dplyr::group_by(year) %>% 
    dplyr::summarise(PopSize_year = sum(PopSize)) %>% 
    dplyr::ungroup()
  ,
  read.csv("inputs/nomis_population_long.csv") %>% 
    # region stratification is not needed
    dplyr::mutate(
      Age = as.numeric(Age),
      ageGroup_s1 = case_when(
        Age < 10 ~ "0-9",
        Age >= 10 ~ "10+"
      )
    ) %>% 
    dplyr::group_by(ageGroup_s1, year) %>% 
    dplyr::summarise(PopSize_s1 = sum(PopSize)) %>% 
    dplyr::ungroup()
  ,
  by = "year"
) %>% 
  dplyr::mutate(
    ageGroup_s1 = factor(ageGroup_s1,
                         levels = c("0-9", "10+")),
    PopProp = round(PopSize_s1/PopSize_year, 2)
  ) %>% 
  # view() %>% 
  glimpse()

# younger people (0-9) was consistently 12% of the total population


# function for HPC run #########################################################
age_validation <- function(n_sts){
  dir_name <- paste0("outputs/genomics/trial_", n_sts, "/")
  dir.create(paste0(dir_name, "/figs"), FALSE, TRUE)
  
  # load model result
  # (based on sir_stochastic_allAge_2post_pmcmc_picts.R, saved in outputs)
  incidence_modelled <-
    read.csv(paste0(dir_name, "incidence_modelled_serotype1.csv")) %>% 
    dplyr::filter(compartment == "model_D") %>% 
    dplyr::mutate(yearWeek = as.Date(yearWeek),
                  year = year(yearWeek)) %>% 
    dplyr::left_join(
      # age proportion in England
      dplyr::left_join(
        read.csv("inputs/nomis_population_long.csv") %>% 
          # region stratification is not needed
          dplyr::group_by(year) %>% 
          dplyr::summarise(PopSize_year = sum(PopSize)) %>% 
          dplyr::ungroup()
        ,
        read.csv("inputs/nomis_population_long.csv") %>% 
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
    dplyr::arrange(yearWeek)
  
  
  # comparison with the real case counts data
  serotype1_data <- read.csv("inputs/serotype1_UKHSA_imperial_date_age_region_MOLIS_sequenced_postThesis_cleaned.csv") %>% 
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
      tidyr::complete(yearWeek, ageGroup6,
                      fill = list(case_data = 0)) %>% 
      
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
  
  # combine df
  png(paste0(dir_name, "figs/model_vs_data.png"),
      width = 24, height = 17, unit = "cm", res = 600)
  p <- dplyr::bind_rows(
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
    geom_vline(aes(xintercept = as.Date("2010-04-01"),
                   colour = "PCV13 (April 2010)"),
               linetype = "dashed") +
    labs(
      title = "Serotype 1: model vs data"
    ) +
    facet_wrap(~ ageGroup6, 
               # scales = "free_y"
    ) +
    theme_bw()+
    theme(legend.position = "bottom")
  
  print(p)
  dev.off()
  
}

# a slight modification for college's HPC
args <- commandArgs(trailingOnly = T)
n_sts <- as.numeric(args[which(args == "--n_steps") + 1])

age_validation(n_sts)