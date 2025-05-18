dat_c <- read.csv("raw_data/12F_Jan_2025_combined_cleaned.csv") %>% 
  dplyr::filter(ageGroup3 != "Unknown") %>% 
  dplyr::mutate(week_date = as.Date(week_date),
                iso_week = paste0(year(week_date), "-W", sprintf("%02d", week(week_date)), "-1"),
                yearWeek =ISOweek::ISOweek2date(iso_week)
  ) %>% 
  dplyr::group_by(yearWeek) %>% 
  dplyr::summarise(count_serotype = sum(counts)) %>% 
  dplyr::ungroup() %>% 
  glimpse()

# load genomic data
gen <- read.csv("raw_data/genomic_data_cleaned.csv") %>% 
  dplyr::filter(!is.na(strain),
                collection_date >= as.Date("2017-08-01")) %>%  # after 2017-08-01
  glimpse()

earlier_ne_df <- read.csv("raw_data/GPSC55_mlesky_cleaned_interpolated_predictedModel_binom.csv") %>% 
  glimpse()

all_GPSC55 <- dplyr::bind_rows(
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
  ,
  earlier_ne_df %>% 
    dplyr::select(yearWeek, predicted_count_GPSC55) %>% 
    dplyr::mutate(yearWeek = as.Date(yearWeek)) %>% 
    dplyr::filter(yearWeek <= as.Date("2017-08-01")) %>% 
    dplyr::rename(count_WGS_GPSC55 = predicted_count_GPSC55)
) %>% 
  dplyr::arrange(yearWeek) %>% 
  glimpse()

rand_ss_counts <- data.frame(
  date = as.Date(c("2001-04-01", "2008-04-01", "2012-04-01", "2015-04-01")),
  child_26 = c(0, 0, 0, 4),
  child_32 = c(1, 3, 5, 1),
  child_55 = c(0, 0, 5, 24),
  adult_26 = c(0, 0, 0, 2),
  adult_32 = c(3, 3, 2, 1),
  adult_55 = c(0, 0, 0, 10)
) %>% 
  dplyr::mutate(
    child_non_55 = child_26 + child_32,
    adult_non_55 = adult_26 + adult_32,
    
    child_total = child_non_55 + child_55,
    adult_total = adult_non_55 + adult_55,
    
    total = child_total + adult_total
  ) %>% 
  dplyr::transmute(
    date = date,
    # assume normally distributed across ages; count/total
    prop_WGS_GPSC55_1 = (child_55*(2/15))/total,
    prop_WGS_GPSC55_2 = ((child_55*(13/15)) + (adult_55*(50/71)))/total, # 15–85
    prop_WGS_GPSC55_3 = (adult_55*(21/71))/total
    
  ) %>% 
  glimpse()

# all recorded dates
new_dates <- seq(from = as.Date(min(earlier_ne_df$yearWeek)),
                 to = as.Date(max(earlier_ne_df$yearWeek)),
                 by = "week") %>% 
  glimpse()

fn_1 <- splinefun(x = as.Date(rand_ss_counts$date),
                  y = rand_ss_counts$prop_WGS_GPSC55_1, method = "natural")
fn_2 <- splinefun(x = as.Date(rand_ss_counts$date),
                  y = rand_ss_counts$prop_WGS_GPSC55_2, method = "natural")
fn_3 <- splinefun(x = as.Date(rand_ss_counts$date),
                  y = rand_ss_counts$prop_WGS_GPSC55_3, method = "natural")

all_GPSC55_ageGroup3 <- dplyr::bind_rows(
  earlier_ne_df %>% 
    dplyr::transmute(
      yearWeek = as.Date(yearWeek),
      p_1 = pmax(0, fn_1(yearWeek)),
      p_2 = pmax(0, fn_2(yearWeek)),
      p_3 = pmax(0, fn_3(yearWeek)),
      
      # proportions; p/total_p
      prop_1 = p_1 / (p_1 + p_2 + p_3),
      prop_2 = p_2 / (p_1 + p_2 + p_3),
      prop_3 = p_3 / (p_1 + p_2 + p_3),
      
      # counts
      count_WGS_GPSC55_1 = prop_1*predicted_count_GPSC55,
      count_WGS_GPSC55_2 = prop_2*predicted_count_GPSC55,
      count_WGS_GPSC55_3 = prop_3*predicted_count_GPSC55,
    ) %>% 
    dplyr::select(
      yearWeek, contains("count_WGS")
    )
  ,
  gen %>% 
    dplyr::filter(strain == "GPSC55") %>% 
    dplyr::mutate(week_date = as.Date(week_date),
                  iso_week = paste0(year(week_date), "-W", sprintf("%02d", week(week_date)), "-1"),
                  yearWeek =ISOweek::ISOweek2date(iso_week)
    ) %>% 
    dplyr::group_by(yearWeek, ageGroup3) %>% 
    dplyr::summarise(count = n()) %>% 
    dplyr::ungroup() %>% 
    tidyr::pivot_wider(
      id_cols = yearWeek,
      names_from = ageGroup3,
      values_from = count
    ) %>% 
    dplyr::rename(
      count_WGS_GPSC55_1 = `<2`,
      count_WGS_GPSC55_2 = `2-64`,
      count_WGS_GPSC55_3 = `65+`
    ) %>% 
    dplyr::mutate(
      count_WGS_GPSC55_1 = as.numeric(count_WGS_GPSC55_1),
      count_WGS_GPSC55_2 = as.numeric(count_WGS_GPSC55_2),
      count_WGS_GPSC55_3 = as.numeric(count_WGS_GPSC55_3)
    )
) %>% 
  # test combine all GPSC55
  dplyr::full_join(
    all_GPSC55,
    by = "yearWeek",
    relationship = "many-to-many"
  ) %>% 
  glimpse()


# test viz GPSC55 stratified by agegroups
ggplot(all_GPSC55_ageGroup3
       , aes(x = yearWeek)) +
  geom_line(aes(y = count_WGS_GPSC55_1), colour = "maroon") +
  geom_line(aes(y = count_WGS_GPSC55_2), colour = "gold2") +
  geom_line(aes(y = count_WGS_GPSC55_3), colour = "darkgreen") +
  geom_line(aes(y = count_WGS_GPSC55), colour = "black") +
  geom_vline(xintercept = as.Date("2017-08-01"), color = "steelblue", linetype = "dashed") +
  scale_x_date(limits = c(as.Date("2012-01-01"), as.Date("2022-06-01")), 
               date_breaks = "1 year",
               date_labels = "%Y") +
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




















earlier_ne_df %>% 
  dplyr::select(yearWeek, predicted_count_GPSC55) %>% 
  dplyr::mutate(yearWeek = as.Date(yearWeek)) %>% 
  dplyr::filter(yearWeek <= as.Date("2017-08-01")) %>% 
  dplyr::rename(count_WGS_GPSC55 = predicted_count_GPSC55) %>% 
  # dplyr::mutate(
  #   count_WGS_GPSC55_1 = ifelse(yearWeek >= as.Date("2010-01-01"), count_WGS_GPSC55/3, ),
  # ) %>% 
  glimpse()




dplyr::bind_rows(
  gen %>% 
    dplyr::filter(strain == "GPSC55") %>% 
    dplyr::mutate(week_date = as.Date(week_date),
                  iso_week = paste0(year(week_date), "-W", sprintf("%02d", week(week_date)), "-1"),
                  yearWeek =ISOweek::ISOweek2date(iso_week)
    ) %>% 
    dplyr::group_by(yearWeek, ageGroup3) %>% 
    dplyr::summarise(count_WGS_GPSC55 = n(), .groups = "drop") %>% 
    dplyr::ungroup() %>% 
    tidyr::pivot_wider(
      names_from = ageGroup3,
      values_from = count_WGS_GPSC55
    ) %>% 
    dplyr::rename(
      count_WGS_GPSC55_1 = `<2`,
      count_WGS_GPSC55_2 = `2-64`,
      count_WGS_GPSC55_3 = `65+`
    )
  ,
  earlier_ne_df %>% 
    dplyr::select(yearWeek, predicted_count_GPSC55) %>% 
    dplyr::mutate(yearWeek = as.Date(yearWeek)) %>% 
    dplyr::filter(yearWeek <= as.Date("2017-08-01")) %>% 
    dplyr::rename(count_WGS_GPSC55 = predicted_count_GPSC55)
) %>% 
  glimpse()

