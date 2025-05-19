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
  dplyr::mutate(
    ageGroup2 = case_when(
      ageGroup6 %in% c("<2", "2-4", "5-14") ~ "child",
      ageGroup6 %in% c("15-44", "45-64", "65+") ~ "adult",
      T ~ NA
    )
  ) %>% 
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

# for rand_ss_counts I've decided to add data points from GPSC55 sampling period;
# 4 pre-GPSC55 data points sampling period was randomly selected & shared by NC
# focused only for GPSC55
rand_ss_counts <- dplyr::bind_rows(
  data.frame(
    date = as.Date(c("2001-04-01", "2008-04-01", "2012-04-01", "2015-04-01")), # random date
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
      
      total_55 = child_55 + adult_55
    ) %>% 
    dplyr::transmute(
      yearWeek = date,
      # 3 age groups assume normally distributed across ages; count/total
      prop_WGS_GPSC55_1 = (child_55*(2/15))/total_55,
      prop_WGS_GPSC55_2 = ((child_55*(13/15)) + (adult_55*(50/71)))/total_55, # 15–85
      prop_WGS_GPSC55_3 = (adult_55*(21/71))/total_55,
      # test view counts
      count_WGS_GPSC55_1 = child_55*(2/15),
      count_WGS_GPSC55_2 = (child_55*(13/15)) + (adult_55*(50/71)),
      count_WGS_GPSC55_3 = adult_55*(21/71),
      
    )
  ,
  dplyr::full_join(
    gen %>% 
      dplyr::filter(strain == "GPSC55") %>% 
      dplyr::mutate(week_date = as.Date(week_date),
                    iso_week = paste0(year(week_date), "-W", sprintf("%02d", week(week_date)), "-1"),
                    yearWeek =ISOweek::ISOweek2date(iso_week)
      ) %>% 
      # 3 ageGroups
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
    ,
    gen %>% 
      dplyr::filter(strain == "GPSC55") %>% 
      dplyr::mutate(week_date = as.Date(week_date),
                    iso_week = paste0(year(week_date), "-W", sprintf("%02d", week(week_date)), "-1"),
                    yearWeek =ISOweek::ISOweek2date(iso_week)
      ) %>% 
      dplyr::group_by(yearWeek) %>% 
      dplyr::summarise(count_55 = n()) %>% 
      dplyr::ungroup()
    ,
    by = "yearWeek"
  ) %>% 
    dplyr::transmute(
      yearWeek = as.Date(yearWeek),
      prop_WGS_GPSC55_1 = count_WGS_GPSC55_1/count_55,
      prop_WGS_GPSC55_2 = count_WGS_GPSC55_2/count_55,
      prop_WGS_GPSC55_3 = count_WGS_GPSC55_3/count_55,
      
      # test view counts
      count_WGS_GPSC55_1 = as.numeric(count_WGS_GPSC55_1),
      count_WGS_GPSC55_2 = as.numeric(count_WGS_GPSC55_2),
      count_WGS_GPSC55_3 = as.numeric(count_WGS_GPSC55_3)
      
    )
) %>% 
  glimpse()

# sanity check for proportions
count <- ggplot(rand_ss_counts
                , aes(x = yearWeek)) +
  geom_line(aes(y = count_WGS_GPSC55_1), colour = "maroon") +
  geom_line(aes(y = count_WGS_GPSC55_2), colour = "orange") +
  geom_line(aes(y = count_WGS_GPSC55_3), colour = "darkgreen") +
  # geom_line(aes(y = prop_WGS_GPSC55), colour = "black") +
  geom_vline(xintercept = as.Date("2017-08-01"), color = "steelblue", linetype = "dashed") +
  scale_x_date(limits = c(as.Date("2001-01-01"), as.Date("2022-06-01")), 
               date_breaks = "1 year",
               date_labels = "%Y") +
  theme_bw() +
  labs(
    title = "GPSC55 Counts + Real Data (3 age groups)",
    y = "GPSC55 counts"
  ) +
  theme(legend.position = c(0.15, 0.85),
        legend.title = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent", colour = "transparent"))


# all recorded dates
new_dates <- seq(from = as.Date(min(earlier_ne_df$yearWeek)),
                 to = as.Date(max(earlier_ne_df$yearWeek)),
                 by = "week") %>% 
  glimpse()

fn_1 <- splinefun(x = as.Date(rand_ss_counts$yearWeek),
                  y = rand_ss_counts$prop_WGS_GPSC55_1, method = "natural")
fn_2 <- splinefun(x = as.Date(rand_ss_counts$yearWeek),
                  y = rand_ss_counts$prop_WGS_GPSC55_2, method = "natural")
fn_3 <- splinefun(x = as.Date(rand_ss_counts$yearWeek),
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
pred <- ggplot(all_GPSC55_ageGroup3
               , aes(x = yearWeek)) +
  geom_line(aes(y = count_WGS_GPSC55_1), colour = "maroon") +
  geom_line(aes(y = count_WGS_GPSC55_2), colour = "orange") +
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


cowplot::plot_grid(count, pred,
                   ncol = 1,
                   labels = c("A", "B"))




# test 2 ageGroups
rand_ss_counts <- dplyr::bind_rows(
  data.frame(
    date = as.Date(c("2001-04-01", "2008-04-01", "2012-04-01", "2015-04-01")), # random date
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
      
      total_55 = child_55 + adult_55
    ) %>% 
    dplyr::transmute(
      yearWeek = date,
      # test 2 ageGroups
      prop_WGS_GPSC55_1 = child_55/total_55,
      prop_WGS_GPSC55_2 = adult_55/total_55,
      # test view counts
      count_WGS_GPSC55_1 = child_55,
      count_WGS_GPSC55_2 = adult_55,
      
    )
  ,
  dplyr::full_join(
    gen %>% 
      dplyr::filter(strain == "GPSC55") %>% 
      dplyr::mutate(week_date = as.Date(week_date),
                    iso_week = paste0(year(week_date), "-W", sprintf("%02d", week(week_date)), "-1"),
                    yearWeek =ISOweek::ISOweek2date(iso_week)
      ) %>% 
      # 2 ageGroups
      dplyr::group_by(yearWeek, ageGroup2) %>% 
      dplyr::summarise(count = n()) %>% 
      dplyr::ungroup() %>% 
      tidyr::pivot_wider(
        id_cols = yearWeek,
        names_from = ageGroup2,
        values_from = count
      ) %>% 
      dplyr::rename(
        count_WGS_GPSC55_1 = `child`,
        count_WGS_GPSC55_2 = `adult`
      ) %>% 
      dplyr::mutate(
        count_WGS_GPSC55_1 = as.numeric(count_WGS_GPSC55_1),
        count_WGS_GPSC55_2 = as.numeric(count_WGS_GPSC55_2)
      )
    ,
    gen %>% 
      dplyr::filter(strain == "GPSC55") %>% 
      dplyr::mutate(week_date = as.Date(week_date),
                    iso_week = paste0(year(week_date), "-W", sprintf("%02d", week(week_date)), "-1"),
                    yearWeek =ISOweek::ISOweek2date(iso_week)
      ) %>% 
      dplyr::group_by(yearWeek) %>% 
      dplyr::summarise(count_55 = n()) %>% 
      dplyr::ungroup()
    ,
    by = "yearWeek"
  ) %>% 
    dplyr::transmute(
      yearWeek = as.Date(yearWeek),
      prop_WGS_GPSC55_1 = count_WGS_GPSC55_1/count_55,
      prop_WGS_GPSC55_2 = count_WGS_GPSC55_2/count_55,
      
      # test view counts
      count_WGS_GPSC55_1 = as.numeric(count_WGS_GPSC55_1),
      count_WGS_GPSC55_2 = as.numeric(count_WGS_GPSC55_2)
      
    )
) %>% 
  glimpse()

# sanity check for proportions
count <- ggplot(rand_ss_counts
                , aes(x = yearWeek)) +
  geom_line(aes(y = count_WGS_GPSC55_1), colour = "maroon") +
  geom_line(aes(y = count_WGS_GPSC55_2), colour = "orange") +
  # geom_line(aes(y = count_WGS_GPSC55_3), colour = "darkgreen") +
  # geom_line(aes(y = prop_WGS_GPSC55), colour = "black") +
  geom_vline(xintercept = as.Date("2017-08-01"), color = "steelblue", linetype = "dashed") +
  scale_x_date(limits = c(as.Date("2001-01-01"), as.Date("2022-06-01")), 
               date_breaks = "1 year",
               date_labels = "%Y") +
  theme_bw() +
  labs(
    title = "GPSC55 Counts + Real Data (2 age groups)",
    y = "GPSC55 counts"
  ) +
  theme(legend.position = c(0.15, 0.85),
        legend.title = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent", colour = "transparent"))



# all recorded dates
new_dates <- seq(from = as.Date(min(earlier_ne_df$yearWeek)),
                 to = as.Date(max(earlier_ne_df$yearWeek)),
                 by = "week") %>% 
  glimpse()

fn_1 <- splinefun(x = as.Date(rand_ss_counts$yearWeek),
                  y = rand_ss_counts$prop_WGS_GPSC55_1, method = "natural")
fn_2 <- splinefun(x = as.Date(rand_ss_counts$yearWeek),
                  y = rand_ss_counts$prop_WGS_GPSC55_2, method = "natural")

all_GPSC55_ageGroup2 <- dplyr::bind_rows(
  earlier_ne_df %>% 
    dplyr::transmute(
      yearWeek = as.Date(yearWeek),
      p_1 = pmax(0, fn_1(yearWeek)),
      p_2 = pmax(0, fn_2(yearWeek)),
      
      # proportions; p/total_p
      prop_1 = p_1 / (p_1 + p_2),
      prop_2 = p_2 / (p_1 + p_2),
      
      # counts
      count_WGS_GPSC55_1 = prop_1*predicted_count_GPSC55,
      count_WGS_GPSC55_2 = prop_2*predicted_count_GPSC55,
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
    dplyr::group_by(yearWeek, ageGroup2) %>% 
    dplyr::summarise(count = n()) %>% 
    dplyr::ungroup() %>% 
    tidyr::pivot_wider(
      id_cols = yearWeek,
      names_from = ageGroup2,
      values_from = count
    ) %>% 
    dplyr::rename(
      count_WGS_GPSC55_1 = `child`,
      count_WGS_GPSC55_2 = `adult`
    ) %>% 
    dplyr::mutate(
      count_WGS_GPSC55_1 = as.numeric(count_WGS_GPSC55_1),
      count_WGS_GPSC55_2 = as.numeric(count_WGS_GPSC55_2)
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
pred <- ggplot(all_GPSC55_ageGroup2
               , aes(x = yearWeek)) +
  geom_line(aes(y = count_WGS_GPSC55_1), colour = "maroon") +
  geom_line(aes(y = count_WGS_GPSC55_2), colour = "orange") +
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


cowplot::plot_grid(count, pred,
                   ncol = 1,
                   labels = c("A", "B"))





