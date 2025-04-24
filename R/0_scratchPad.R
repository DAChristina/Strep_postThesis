# test spline

ne_55 <- read.csv("raw_data/GPSC55_mlesky_cleaned.csv") %>% 
  # dplyr::filter(date >= min(dat_c$week_date)) %>% 
  dplyr::mutate(date = as.Date(date)) %>% 
  # view() %>% 
  glimpse()

ne_55_aveYear <- ne_55 %>% 
  dplyr::group_by(year) %>% 
  dplyr::summarise(Ne = mean(Ne),
                   Ne_lb = mean(Ne_lb),
                   Ne_ub = mean(Ne_ub)) %>% 
  mutate(date = ymd(paste0(year, "-01-01"))) %>% 
  glimpse()

weekly_data <- ne_55_aveYear %>%
  rowwise() %>%
  mutate(
    date_seq = list(seq(date, date + years(1) - days(1), by = "week")),
    Ne = Ne / length(date_seq),
    Ne_lb = Ne_lb / length(date_seq),
    Ne_ub = Ne_ub / length(date_seq)
  ) %>%
  unnest(date_seq) %>%
  # select(date = date_seq, value = value_per_week) %>%
  mutate(type = "Weekly") %>% 
  glimpse()



# ne_55_distributed <- ne_55 %>% 
#   dplyr::mutate(date = as.Date(date)) %>%
#   dplyr::rowwise() %>%
#   dplyr::mutate(
#     days = list(seq.Date(date, date + years(1) - days(1), by = "day")),
#     day_Ne = Ne / length(days),
#     day_Ne_lb = Ne_lb / length(days),
#     day_Ne_ub = Ne_ub / length(days)
#   ) %>%
#   tidyr::unnest(c(days)) %>%
#   dplyr::rename(date_days = days) %>%
#   glimpse()

# test distribute weekly
ne_55_distributed <- ne_55 %>% 
  dplyr::mutate(date = as.Date(date)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    weeks = list(seq.Date(date, date + years(1) - days(1), by = "week")),
    day_Ne = Ne / length(weeks),
    day_Ne_lb = Ne_lb / length(weeks),
    day_Ne_ub = Ne_ub / length(weeks)
  ) %>%
  # tidyr::unnest(c(weeks)) %>%
  dplyr::rename(date_weeks = weeks) %>% 
  distinct(date_weeks, .keep_all = T) %>% 
  glimpse()

ggplot(ne_55_distributed, aes(x = date, y = day_Ne)) +
  geom_line() +
  geom_vline(xintercept = as.Date("2017-07-01"), color = "steelblue", linetype = "dashed") +
  # geom_label(aes(x = as.Date("2017-07-01"), y = 40, label = "Jul-17"),
  #            fill = "white", color = "black") +
  # geom_point(data = test, aes(x = yearWeek, y = day_Ne), color = "steelblue", size = 1) +
  geom_ribbon(data = ne_55_distributed,
              aes(x = date, ymin = day_Ne_lb, ymax = day_Ne_ub),
              inherit.aes = FALSE,
              fill = "steelblue", alpha = 0.3
  ) +
  scale_x_date(limits = c(min(as.Date(gen$collection_date)), max(as.Date(dat_c$week_date))),
               date_breaks = "1 year",
               date_labels = "%Y") +
  theme_bw()


# test weekly, monthly & annually
# ne_55_monthly <- ne_55_distributed %>%
#   dplyr::group_by(yearMonth) %>%
#   dplyr::summarise(day_Ne = sum(day_Ne), .groups = "drop") %>% 
#   glimpse()
# 
ne_55_weekly <- weekly_data %>%
  dplyr::mutate(
    iso_week = paste0(year(date), "-W", sprintf("%02d", week(date)), "-1"),
    yearWeek =ISOweek::ISOweek2date(iso_week)
  ) %>%
  # dplyr::group_by(yearWeek) %>%
  # dplyr::summarise(week_Ne = sum(day_Ne), .groups = "drop") %>%
  glimpse()

# test combine interpolated data
test <- allAges_weekly %>% 
  select(-contains("Ne")) %>% 
  dplyr::left_join(
    # ne_55_weekly
    ne_55_weekly %>%
      dplyr::mutate(
        iso_week = paste0(year(date), "-W", sprintf("%02d", week(date)), "-1"),
        yearWeek =ISOweek::ISOweek2date(iso_week)
      )
    ,
    by = c("yearWeek")
  ) %>% 
  glimpse()



ggplot(test, aes(x = yearWeek, y = count, color = type.x)) +
  geom_line() +
  geom_vline(xintercept = as.Date("2017-07-01"), color = "steelblue", linetype = "dashed") +
  # geom_label(aes(x = as.Date("2017-07-01"), y = 40, label = "Jul-17"),
  #            fill = "white", color = "black") +
  geom_point(data = test, aes(x = yearWeek, y = Ne), color = "steelblue", size = 1) +
  geom_ribbon(data = test,
              aes(x = yearWeek, ymin = Ne_lb, ymax = Ne_ub),
              inherit.aes = FALSE,
              fill = "steelblue", alpha = 0.3
  ) +
  scale_x_date(limits = c(min(as.Date(gen$collection_date)), max(as.Date(dat_c$week_date))),
               date_breaks = "1 year",
               date_labels = "%Y") +
  theme_bw() +
  theme(legend.position = c(0.1, 0.85),
        legend.title = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent", color = "transparent"))


# report 1






# report 2











































new_dates <- seq(from = as.Date(min(ne_55$date)), to = as.Date(max(ne_55$date)), by = "day") %>% 
  glimpse()

spline_fit_Ne <- stats::spline(x = as.numeric(ne_55$date),
                               y = ne_55$Ne,
                               xout = as.numeric(new_dates))$y
spline_fit_Ne_lo <- stats::spline(x = as.numeric(ne_55$date),
                               y = ne_55$Ne_lb,
                               xout = as.numeric(new_dates))$y
spline_fit_Ne_up <- stats::spline(x = as.numeric(ne_55$date),
                               y = ne_55$Ne_ub,
                               xout = as.numeric(new_dates))$y

interpolated_df <- data.frame(
  date = as.Date(new_dates),
  itr_Ne = spline_fit_Ne
) %>% 
  dplyr::full_join(
    data.frame(
      date = as.Date(new_dates),
      itr_Ne_lo = spline_fit_Ne_lo
    )
      ,
      by = "date"
    ) %>% 
  dplyr::full_join(
    data.frame(
      date = as.Date(new_dates),
      itr_Ne_up = spline_fit_Ne_up
    ),
    by = "date"
  ) %>% 
  dplyr::mutate(
    iso_week = paste0(year(date), "-W", sprintf("%02d", week(date)), "-1"),
    yearWeek =ISOweek::ISOweek2date(iso_week)
  ) %>% 
  glimpse()

ggplot() +
  geom_point(data = interpolated_df, aes(x = date, y = itr_Ne), color = "steelblue", size = 1) +
  geom_ribbon(data = interpolated_df,
              aes(x = date, ymin = itr_Ne_lo, ymax = itr_Ne_up),
              inherit.aes = FALSE,
              fill = "steelblue", alpha = 0.3
  ) +
  theme_bw()






