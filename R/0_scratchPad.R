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

################################################################################
# report 1
weekly_data <- ne_55_aveYear %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    # weekly seq
    weekly_seq = list(seq(date, date + years(1) - days(1), by = "week")),
    Ne = Ne / length(weekly_seq),
    Ne_lb = Ne_lb / length(weekly_seq),
    Ne_ub = Ne_ub / length(weekly_seq)
  ) %>%
  tidyr::unnest(weekly_seq) %>%
  dplyr::mutate(
    iso_week = paste0(year(weekly_seq), "-W", sprintf("%02d", week(weekly_seq)), "-1"),
    yearWeek = ISOweek::ISOweek2date(iso_week)
  ) %>% 
  ungroup() %>%
  glimpse()
plot(weekly_data$yearWeek, weekly_data$Ne)

# test combine interpolated data
test <- allAges_weekly %>% 
  select(-contains("Ne")) %>% 
  dplyr::left_join(
    # ne_55_weekly
    weekly_data %>%
      dplyr::select(yearWeek, contains("Ne"))
    ,
    by = c("yearWeek")
  ) %>% 
  glimpse()

ggplot(test, aes(x = yearWeek, y = count, color = type)) +
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


################################################################################
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

p1 <- ggplot() +
  geom_point(data = interpolated_df, aes(x = date, y = itr_Ne), color = "steelblue", size = 1) +
  geom_ribbon(data = interpolated_df,
              aes(x = date, ymin = itr_Ne_lo, ymax = itr_Ne_up),
              inherit.aes = FALSE,
              fill = "steelblue", alpha = 0.3
  ) +
  annotate(
    "rect",
    xmin = min(as.Date(gen$collection_date)), xmax = as.Date("2021-03-20"),
    ymin = 0, ymax = max(interpolated_df$itr_Ne)+50, na.rm = TRUE,
    alpha = 0.1, fill = NA, colour = "skyblue"
  ) +
  theme_bw()


p2 <- ggplot() +
  geom_point(data = interpolated_df, aes(x = date, y = itr_Ne), color = "steelblue", size = 1) +
  geom_ribbon(data = interpolated_df,
              aes(x = date, ymin = itr_Ne_lo, ymax = itr_Ne_up),
              inherit.aes = FALSE,
              fill = "steelblue", alpha = 0.3
  ) +
  scale_x_date(limits = c(min(as.Date(gen$collection_date)), as.Date("2021-03-20")),
               date_breaks = "1 year",
               date_labels = "%Y") +
  theme_bw()

cowplot::plot_grid(p1, p2,
                   ncol = 1,
                   labels = c("A", "B"))


