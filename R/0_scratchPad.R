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

# test 
ggplot(test, aes(x = yearWeek, y = count, color = type)) +
  geom_line() +
  geom_vline(xintercept = as.Date("2017-07-01"), color = "steelblue", linetype = "dashed") +
  # geom_label(aes(x = as.Date("2017-07-01"), y = 40, label = "Jul-17"),
  #            fill = "white", color = "black") +
  # geom_point(data = test, aes(x = yearWeek, y = Ne), color = "steelblue", size = 1) +
  # geom_ribbon(data = test,
  #             aes(x = yearWeek, ymin = Ne_lb, ymax = Ne_ub),
  #             inherit.aes = FALSE,
  #             fill = "steelblue", alpha = 0.3
  # ) +
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
new_dates <- seq(from = as.Date(min(ne_55$date)), to = as.Date(max(ne_55$date)), by = "1 week") %>% 
  glimpse()

fn_Ne    <- splinefun(x = numeric_dates, y = ne_55$Ne, method = "natural")
fn_Ne_lo <- splinefun(x = numeric_dates, y = ne_55$Ne_lb, method = "natural")
fn_Ne_up <- splinefun(x = numeric_dates, y = ne_55$Ne_ub, method = "natural")

spline_fit_Ne    <- fn_Ne(numeric_new_dates)
spline_fit_Ne_lo <- fn_Ne_lo(numeric_new_dates)
spline_fit_Ne_up <- fn_Ne_up(numeric_new_dates)

time_diffs_weeks <- as.numeric(diff(new_dates)) / 7
change_Ne <- abs(diff(spline_fit_Ne)) / time_diffs_weeks
change_Ne_lo <- abs(diff(spline_fit_Ne_lo)) / time_diffs_weeks
change_Ne_up <- abs(diff(spline_fit_Ne_up)) / time_diffs_weeks

interpolated_df <- tibble(
  date = new_dates[-1],
  itr_Ne = spline_fit_Ne[-1],
  itr_Ne_lo = spline_fit_Ne_lo[-1],
  itr_Ne_up = spline_fit_Ne_up[-1],
  time_diffs_weeks = time_diffs_weeks,
  change_Ne = change_Ne,
  change_Ne_lo = change_Ne_lo,
  change_Ne_up = change_Ne_up
) %>%
  mutate(
    iso_week = paste0(year(date), "-W", sprintf("%02d", week(date)), "-1"),
    yearWeek = ISOweek::ISOweek2date(iso_week)
  ) %>%
  glimpse()

test <- allAges_weekly %>% 
  select(-contains("Ne")) %>% 
  dplyr::left_join(
    interpolated_df %>%
      dplyr::select(yearWeek, contains("Ne"))
    ,
    by = c("yearWeek")
  ) %>% 
  glimpse()

################################################################################
selected_GPSC55 <- test %>% 
  dplyr::filter(type == "count_WGS55") %>% 
  glimpse()

model <- mgcv::gam(count ~ s(change_Ne), data = selected_GPSC55)
earlier_ne_df <- interpolated_df %>%
  dplyr::filter(yearWeek < as.Date("2018-07-01")) %>% # midpoint; they started intensively sequenced GPSC55 ("2017-07-01")
  glimpse()

earlier_ne_df <- earlier_ne_df %>% 
  dplyr::mutate(predicted_GPSC55 = predict(model, newdata = earlier_ne_df)) %>% 
  glimpse()

# test change_Ne vs. WGS_GPSC55 counts
ggplot(selected_GPSC55, aes(x = change_Ne, y = count)) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x), se = TRUE) +
  theme_bw()

# test GPSC55 previous WGS
combined <- dplyr::bind_rows(
  selected_GPSC55 %>% 
    filter(type == "count_WGS55") %>% 
    select(yearWeek, count) %>%
    mutate(source = "GPSC55 WGS data")
  ,
  earlier_ne_df %>% 
    select(yearWeek, predicted_GPSC55) %>%
    rename(count = predicted_GPSC55) %>% 
    mutate(source = "GAM model")
) %>% 
  # weird array conversion
  dplyr::mutate(
    count = as.data.frame(count)
    ) %>% 
  unnest(cols = count) %>% 
  glimpse()


ggplot(combined, aes(x = yearWeek, y = count, color = source)) +
  geom_line() +
  # geom_point(size = 0.5, alpha = 0.6) +
  labs(
    color = "source"
  ) +
  scale_x_date(limits = c(min(as.Date(dat_c$week_date)), max(as.Date(dat_c$week_date))),
               date_breaks = "1 year",
               date_labels = "%Y") +
  theme_bw() +
  theme(legend.position = c(0.1, 0.85),
        legend.title = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent", color = "transparent"))






test <- allAges_weekly %>% 
  select(-contains("Ne")) %>% 
  dplyr::left_join(
    interpolated_df %>%
      dplyr::select(yearWeek, contains("Ne"))
    ,
    by = c("yearWeek")
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

ggplot(test, aes(x = yearWeek, y = count, color = type)) +
  geom_line() +
  geom_vline(xintercept = as.Date("2017-07-01"), color = "steelblue", linetype = "dashed") +
  # geom_label(aes(x = as.Date("2017-07-01"), y = 40, label = "Jul-17"),
  #            fill = "white", color = "black") +
  geom_point(data = test, aes(x = yearWeek, y = change_Ne), color = "steelblue", size = 0.5) +
  # geom_ribbon(data = test,
  #             aes(x = yearWeek, ymin = change_Ne_lo, ymax = change_Ne_up),
  #             inherit.aes = FALSE,
  #             fill = "steelblue", alpha = 0.3
  # ) +
  scale_x_date(limits = c(min(as.Date(dat_c$week_date)), max(as.Date(dat_c$week_date))),
               date_breaks = "1 year",
               date_labels = "%Y") +
  theme_bw() +
  theme(legend.position = c(0.1, 0.85),
        legend.title = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent", color = "transparent"))

