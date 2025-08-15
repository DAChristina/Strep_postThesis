# seasonality check of 12F and GPSC55; pre- and post-2017

library(tidyverse)

# load gen data
gen <-  read.csv("raw_data/gen_lw/genomic_epi_12F.csv") %>% 
  dplyr::mutate(strain = if_else(GPSC == 55, "GPSC55", "non55"),
                ageGroup3 = case_when(
                  Age < 2 ~ "<2",
                  Age >= 2 & Age < 65 ~ "2-64",
                  Age >= 65 ~ "65+",
                  is.na(Age) ~ "Unknown"
                ),
                ageGroup6 = case_when(
                  Age < 2 ~ "<2",
                  Age >= 2 & Age < 5 ~ "2-4",
                  Age >= 5 & Age < 15 ~ "5-14",
                  Age >= 15 & Age < 45 ~ "15-44",
                  Age >= 45 & Age < 65 ~ "45-64",
                  Age >= 65 ~ "65+",
                  is.na(Age) ~ "Unknown"
                ),
                collection_date = lubridate::dmy(collection_date),
                month = lubridate::month(collection_date),
                # year = if_else(month < 7, YEAR + 0.5, YEAR + 1.5), # this is epi year
                year = lubridate::year(collection_date),
                month = lubridate::month(collection_date),
                yearMonth = as.Date(paste0(format(collection_date, "%Y-%m"), "-01")), # year-month as date
                week_date = lubridate::floor_date(collection_date, "week")
  ) %>% 
  dplyr::select(-ngsid, -Age, -YEAR) %>% 
  # dplyr::filter(!is.na(strain), year > 2018) %>% 
  # view() %>% 
  dplyr::filter(strain == "GPSC55") %>% 
  dplyr::mutate(week_date = as.Date(week_date),
                iso_week = paste0(year(week_date), "-W", sprintf("%02d", week(week_date)), "-1"),
                yearWeek = ISOweek::ISOweek2date(iso_week)
  ) %>% 
  dplyr::group_by(yearMonth) %>% 
  dplyr::summarise(count_GPSC55 = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(yearMonth = as.Date(yearMonth)) %>% 
  glimpse()

ggplot(gen %>% 
         dplyr::filter(yearMonth >= as.Date("2017-08-01"))
       ,
       aes(x = yearMonth)) +
  # geom_line(aes(y = count_12F, colour = "12F")) +
  geom_line(aes(y = count_GPSC55, colour = "GPSC55")) +
  # geom_line(aes(y = itr_Ne, colour = "Ne"), size = 1.5) +
  # geom_line(aes(y = centred_Ne, colour = "Ne (centred)"), size = 1.5) +
  geom_vline(xintercept = as.Date("2016-01-01"), color = "steelblue", linetype = "dashed") +
  geom_vline(xintercept = as.Date("2017-09-01"), color = "steelblue", linetype = "dashed") +
  scale_colour_manual(values = c("12F" = "steelblue",
                                 "GPSC55" = "maroon",
                                 "Ne" = "gold2",
                                 "Ne (centred)" = "orange")) +
  # geom_rect(data = blocks, aes(xmin = xmin, xmax = xmax,
  #                              ymin = 0, ymax = Inf),
  #           fill = "darkred",
  #           inherit.aes = F, alpha = 0.2
  # ) +
  theme_bw() +
  scale_y_log10() +
  scale_x_date(date_breaks = "1 month",
               date_labels = "%m") +
  labs(title = "Data for Model Inference") +
  theme(legend.position = c(0.2, 0.85),
        legend.title = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent", color = "transparent"))

serotype_level <- read.csv("raw_data/12F_Jan_2025_combined_cleaned.csv") %>% 
  dplyr::mutate(week_date = as.Date(week_date),
                iso_week = paste0(year(week_date), "-W", sprintf("%02d", week(week_date)), "-1"),
                yearWeek =ISOweek::ISOweek2date(iso_week),
                
                month = lubridate::month(week_date),
                # year = if_else(month < 7, YEAR + 0.5, YEAR + 1.5), # this is epi year
                year = lubridate::year(week_date),
                month = lubridate::month(week_date),
                yearMonth = as.Date(paste0(format(week_date, "%Y-%m"), "-01")), # year-month as date
                week_date = lubridate::floor_date(week_date, "week")
  ) %>% 
  dplyr::group_by(yearMonth) %>% 
  dplyr::summarise(count_12F = sum(counts)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(yearMonth = as.Date(yearMonth),
                year = lubridate::year(yearMonth),) %>% 
  glimpse()

blocks <- data.frame(
  xmin = as.Date(paste0((min(serotype_level$year)-1):(max(serotype_level$year)-1), "-12-01")),  # Start December
  xmax = as.Date(paste0((min(serotype_level$year)):(max(serotype_level$year)), "-03-01")),  # End March
  year = (min(serotype_level$year) - 1):(max(serotype_level$year) - 1)
) %>% 
  glimpse()

ggplot(serotype_level #%>% 
         # dplyr::filter(yearMonth >= as.Date("2017-08-01"))
       ,
       aes(x = yearMonth)) +
  geom_line(aes(y = count_12F, colour = "12F")) +
  # geom_line(aes(y = count_GPSC55, colour = "GPSC55")) +
  # geom_line(aes(y = itr_Ne, colour = "Ne"), size = 1.5) +
  # geom_line(aes(y = centred_Ne, colour = "Ne (centred)"), size = 1.5) +
  geom_vline(xintercept = as.Date("2016-01-01"), color = "steelblue", linetype = "dashed") +
  geom_vline(xintercept = as.Date("2017-09-01"), color = "steelblue", linetype = "dashed") +
  scale_colour_manual(values = c("12F" = "steelblue",
                                 "GPSC55" = "maroon",
                                 "Ne" = "gold2",
                                 "Ne (centred)" = "orange")) +
  geom_rect(data = blocks, aes(xmin = xmin, xmax = xmax,
                               ymin = 0, ymax = Inf),
            fill = "darkred",
            inherit.aes = F, alpha = 0.2
  ) +
  theme_bw() +
  scale_y_log10() +
  scale_x_date(date_breaks = "1 year",
               date_labels = "%Y",
               limits = c(as.Date("2002-12-01"), as.Date("2023-01-01"))) +
  labs(title = "Data for Model Inference") +
  theme(legend.position = c(0.2, 0.85),
        legend.title = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent", color = "transparent"))



################################################################################
# 29 July 2025 discussion:
# If you plot the pre-2020 data by proportion of annual cases occurring in each week,
# we could fit a more sensible seasonality function using something like a Fourier transform I think
# (Plot for all years on top of one another, with x axis as week of the year)

data <- readRDS("inputs/pmcmc_data_week_allAge.rds") %>% 
  dplyr::mutate(year = year(yearWeek)) %>% 
  dplyr::group_by(year) %>% 
  dplyr::summarise(annual_count_55 = sum(count_WGS_GPSC55, na.rm = T), .groups = "drop") %>% 
  dplyr::right_join(
    readRDS("inputs/pmcmc_data_week_allAge.rds") %>% 
      dplyr::mutate(
        year = year(yearWeek)
      )
    ,
    by = "year"
  ) %>% 
  dplyr::mutate(
    year = year(yearWeek),
    yearWeek_count = week(yearWeek),
    week_shifted = ((yearWeek_count-26)%%52)+1,
    weekly_annual_proportion_55 = count_WGS_GPSC55/annual_count_55) %>% 
  dplyr::filter(year <= 2020,
                # !is.na(weekly_annual_proportion_55)
                ) %>% 
  dplyr::arrange(year, week_shifted) %>%
  # view() %>% 
  glimpse()

data_max <- data %>% 
  dplyr::slice_max(weekly_annual_proportion_55, n = 5) %>% 
  glimpse()

# counts per yearWeek test
ggplot(data #%>% 
         #filter(year != 2020)
       ,
       aes(x = yearWeek, y = count_WGS_GPSC55,
           # color = factor(year)
           )) +
  geom_line(size = 1) +
  theme_bw() +
  # scale_y_log10() +
  theme(legend.position = c(0.1, 0.85),
        legend.title = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent", color = "transparent"))

# count
ggplot(data #%>% 
         #filter(year != 2020)
       ,
       aes(x = week_shifted, y = count_WGS_GPSC55, color = factor(year))) +
  geom_line(size = 1) +
  theme_bw() +
  # scale_y_log10() +
  scale_x_continuous(breaks = 0:54,
                     # limits = c(0, 15)
  ) +
  theme(legend.position = c(0.1, 0.85),
        legend.title = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent", color = "transparent"))

# proportion
ggplot(data #%>% 
         #filter(year != 2020)
       ,
       aes(x = week_shifted, y = weekly_annual_proportion_55, color = factor(year))) +
  geom_line(size = 1) +
  theme_bw() +
  # scale_y_log10() +
  scale_x_continuous(breaks = 0:54,
                     # limits = c(0, 15)
                     ) +
  theme(legend.position = c(0.1, 0.85),
        legend.title = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent", color = "transparent"))


data <- data %>% 
  dplyr::filter(!is.na(weekly_annual_proportion_55))

fft_result <- fft(data$weekly_annual_proportion_55)
amplitudes <- Mod(fft_result)
length(amplitudes)
phases <- Arg(fft_result)
phases

reconstructed <- Re(fft(fft_result, inverse = TRUE) / length(data$weekly_annual_proportion_55))
plot(reconstructed, type = "l")  # Now the peak should be shifted!

dominant_freq_index <- which.max(amplitudes[2:(n/2)]) + 1  # skip DC
cycle_length <- 1 / (dominant_freq_index / n)
phase <- phases[dominant_freq_index]

cat("Cycle length:", cycle_length, "\nPhase (radians):", phase)
phase_shift_in_weeks <- (phase / (2 * pi)) * cycle_length


n <- length(data$weekly_annual_proportion_55)
shift_weeks <- 26
prop_shifted <- dplyr::lead(data$weekly_annual_proportion_55, shift_weeks)
prop_shifted[is.na(prop_shifted)] <- tail(data$weekly_annual_proportion_55, shift_weeks)

frequencies <- (0:(n-1)) / n
length(frequencies)

# store it as df
fft_df <- data.frame(
  frequency = frequencies,
  # cycle_length = 1/frequency,
  amplitudes = amplitudes
) %>% 
  # dplyr::filter(frequency > 0, frequency <= 0.5) %>% 
  # view() %>% 
  glimpse()

plot(fft_df$frequency[-1], fft_df$amplitudes[-1]^2, type = "h",
     main = "Periodogram of Weekly Data")

# trial max value
fft_df_max5 <- fft_df %>% 
  dplyr::slice_max(amplitudes, n = 5) %>% 
  dplyr::mutate(
    cycle_length = 1/frequency,
  ) %>% 
  glimpse()


################################################################################
# test dispersion (sort of IQR & variance) for weekly, monthy, season and annual data
data %>% 
  # dplyr::filter(year <= 2019) %>% # omit pandemic era
  dplyr::mutate(week_count = isoweek(yearWeek),
                month = month(yearWeek)
  ) %>% 
  dplyr::group_by(year, week_count) %>% 
  dplyr::summarise(
    count = sum(count_WGS_GPSC55, na.rm = TRUE), .groups = "drop",
    # median = median(count_WGS_GPSC55, na.rm = TRUE),
    # mean = mean(count_WGS_GPSC55, na.rm = TRUE),
    # q1 = quantile(count_WGS_GPSC55, 0.25, na.rm = TRUE),
    # q3 = quantile(count_WGS_GPSC55, 0.75, na.rm = TRUE),
    # variance = var(count_WGS_GPSC55, na.rm = TRUE),
    # VMR = variance/mean
  ) %>% 
  ggplot(.,
         aes(x = factor(week_count), y = count)) +
  geom_boxplot(fill = "skyblue") +
  theme_bw()

data %>% 
  dplyr::mutate(week_count = isoweek(yearWeek),
                month = month(yearWeek)
  ) %>% 
  dplyr::group_by(year, month) %>% 
  dplyr::summarise(
    count = sum(count_WGS_GPSC55, na.rm = TRUE), .groups = "drop",
  ) %>% 
  ggplot(.,
         aes(x = factor(month), y = count)) +
  geom_boxplot(fill = "skyblue") +
  theme_bw()

data %>% 
  dplyr::mutate(week_count = isoweek(yearWeek),
                month = month(yearWeek),
                season = case_when(
                  week_count >= 9  & week_count <= 21 ~ "spring",
                  week_count >= 22 & week_count <= 34 ~ "summer",
                  week_count >= 35 & week_count <= 47 ~ "autumn",
                  week_count >= 48 | week_count <= 8  ~ "winter"
                )
  ) %>% 
  dplyr::group_by(year, season) %>% 
  dplyr::summarise(
    count = sum(count_WGS_GPSC55, na.rm = TRUE), .groups = "drop",
  ) %>% 
  ggplot(.,
         aes(x = factor(season,
                        levels = c("spring", "summer", "autumn", "winter")),
             y = count)) +
  geom_boxplot(fill = "skyblue") +
  theme_bw()

data %>% 
  dplyr::group_by(year) %>% 
  dplyr::summarise(
    count = sum(count_WGS_GPSC55, na.rm = TRUE), .groups = "drop",
  ) %>% 
  ggplot(.,
         aes(x = factor(year), y = count)) +
  geom_boxplot(fill = "skyblue") +
  theme_bw()



