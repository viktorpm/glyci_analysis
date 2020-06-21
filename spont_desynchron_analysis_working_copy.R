### LOADING PACKAGES ############

library(R.matlab)
library(tidyverse)
library(reshape2)
library(signal)
library(bspec) ### power spectrum
library(WaveletComp) ### wavelet
library(diptest) ### to test distribution uni/multimodality (ISI)
library(ggrepel)
library(exactRankTests) ### wilcox.test for ties
library(ggsignif)
library(lumberjack)


file_list <- list.files(
  path = "data",
  pattern = "*.mat", full.names = F, recursive = F
)

SyncDesyncAnalysis <- function(file) {
#browser()
file_to_load <- file
filename <- as.character(substring(file_to_load, 1, nchar(file_to_load) - 4))
raw.rec <- readMat(file.path("data", file_to_load))

### Action potentials ----------------------------------------------------------------

### contents of ap channel
### title: channel name
### comment: notes about the channel
### resolution: 2.5e-06
### interval: time between two data points in sec (5e-05)
### scale: scaling factor (y scale), real = (data * scale) + offset
### offset: scale offset
### units: unit of the channel (volt)
### length: number of APs
### items: No. points easch AP is represented by
### trigger:
### traces:
### times: times of the APs (first data point of the waveform)
### codes: marker codes of the APs
### values: matrix containing the waveforms of the APs (described by a given No. points)

### number of data points to the peak of the AP from its first data point
points_to_peak <- which(raw.rec$ap[, , 1]$values[150, ] ==
  max(raw.rec$ap[, , 1]$values[150, ])) %>%
  as.numeric()

### time of the peak of the APs after its first point
raw.rec$ap[, , 1]$interval * points_to_peak

ap <- raw.rec$ap[, , 1]$times %>% as.double()
ap_peaks <- ap + c(raw.rec$ap[, , 1]$interval * points_to_peak)

# ap_peaks_tibble <- tibble(peak_times = (ap + c(raw.rec$ap[, , 1]$interval * points_to_peak)))



### EEG ----------------------------------------------------------

### contents of EEG channel
### title: channel name
### comment:
### interval:
### scale: scaling factor (y scale), real = (data * scale) + offset
### offset: scale offset
### units:
### start:
### length:
### values:
EEG <- raw.rec$EEG[, , 1]$values

# samp_rate_original <- 1000/(raw.rec$EEG[,,1]$interval*1000) %>% as.double()
# rec_length_original <- (raw.rec$EEG[,,1]$length/samp_rate_original) %>% as.double()
# times_original <- seq(from = 0,
#                       to = rec_length_original - raw.rec$EEG[, , 1]$interval %>% c(),
#                       by = raw.rec$EEG[, , 1]$interval %>% c())
#
# rec_original <- tibble(
#   time=times_original,
#   EEG=EEG_scaled[,1],
#   AP=NA)
# rec_original$AP[c(lapply(ap_peaks$peak_times,
#                          function(x) {which(abs(times_original-x) == min(abs(times_original-x)))}) %>%
#                     unlist())] = 1
# rec_original$AP %>% sum(na.rm = T)



EEG_scaled <- (EEG * as.double(raw.rec$EEG[, , 1]$scale)) + as.double(raw.rec$EEG[, , 1]$offset)
samp_rate_original <- 1000 / (raw.rec$EEG[, , 1]$interval * 1000) %>% as.double()
rec_length_original <- (raw.rec$EEG[, , 1]$length / samp_rate_original) %>% as.double()



### Filtering ---------------
### function parameters:
### n: order
### W: freq (digital filters: between 0 and 1, 1 is the Nyquist freq. eg: 20/10000 from 20 Hz)


### constructing band pass filter
ButterBP <- signal::butter(
  n = 2,
  W = c(1 / (samp_rate_original / 2), 10 / (samp_rate_original / 2)),
  type = "pass",
  plane = "z"
)

### constructing high pass filter
# ButterHP <- butter(
#   n = 3,
#   W = 200 / (samp_rate / 2),
#   type = "high",
#   plane = "z"
# )

EEG_scaled_bp <- signal::filter(ButterBP, EEG_scaled) %>%
  as.double() %>%
  scale()


### high pass filtering EEG
# EEG_hp <- signal::filter(ButterHP, EEG_scaled) %>% as.double() %>% scale()
# ds_fun_out_hp <- downSamp(data = EEG_hp, ds_factor = 8, samp_rate = 20000)
# plot(ds_fun_out_hp$data[(6*2500):(11*2500)], type = "l")



### downsampling and standardizing in R -----------------------
source(file.path("supplementary_functions", "downSamp.R"))

ds_fun_out <- downSamp(
  data = EEG_scaled_bp,
  ds_factor = 512,
  samp_rate = 20000
)

EEG_ds_scaled <- ds_fun_out$data %>% scale()
samp_rate_ds <- ds_fun_out$data_points / ds_fun_out$r_length
rec_length_ds <- ds_fun_out$r_length
interval_ds <- 1 / samp_rate_ds



# source(file.path("supplementary_functions", "downSamp2.R"))
# ds_fun_out <- downSamp2(
#   data = EEG_scaled,
#   samp_rate_orig = 20000,
#   samp_rate_new = 40
# )
# EEG_ds_scaled <- ds_fun_out$data %>% scale()
# rec_length_ds <- ds_fun_out$r_length_ds
# samp_rate_ds <- ds_fun_out$data_points / ds_fun_out$r_length_ds
# interval_ds <- 1 / samp_rate_ds


### times of data points
times <- seq(from = 0, to = rec_length_ds - interval_ds, by = interval_ds)

### Power spectrum ------
### calculating Power Spectral Density to find the slow oscillation frequency
PSD <- bspec::welchPSD(EEG_ds_scaled, seglength = 256, windowfun = hannwindow)


### finding peaks on PSD
### returns a list with indices [[1]] and "y" values [[2]] of peaks and troughs ([[3]] [[4]])
source(file.path("supplementary_functions", "findpeaks.R"))
peaks <- FindPeaks(PSD$power, bw = 2)
plot(PSD$power, type = "h")
points(peaks[[1]], peaks[[2]], col = "red")

PSD$frequency[peaks[[1]][1]]

### frequency of the first and second peak
first_peak <- (PSD$frequency * samp_rate_ds)[peaks[[1]][1]]
second_peak <- (PSD$frequency * samp_rate_ds)[peaks[[1]][2]]

plot(
  PSD$frequency * samp_rate_ds,
  PSD$power,
  type = "s",
  # lwd = 2,
  xlim = c(0, 20),
  xlab = "Time",
  ylab = "Power"
)
abline(v = first_peak, col = "red")
text(x = first_peak + first_peak / 2, y = max(PSD$power) - 10, round(first_peak, 3), col = "red")

# dev.off()


### calculating moving SDs and means
source(file.path("supplementary_functions", "slideFunct.R"))

### window size: first_peak * samp_rate_ds

slide_sd <- slideFunct(
  data = EEG_ds_scaled,
  window = round(1 / first_peak * samp_rate_ds),
  step = round(1 / first_peak / 2 * samp_rate_ds),
  type = "sd"
) %>%
  scale() %>%
  c() %>%
  rep(each = (length(EEG_ds_scaled) / length(.)) %>% round())
length(slide_sd) <- length(times)

slide_mean <- slideFunct(
  data = EEG_ds_scaled,
  window = round(1 / first_peak * samp_rate_ds),
  step = round(1 / first_peak / 2 * samp_rate_ds),
  type = "mean"
) %>%
  rep(each = (length(EEG_ds_scaled) / length(.)) %>% round())
length(slide_mean) <- length(times)

# slide_cumsum <- slideFunct(
#   data = EEG_ds_scaled,
#   window = round(2 * samp_rate_ds),
#   step = round(2/2 * samp_rate_ds),
#   type = "cumsum"
# ) %>%
#   rep(each = (length(EEG_ds_scaled) / length(.)) %>% round())
# length(slide_mean) <- length(times)




### finding threshold for sync/desync periods (sd_threshold)
par(mfrow = c(2, 1))

EEG_ds_scaled %>%
  ts(
    start = 0,
    end = rec_length_ds - interval_ds,
    deltat = 1 / samp_rate_ds
  ) %>%
  plot(type = "l", col = "gray")
slide_sd %>%
  ts(
    start = 0,
    end = rec_length_ds - interval_ds,
    deltat = 1 / samp_rate_ds
  ) %>%
  points(type = "l")
slide_mean %>%
  ts(
    start = 0,
    end = rec_length_ds - interval_ds,
    deltat = 1 / samp_rate_ds
  ) %>%
  points(type = "l", col = "blue")
abline(h = slide_sd %>% quantile(na.rm = T), col = "red")
text(x = 0, y = c(slide_sd %>% quantile(na.rm = T)), labels = slide_sd %>% quantile(na.rm = T) %>% names(), col = "red")
abline(
  h = (slide_sd %>% median(na.rm = T) + slide_sd %>% quantile(na.rm = T) %>% `[`(2)) / 2,
  lty = 2,
  col = "magenta",
  lwd = 2
)

slide_sd %>% hist(breaks = "FD")
abline(v = slide_sd %>% quantile(na.rm = T), col = "red")
text(x = c(slide_sd %>% quantile(na.rm = T)), y = 50, labels = slide_sd %>% quantile(na.rm = T) %>% names(), col = "red")

par(mfrow = c(1, 1))

### constructing data frame with eeg values, moving SDs and means

sd_threshold <- ((slide_sd %>% median(na.rm = T) + slide_sd %>% quantile(na.rm = T) %>% `[`(2)) / 2) %>% as.numeric()



replace_index <- sapply(
  ap_peaks,
  function(x) {
    ((times - x) %>% abs() %>% min() == (times - x) %>% abs()) %>% which()
  }
)
replace_index %>%
  duplicated() %>%
  plyr::count()

repeat{
  duplicate_check <- any(replace_index %>% duplicated())
  if (duplicate_check == T) {
    replace_index[duplicated(replace_index)] <- replace_index[duplicated(replace_index)] + 1
  } else {
    break()
  }
}

#### burst threshold --------
### burst thershold based on ISI peaks in xlim = c(0,1 / first_peak) window, FILTERS DATA INSIDE
### built in dip test to test ISI uni/multimodality (all ISIs are used)
source(file.path("supplementary_functions", "BurstThresholdDetect.R"))

burst_threshold_detect <- BurstThresholdDetect(
  hist_data = diff(ap_peaks),
  histbreaks = "FD", ### Freedman-Diaconis rule for optimal bin-width
  isi_range = c(0, 1 / first_peak)
)

burst_threshold_detect %>% unlist()
burst_threshold_isi <- burst_threshold_detect[[1]]
length(burst_threshold_isi) = 1

EEG_ds_df <- tibble(
  ### only works with wavelet if it is a matrix!
  times = times,
  eeg_values = as.matrix(EEG_ds_scaled) %>%
    ts(
      start = 0,
      end = rec_length_ds - interval_ds,
      deltat = 1 / samp_rate_ds
    ),
  first_peak,
  sd_threshold,
  burst_threshold_isi,
  nbins = burst_threshold_detect$bins_n,
  moving_sd = slide_sd,
  moving_average = slide_mean,
  levels = ifelse(moving_sd < sd_threshold, yes = 0, no = 1),
  eeg_state = ifelse(levels == 1, yes = "sync", no = "desync"),
  AP = NA,
  # ap_peak_times = ap_peak_times,
  ID = file_to_load
) %>%
  dplyr::mutate(AP = replace(
    x = AP,
    list = c(replace_index),
    values = 1
  ))


EEG_ds_df <- left_join(
  EEG_ds_df,
  EEG_ds_df %>% dplyr::filter(AP == 1) %>%
    dplyr::mutate(
      ap_peak_times = ap_peaks,
      isi = c(diff(ap_peaks), NA),
      burst_isi = ifelse(test = lag(isi > burst_threshold_isi), yes = 1, no = 0)
    )
) %>%
  dplyr::mutate(lag_levels = abs(levels - dplyr::lag(levels, default = 0)) %>% cumsum()) %>%
  dplyr::group_by(lag_levels) %>%
  dplyr::mutate(state_length = times[length(times)] - times[1]) %>%
  dplyr::mutate(state_start = times[1], state_end = times[length(times)]) %>%
  ungroup() %>%
  select(-lag_levels) %>%
  dplyr::mutate(levels2 = ifelse(levels == 0 & state_length <= first_peak * 2, yes = abs(levels - 1), no = levels)) %>%
  dplyr::mutate(levels2 = ifelse(levels == 1 & state_length <= first_peak * 2, yes = abs(levels - 1), no = levels))



return(EEG_ds_df)

}

all_eeg_df <- lapply(file_list, SyncDesyncAnalysis)

EEG_ds_df <- purrr::reduce(all_eeg_df, bind_rows)

# all_eeg_df[[1]]
# EEG_ds_df


EEG_ds_df %>% dplyr::group_by(ID) %>% dplyr::select(sd_threshold) %>% unique()
EEG_ds_df %>% dplyr::group_by(ID) %>% dplyr::select(nbins) %>% unique()
EEG_ds_df %>% dplyr::group_by(ID) %>% dplyr::select(burst_threshold_isi) %>% unique()

gg_isi_periods <- ggplot(
  data = EEG_ds_df %>%
    dplyr::filter(
      AP == 1,
      !is.na(eeg_state),
      state_length >= 1 / first_peak * 2,
      isi > 0, isi < 1 / first_peak
    ), ### NAs in eeg_state because moving_sd has NAs at the end
  mapping = aes(x = isi)
) +
  geom_histogram(aes(
    fill = eeg_state,
    color = eeg_state
  ),
  alpha = .3,
  #bins = burst_threshold_detect$bins_n
  ) +
  # geom_vline(
  #   xintercept = burst_threshold_isi) +
  # annotate(
  #   geom = "text",
  #   x = burst_threshold_isi - burst_threshold_isi / 15,
  #   y = Inf,
  #   label = paste0("Burst threshold = ", burst_threshold_isi),
  #   angle = 90,
  #   hjust = 1.5
  # ) +
  facet_wrap(~ID)
# scale_y_log10() +
# geom_histogram(aes(y = stat(count / sum(count)),
#                    fill = eeg_state,
#                    color = eeg_state),
#                alpha = .6,
#                bins = burst_threshold_detect$bins_n) +
# geom_histogram(aes(y=..count../sum(..count..))) +
# geom_histogram(aes(y=..density.., fill = eeg_state),
#                alpha=0.6,
#                bins = burst_threshold_detect$bins_n,
#                position="identity") +



gg_isi_all <- ggplot(
  data = ap_peaks %>% diff() %>% tibble(isi = .),
  mapping = aes(x = isi)
) +
  geom_histogram()

ggpubr::ggarrange(gg_isi_all, gg_isi_periods, ncol = 2)


left_join(
EEG_ds_df %>%
  dplyr::filter(!is.na(eeg_state), state_length >= 1 / first_peak * 2) %>%
  dplyr::group_by(ID,eeg_state) %>%
  dplyr::summarise(No_APs = length(ap_peak_times)),

EEG_ds_df %>%
  dplyr::filter(!is.na(eeg_state), state_length >= 1 / first_peak * 2) %>%
  dplyr::group_by(ID,eeg_state) %>%
  dplyr::summarise(state_length_sum = sum(state_length %>% unique(), na.rm = T))
) %>% 
  ggplot(aes(x = eeg_state, y = No_APs/state_length_sum, group = ID)) +
  geom_point() +
  geom_line() +
  geom_signif(
    comparisons = list(c("desync", "sync")),
    map_signif_level = TRUE,
    test = "wilcox.test",
    margin_top = .1,
    color = "darkgray",
    
  ) +
  geom_text_repel(
    data = . %>% dplyr::filter(eeg_state == "sync"),
    aes(label = ID),
    direction = "y",
    hjust = -.5)

left_join(
EEG_ds_df %>%
  dplyr::filter(!is.na(eeg_state), state_length >= 1 / first_peak * 2) %>%
  dplyr::group_by(ID,eeg_state) %>%
  dplyr::summarise(No_clusters = sum(burst_isi, na.rm = T)),
EEG_ds_df %>%
  dplyr::filter(!is.na(eeg_state), state_length >= 1 / first_peak * 2) %>%
  dplyr::group_by(ID,eeg_state) %>%
  dplyr::summarise(state_length_sum = sum(state_length %>% unique(), na.rm = T))
) %>% 
  ggplot(aes(x = eeg_state, y = No_clusters/state_length_sum, group = ID)) +
  geom_point(aes(col = ID),size = 3) +
  geom_line(aes(col = ID)) +
  geom_boxplot(
    aes(x = eeg_state, y = No_clusters/state_length_sum), 
    inherit.aes = F,
    alpha = 0, width = 0.2, lwd = 1, fatten = 1.2) +
  geom_signif(
    comparisons = list(c("desync", "sync")),
    map_signif_level = TRUE,
    test = "wilcox.test",
    margin_top = .1,
    color = "darkgray",
    
  ) 
  # geom_text_repel(
  #   data = . %>% dplyr::filter(eeg_state == "sync"),
  #   aes(label = ID),
  #   direction = "y",
  #   hjust = -.5)


EEG_ds_df %>%
  dplyr::filter(!is.na(eeg_state), state_length >= 1 / first_peak * 2) %>%
  dplyr::group_by(ID,eeg_state) %>%
  dplyr::summarise(state_length_sum = sum(state_length %>% unique(), na.rm = T))



### PLOT: eeg, APs (last 80 seconds of the recording) -----------
plotting_region <- c(times[length(times)] - times[length(times)], times[length(times)])
plotting_region <- c(90, 150)


# plotting_region <- c(all_eeg_df[[4]]$times[length(all_eeg_df[[4]]$times)] - 80, all_eeg_df[[4]]$times[length(all_eeg_df[[4]]$times)])



gg_color_hue <- function(n) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

### rect_index: stores the coordinates of the rectangles marking the omitted periods
rect_index <- tibble(
  state_start = subset(EEG_ds_df, state_length <= 1 / first_peak * 2) %>%
    dplyr::select(state_start, state_end) %>% pull(state_start) %>% unique(),
  state_end = subset(EEG_ds_df, state_length <= 1 / first_peak * 2) %>%
    dplyr::select(state_start, state_end) %>% pull(state_end) %>% unique()
)


EEG_ds_df %>% 
  dplyr::filter(times > plotting_region[1], times < plotting_region[2]) %>% 
  dplyr::mutate(lag_levels = abs(levels - dplyr::lag(levels, default = 0)) %>% cumsum()) %>% 
  dplyr::group_by(lag_levels)




eeg_sum_plot <- ggplot(
  data = dplyr::filter(EEG_ds_df, times > plotting_region[1], times < plotting_region[2]),
  mapping = aes(x = times, y = eeg_values)
) +
  geom_line(aes(color = "gray")) +
  geom_line(
    mapping = aes(x = times, y = moving_sd, color = gg_color_hue(2)[1]), size = 1
  ) +
  geom_line(
    mapping = aes(x = times, y = levels + 3, color = gg_color_hue(2)[2]), size = 1
  ) +
  # geom_segment(
  #   mapping = aes(color = "black"),
  #   x = plotting_region[1] - 2, xend = plotting_region[2] + 2,
  #   y = sd_threshold, yend = sd_threshold,
  #   linetype = "dashed"
  # ) +
  geom_point(
    data = . %>% dplyr::filter(!is.na(ap_peak_times)),
    mapping = aes(x = ap_peak_times, y = 5),
    shape = "|",
    size = 4,
    color = "dark gray"
  ) +
  # geom_point(
  #   data = EEG_ds_df,
  #   # %>%
  #   #   dplyr::filter(
  #   #     peak_times > plotting_region[1],
  #   #     peak_times < plotting_region[2]
  #   #   ),
  #   mapping = aes(x = times, y = AP + 7)
  # ) +
  geom_point(
    data = . %>% dplyr::filter(burst_isi == 1),
    mapping = aes(x = times, y = burst_isi + 4),
    color = "black",
    shape = "|",
    size = 5
  ) +
  # geom_rect(
  #   data = dplyr::filter(rect_index, state_start > plotting_region[1], state_start < plotting_region[2]),
  #   inherit.aes = F,
  #   mapping = aes(xmin = state_start, xmax = state_end, ymin = -Inf, ymax = +Inf), fill = "pink", alpha = 0.5
  # ) +
  ylab("EEG [au]") +
  xlab("Time") +

  ### creating legend
  # scale_color_identity(
  #   name = "Lines",
  #   guide = guide_legend(override.aes = list(shape = c(rep(NA, 4)))),
  #   breaks = c("gray", gg_color_hue(2)[1], gg_color_hue(2)[2], "black"),
  #   labels = c("EEG", "Moving SD", "Levels", "SD threshold")
  # ) +
  theme_minimal() +
  facet_wrap(~ID, ncol = 3)

eeg_sum_plot





### Analyzing firing based on EEG periods --------------------------------------------------



### PLOT: ISI ------------------------
### calculating bw based on Freedman-Diaconis rule (max-min)/h (h: bin-width)
sync_isi_data <- ap_peaks_tibble %>%
  dplyr::filter(eeg_state == "sync") %>%
  pull(isi)
desync_isi_data <- ap_peaks_tibble %>%
  dplyr::filter(eeg_state == "desync") %>%
  pull(isi)

bins_isi_sync <- diff(range(sync_isi_data, na.rm = T)) /
  (2 * IQR(sync_isi_data, na.rm = T) /
    length(sync_isi_data)^(1 / 3))

bins_isi_desync <- diff(range(desync_isi_data, na.rm = T)) /
  (2 * IQR(desync_isi_data, na.rm = T) /
    length(desync_isi_data)^(1 / 3))

y_norm <- F
ggplot() +
  geom_histogram(
    data = ap_peaks_tibble %>%
      dplyr::filter(eeg_state == "sync"),
    mapping = aes(
      x = isi,
      y = if (y_norm == F) {
        ..count..
      } else {
        ..count.. / EEG_STATE_length %>%
          dplyr::filter(eeg_state == "sync") %>%
          pull(state_length)
      },
      fill = eeg_state
    ),
    bins = bins_isi_sync
  ) +
  geom_histogram(
    data = ap_peaks_tibble %>%
      dplyr::filter(eeg_state == "desync"),
    mapping = aes(
      x = isi,
      y = if (y_norm == F) {
        ..count..
      } else {
        ..count.. / EEG_STATE_length %>%
          dplyr::filter(eeg_state == "desync") %>%
          pull(state_length)
      },
      fill = eeg_state
    ),
    bins = bins_isi_sync, ### same number of bins on the two ISIs!!
    alpha = 0.5
  ) +
  ylab("Count") +
  xlab("Time (s)")

# ggsave(file.path("output_data", paste0(filename, c("_ISI.png"))),
#        width = 24,
#        height = 18,
#        units = "cm",
#        dpi = 300
# )

# jpeg(file.path("output_data", paste0(filename, c("_isi_with_thresholds.jpg"))),
#      width = 800,
#      height = 600
# )

hist(ap_peaks_tibble %>% dplyr::filter(eeg_state == "sync") %>% pull(isi),
  breaks = "FD"
)

### burst threshold calculated from ISI and its legend
abline(v = burst_threshold_isi[1])
text(
  x = par("usr")[2] - par("usr")[2] / 4,
  y = par("usr")[4] - par("usr")[4] / 10,
  paste0("Burst threshold ISI: ", burst_threshold_isi[1], " s"),
  pos = 2
)

### burst threshold calculated from PSD and its legend
abline(v = burst_threshold_power, col = "red")
text(
  x = par("usr")[2] - par("usr")[2] / 4,
  y = (par("usr")[4] - par("usr")[4] / 10) - par("usr")[4] / 10 / 2,
  paste0("Burst threshold PSD: ", round(burst_threshold_power, 3), " s"),
  col = "red",
  pos = 2
)

### ISI burst threshold calculation window: xlim = c(0,1)
rect(0, 0, 1, par("usr")[4], lty = 2)
text(
  x = par("usr")[2] - par("usr")[2] / 4,
  y = (par("usr")[4] - par("usr")[4] / 10) - 2 * (par("usr")[4] / 10 / 2),
  "---- burst detection window 0 - 1 s",
  pos = 2
)


if (burst_threshold_detect$clustered == F) {
  text(par("usr")[2] - par("usr")[2] / 4,
    y = (par("usr")[4] - par("usr")[4] / 10) - 3 * (par("usr")[4] / 10 / 2),
    pos = 2,
    "Hist. dip test: NOT CLUSTERED",
    col = "red"
  )
} else {
  text(par("usr")[2] - par("usr")[2] / 4,
    y = (par("usr")[4] - par("usr")[4] / 10) - 3 * (par("usr")[4] / 10 / 2),
    pos = 2,
    "Hist. dip test: CLUSTERED",
    col = "green"
  )
}

if (burst_threshold_detect$clustered_d_log == F) {
  text(par("usr")[2] - par("usr")[2] / 4,
    y = (par("usr")[4] - par("usr")[4] / 10) - 4 * (par("usr")[4] / 10 / 2),
    pos = 2,
    "log density dip test: NOT CLUSTERED",
    col = "red"
  )
} else {
  text(par("usr")[2] - par("usr")[2] / 4,
    y = (par("usr")[4] - par("usr")[4] / 10) - 4 * (par("usr")[4] / 10 / 2),
    pos = 2,
    "log density dip test: CLUSTERED",
    col = "green"
  )
}
dev.off()



### AUTOCORRELOGRAM ------------------------------------------------------------

### rel_time matrix creator
RTM_creator <- function(x) {
  mat_spike <- matrix(x,
    nrow = length(x),
    ncol = length(x), byrow = T
  )
  mat_reltimes <- mat_spike - t(mat_spike)
  # mat_reltimes <- data.frame(mat_reltimes)
  mat_reltimes[mat_reltimes == 0] <- NA
  return(mat_reltimes)
}

### calculating rel_time matrices

sync_AC <- RTM_creator(
  ap_peaks_tibble %>%
    dplyr::filter(eeg_state == "sync") %>%
    pull(peak_times)
) %>%
  melt() %>%
  dplyr::mutate(eeg_state = "sync")

desync_AC <- RTM_creator(
  ap_peaks_tibble %>%
    dplyr::filter(eeg_state == "desync") %>%
    pull(peak_times)
) %>%
  melt() %>%
  dplyr::mutate(eeg_state = "desync")

AC <- rbind(sync_AC, desync_AC) %>% ### DOES NOT WORK WITH bind_rows
  transmute(time = value, eeg_state = eeg_state)


### plotting AC

# sync_AC_data <- AC %>%
#   dplyr::filter(eeg_state == "sync") %>%
#   dplyr::filter(time > -1, time < 1) %>%
#   pull(time)
#
# desync_AC_data <- AC %>%
#   dplyr::filter(eeg_state == "desync") %>%
#   dplyr::filter(time > -1, time < 1) %>%
#   pull(time)
#
#
# bins_AC_sync <- diff(range(sync_AC_data, na.rm = T)) /
#   (2 * IQR(sync_AC_data, na.rm = T) /
#      length(sync_AC_data)^(1/3))
#
#
# bins_AC_desync <- diff(range(desync_AC_data, na.rm = T)) /
#   (2 * IQR(desync_AC_data, na.rm = T) /
#      length(desync_AC_data)^(1/3))

### PLOT: ac ----------------
ggplot() +
  geom_histogram(
    data = AC %>%
      dplyr::filter(eeg_state == "sync") %>%
      dplyr::filter(time > -1, time < 1),
    mapping = aes(x = time, fill = eeg_state), bins = 500
  ) +
  geom_histogram(
    data = AC %>%
      dplyr::filter(eeg_state == "desync") %>%
      dplyr::filter(time > -1, time < 1),
    aes(x = time, fill = eeg_state), bins = 500
  ) +
  scale_fill_brewer(type = "div", palette = "Paired")

# ggsave(file.path("output_data", paste0(filename, c("_AC.png"))),
#        width = 24,
#        height = 18,
#        units = "cm",
#        dpi = 300
# )


### WAVELET -------------------------------------------------------------------

time_window <- c(90, 150)

### calculating wavelet
wave <- analyze.wavelet(EEG_ds_df %>%
  dplyr::filter(times > time_window[1], times < time_window[2]),
"eeg_values",
loess.span = 0,
dt = 1 / samp_rate_ds, ### 1/sampling rate (number of intervals/time unit)
dj = 1 / samp_rate_ds / 2,
# lowerPeriod = ,
# upperPeriod = 2,
make.pval = T,
method = "ARIMA",
n.sim = 1
)

### calculating frequency from period
# wave$Freq = c(seq(1,1, length.out = length(wave$Period)))/wave$Period

powers <- wave$Power
freqs <- (1 / wave$Period) # %>% round(2)
period <- (wave$Period) %>% round(2)


powers <- powers %>%
  melt() %>%
  group_by(Var2) %>%
  dplyr::mutate(freqs = freqs) %>%
  ungroup()

### distance between freqs: logarithmic!!! in ggplot it only works if I set y scale to log10
### y axis can be the index of each value (in that case I have to reverse y scale and add the freq values manually)
powers$freqs %>%
  unique() %>%
  diff() %>%
  plot()

powers %>% View()
### PLOT: wavelet --------------------------------------------------------------

ggplot() +
  ### wavelet
  geom_raster(
    data = powers,
    mapping = aes(
      x = Var2 / samp_rate_ds + time_window[1],
      y = freqs, ### Var1: index of values, normal y scale but reversed!!!
      fill = scale(value)
    )
  ) +

  ### EEG
  geom_line(
    data = EEG_ds_df %>%
      dplyr::filter(times > time_window[1], times < time_window[2]),
    mapping = aes(x = times, y = -.5 * eeg_values + 5), color = "white"
  ) +

  ### APs
  geom_point(
    data = EEG_ds_df %>%
      dplyr::filter(times > time_window[1], times < time_window[2]),
    mapping = aes(x = ap_peak_times, y = 10),
    shape = "|",
    color = "white", size = 5
  ) +

  ### first spikes of clusters defined by ISI
  geom_point(
    data = EEG_ds_df %>%
      dplyr::filter(ap_peak_times > time_window[1], ap_peak_times < time_window[2]) %>%
      dplyr::filter(burst_isi == 1),
    mapping = aes(x = ap_peak_times, y = 10),
    shape = "|",
    color = "red", size = 7
  ) +

  ### eeg state (sync: high, desync: low)
  geom_line(
    data = EEG_ds_df %>%
      dplyr::filter(times > time_window[1], times < time_window[2]),
    mapping = aes(x = times, y = .5 * levels + 7), color = "white"
  ) +

  ### PSD first peak
  geom_segment(
    aes(
      x = time_window[1], xend = time_window[2],
      y = first_peak, yend = first_peak
    ),
    col = "gray",
    linetype = "dashed"
  ) +
  # geom_hline(yintercept = first_peak, size = 1, linetype = "dashed", col = "gray") +

  ### theme and scale settings
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(margin = margin(t = 0, r = -30, b = 0, l = 20)),
    axis.text.x = element_text(margin = margin(t = -20, r = 0, b = 20, l = 0)),
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.ticks.length = unit(5, "mm")
  ) +
  xlab("Time (s)") +
  scale_y_log10() +
  # scale_y_reverse(
  #   name = "Frequency (Hz)",
  #   breaks = c(1, 200, 400, 600),
  #   labels = c(
  #     freqs[1],
  #     freqs[200],
  #     freqs[400],
  #     freqs[600]
  #   )
  # ) +
  scale_fill_distiller(palette = "RdGy", name = "Power")


# ggsave(file.path("output_data", paste0(filename, c("_wavelet.png"))),
#        width = 24,
#        height = 18,
#        units = "cm",
#        dpi = 300
# )


### write data to csv
file_exist_test <- file.exists(file.path("output_data", "sync_desync_data.csv"))

write_csv(as_tibble(left_join(EEG_STATE_ap, EEG_STATE_cluster) %>%
  left_join(EEG_STATE_cluster_PSD) %>%
  left_join(EEG_STATE_length) %>%
  dplyr::mutate(ID = file_to_load) %>%
  dplyr::mutate(sd_threshold)) %>%
  dplyr::mutate(burst_threshold_ISI = burst_threshold_isi[1]) %>%
  dplyr::mutate(burst_threshold_PSD = burst_threshold_power) %>%
  dplyr::mutate(clustered = burst_threshold_detect$clustered) %>%
  dplyr::mutate(clustered_d_log = burst_threshold_detect$clustered_d_log) %>%
  dplyr::mutate(MFR = No_APs / state_length),
file.path("output_data", "sync_desync_data.csv"),
append = T,
col_names = !file_exist_test
)

out <- list(
  ID = file_to_load,
  AP_waveforms = raw.rec$ap[, , 1]$values,
  ap_peaks_tibble = ap_peaks_tibble,
  points_to_peak = points_to_peak,
  EEG_ds_df = EEG_ds_df,
  EGG_plot = eeg_sum_plot
)
return(out)
# }


### run on multiple files
tmp <- lapply(file_list, SyncDesyncAnalysis, sd_threshold = -1)
tmp[[1]]$EGG_plot + ggtitle(tmp[[1]]$ID)
