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

file_list <- list.files(
  path = "data",
  pattern = "*.mat", full.names = F, recursive = F
)
sd_threshold = 1
SyncDesyncAnalysis <- function(file, sd_threshold) {
  file_to_load <- file_list[1]
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
  ap_peaks <- tibble(peak_times = (ap + c(raw.rec$ap[, , 1]$interval * points_to_peak)))
  
  
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
  
  # samp_rate <- 1000/(raw.rec$EEG[,,1]$interval*1000) %>% as.double()
  # rec_length_ds <- (raw.rec$EEG[,,1]$length/samp_rate) %>% as.double()
  
  EEG_scaled <- (EEG * as.double(raw.rec$EEG[, , 1]$scale)) + as.double(raw.rec$EEG[, , 1]$offset)
  
  ### downsampling and standardizing in R -----------------------
  source(file.path("supplementary_functions", "downSamp.R"))
  
  ds_fun_out <- downSamp(data = EEG_scaled, 
                         ds_factor = 512, 
                         samp_rate = 20000)
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
  
  
  ### Filtering ---------------
  ### function parameters:
  ### n: order
  ### W: freq (digital filters: between 0 and 1, 1 is the Nyquist freq. eg: 20/10000 from 20 Hz)
  
  ButterBP <- butter(
    n = 3,
    W = c(1 / (samp_rate_ds / 2), 10 / (samp_rate_ds / 2)),
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
  
  EEG_ds_scaled <- signal::filter(ButterBP, EEG_ds_scaled) %>% as.double() %>% scale()
  
  ### high pass filtering EEG
  # EEG_hp <- signal::filter(ButterHP, EEG_scaled) %>% as.double() %>% scale()
  # ds_fun_out_hp <- downSamp(data = EEG_hp, ds_factor = 8, samp_rate = 20000)
  # plot(ds_fun_out_hp$data[(6*2500):(11*2500)], type = "l")
  
  
  ### times of data points
  times <- seq(from = 0, to = rec_length_ds - interval_ds, by = interval_ds)
  
  
  ### calculating moving SDs and means
  source(file.path("supplementary_functions", "slideFunct.R"))
  
  slide_sd <- slideFunct(
    data = EEG_ds_scaled,
    window = 1 * round(samp_rate_ds),
    step = 1 * round(samp_rate_ds / 2),
    type = "sd"
  ) %>% scale()
  
  slide_mean <- slideFunct(
    data = EEG_ds_scaled,
    window = 1 * round(samp_rate_ds),
    step = 1 * round(samp_rate_ds / 2),
    type = "mean"
  ) # %>% scale()
  
  
  ### No. datapoints the slide_sd/slide_mean was calculated from
  (length(EEG_ds_scaled) / length(slide_sd)) %>% round()
  slide_sd_rep <- rep(slide_sd, each = (length(EEG_ds_scaled) / length(slide_sd)) %>% round())
  
  slide_mean_rep <- rep(slide_mean, each = (length(EEG_ds_scaled) / length(slide_mean)) %>% round())
  
  ### length difference of the vectors
  length_diff <- length(EEG_ds_scaled) - length(slide_sd_rep)
  # ((rec_length_ds - as.integer(rec_length_ds)) * samp_rate_ds) %>% round()
  
  
  ### filling length difference
  slide_sd_rep <- c(slide_sd_rep, rep(
    slide_sd_rep[length(slide_sd_rep)],
    length_diff
  ))
  slide_mean_rep <- c(slide_mean_rep, rep(
    slide_mean_rep[length(slide_mean_rep)],
    length_diff
  ))
  
  ### constructing data frame with eeg values, moving SDs and means
  SD_THRESHOLD <- sd_threshold %>% as.numeric()
  EEG_ds_df <- as.matrix(EEG_ds_scaled) %>% ### only works with wavelet if it is a matrix!
    ts(
      start = 0,
      end = rec_length_ds - interval_ds,
      deltat = 1 / samp_rate_ds
    ) %>%
    data.frame() %>%
    rename(eeg_values = Series.1) %>%
    mutate(times) %>%
    mutate(sd = slide_sd_rep) %>%
    mutate(mean = slide_mean_rep) %>%
    mutate(levels = 1) %>%
    mutate(levels = replace(levels, sd < SD_THRESHOLD, 0)) %>%
    mutate(ID = "gv120")
  
  
  
  ### Find burst_threshold based on power spectrum (oscillation frequency of the EEG)
  
  ### calculating Power Spectral Density
  PSD <- bspec::welchPSD(EEG_ds_df %>% dplyr::filter(levels == 1) %>% pull(eeg_values),
                         seglength = 256, windowfun = hannwindow
  )
  
  ### finding peaks on PSD
  ### returns a list with indices [[1]] and values [[2]] of peaks and troughs [[3]] [[4]]
  source(file.path("supplementary_functions", "findpeaks.R"))
  peaks <- FindPeaks(PSD$power, bw = 3)
  plot(PSD$power, type = "h")
  points(peaks[[1]], peaks[[2]], col = "red")
  
  ### frequency of the second peak
  second_peak <- (PSD$frequency * samp_rate_ds)[peaks[[1]][1]]
  
  # jpeg(file.path("output_data", paste0(filename, c("_burst_threshold_PSD.jpg"))),
  #      width = 800,
  #      height = 600
  # )
  
  plot(
    PSD$frequency * samp_rate_ds,
    PSD$power,
    type = "h",
    lwd = 2,
    xlim = c(0, 20),
    xlab = "Time",
    ylab = "Power"
  )
  abline(v = second_peak, col = "red")
  text(x = second_peak + 1, y = max(PSD$power) - 10, round(second_peak, 3), col = "red")
  
  dev.off()
  
  ### calculating time from frequency in miliseconds and in seconds
  burst_threshold_power <- 1000 / second_peak / 1000
  
  
  ### PLOT: eeg, APs (last 80 seconds of the recording) -----------
  plotting_region <- c(times[length(times)] - 80, times[length(times)])
  
  ggplot(data = EEG_ds_df, mapping = aes(x = times, y = eeg_values)) +
    geom_line() +
    xlim(plotting_region[1], plotting_region[2]) +
    geom_line(data = EEG_ds_df, mapping = aes(x = times, y = sd), color = "red") +
    geom_line(data = EEG_ds_df, mapping = aes(x = times, y = mean), color = "blue") +
    geom_line(data = EEG_ds_df, mapping = aes(x = times, y = levels + 5), color = "purple") +
    geom_point(
      data = ap_peaks %>%
        dplyr::filter(
          peak_times > plotting_region[1],
          peak_times < plotting_region[2]
        ),
      mapping = aes(x = peak_times, y = 7),
      shape = "|",
      color = "black", size = 4
    )
  
  # ggsave(file.path("output_data", paste0(filename, c("_sync_desync_periods.png"))),
  #        width = 24,
  #        height = 18,
  #        units = "cm",
  #        dpi = 300
  # )
  
  
  
  ### Analyzing firing based on EEG periods --------------------------------------------------
  
  ### finding start and end indices of sync/desync eeg periods
  starts <- which(EEG_ds_df %>% pull(levels) - (EEG_ds_df %>% pull(levels) %>% lag()) != 0)
  starts <- c(1, starts)
  
  ends <- which(EEG_ds_df %>% pull(levels) - (EEG_ds_df %>% pull(levels) %>% lag()) != 0) - 1
  ends <- c(ends, length(EEG_ds_df$levels))
  
  if (EEG_ds_df$levels[1] == 1) { ### starting with sync period
    starts_ends <- cbind(starts, ends) %>%
      as.tibble() %>%
      mutate(eeg_period = rep(
        c("sync", "desync"),
        length(starts)
      )[1:length(starts)])
  } else { ### starting with sync period
    starts_ends <- cbind(starts, ends) %>%
      as.tibble() %>%
      mutate(eeg_period = rep(
        c("desync", "sync"),
        length(starts)
      )[1:length(starts)])
  }
  
  starts_ends %>%
    group_by(eeg_period) %>%
    summarise(length(starts))
  
  
  ### constructing eeg_periods tibble to store start/end times of sync/desync periods
  
  if (starts_ends$starts %>% length() %% 2 == 0) { ### odd and even cases
    eeg_periods <- tibble(
      sync_start = EEG_ds_df$times[starts_ends$starts[seq(
        1,
        length(starts_ends$starts),
        2
      )]                                   ],
      sync_end = EEG_ds_df$times[starts_ends$ends[seq(
        1,
        length(starts_ends$starts),
        2
      )]                                 ],
      desync_start = EEG_ds_df$times[starts_ends$starts[seq(
        2,
        length(starts_ends$starts),
        2
      )]                                     ],
      desync_end = EEG_ds_df$times[starts_ends$ends[seq(
        2,
        length(starts_ends$starts),
        2
      )]                                   ]
    ) %>%
      mutate(
        sync_length = sync_end - sync_start,
        desync_length = desync_end - desync_start
      )
  } else {
    eeg_periods <- tibble(
      sync_start = EEG_ds_df$times[starts_ends$starts[seq(
        1,
        length(starts_ends$starts),
        2
      )]                                   ],
      sync_end = EEG_ds_df$times[starts_ends$ends[seq(
        1,
        length(starts_ends$starts),
        2
      )]                                 ],
      desync_start = c(
        EEG_ds_df$times[starts_ends$starts[seq(
          2,
          length(starts_ends$starts),
          2
        )]                                       ],
        NA
      ),
      desync_end = c(
        EEG_ds_df$times[starts_ends$ends[seq(
          2,
          length(starts_ends$starts),
          2
        )]                                     ],
        NA
      )
    ) %>%
      mutate(
        sync_length = sync_end - sync_start,
        desync_length = desync_end - desync_start
      )
  }
  
  
  
  ### adding columns to ap_peaks
  ap_peaks <- ap_peaks %>%
    dplyr::mutate(eeg_state = "desync") %>%
    dplyr::mutate(isi = c(diff(peak_times), NA)) %>%
    dplyr::mutate(burst = 0) %>%
    dplyr::mutate(burst_PSD = 0)
  
  ### finding sync periods based on conditions
  for (i in 1:length(eeg_periods$sync_start)) {
    ap_peaks <- ap_peaks %>%
      dplyr::mutate(eeg_state = replace(
        eeg_state,
        peak_times > eeg_periods$sync_start[i] & peak_times < eeg_periods$sync_end[i],
        values = "sync"
      ))
  }
  
  ### burst thershold based on ISI peaks in xlim = c(0,1) window, FILTERS DATA INSIDE
  ### built in dip test to test ISI uni/multimodality (all ISIs are used)
  source(file.path("supplementary_functions", "BurstThresholdDetect.R"))
  
  burst_threshold_detect <- BurstThresholdDetect(
    hist_data = ap_peaks %>%
      dplyr::filter(eeg_state == "sync" | eeg_state == "desync") %>%
      pull(isi),
    histbreaks = "FD" ### Freedman-Diaconis rule for optimal bin-width
  )
  
  burst_threshold_isi <- burst_threshold_detect[[1]]
  
  
  ap_peaks <- ap_peaks %>%
    dplyr::mutate(burst = replace(burst, lag(isi > burst_threshold_isi[1]), 1)) %>%
    dplyr::mutate(burst_PSD = replace(
      burst_PSD,
      lag(isi > burst_threshold_power),
      1
    )) ### lag(): to find the firs APs in a cluster (burst) instead of the last one
  
  
  ### save to file at the end:
  EEG_STATE_ap <- ap_peaks %>%
    group_by(eeg_state) %>%
    summarise(length(peak_times)) %>%
    rename(No_APs = "length(peak_times)")
  
  EEG_STATE_cluster <- ap_peaks %>%
    group_by(eeg_state) %>%
    summarise(sum(burst)) %>%
    rename(No_clusters = "sum(burst)")
  
  EEG_STATE_cluster_PSD <- ap_peaks %>%
    group_by(eeg_state) %>%
    summarise(sum(burst_PSD)) %>%
    rename(No_clusters_PSD = "sum(burst_PSD)")
  
  EEG_STATE_length <- tibble(
    eeg_state = c("desync", "sync"),
    state_length = c(
      eeg_periods$desync_length %>% sum(na.rm = T),
      eeg_periods$sync_length %>% sum(na.rm = T)
    )
  )
  
  ### PLOT: ISI ------------------------
  ### calculating bw based on Freedman-Diaconis rule (max-min)/h (h: bin-width)
  sync_isi_data <- ap_peaks %>% dplyr::filter(eeg_state == "sync") %>% pull(isi)
  desync_isi_data <- ap_peaks %>% dplyr::filter(eeg_state == "desync") %>% pull(isi)
  
  bins_isi_sync <- diff(range(sync_isi_data, na.rm = T)) /
    (2 * IQR(sync_isi_data, na.rm = T) /
       length(sync_isi_data)^(1 / 3))
  
  bins_isi_desync <- diff(range(desync_isi_data, na.rm = T)) /
    (2 * IQR(desync_isi_data, na.rm = T) /
       length(desync_isi_data)^(1 / 3))
  
  y_norm <- F
  ggplot() +
    geom_histogram(
      data = ap_peaks %>%
        dplyr::filter(eeg_state == "sync"),
      mapping = aes(x = isi, 
                    y = if(y_norm == F) {
                      ..count..
                    } else { 
                      ..count../EEG_STATE_length %>% 
                        dplyr::filter(eeg_state == "sync") %>% 
                        pull(state_length)
                    },
                    fill = eeg_state),
      bins = bins_isi_sync
    ) +
    geom_histogram(
      data = ap_peaks %>%
        dplyr::filter(eeg_state == "desync"),
      mapping = aes(x = isi, 
                    y = if(y_norm == F) {
                      ..count..
                    } else { 
                      ..count../EEG_STATE_length %>% 
                        dplyr::filter(eeg_state == "desync") %>% 
                        pull(state_length)
                    }, 
                    fill = eeg_state),
      bins = bins_isi_sync ### same number of bins on the two ISIs!!
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
  
  hist(ap_peaks %>% dplyr::filter(eeg_state == "sync") %>% pull(isi),
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
    ap_peaks %>%
      dplyr::filter(eeg_state == "sync") %>%
      pull(peak_times)
  ) %>%
    melt() %>%
    mutate(eeg_state = "sync")
  
  desync_AC <- RTM_creator(
    ap_peaks %>%
      dplyr::filter(eeg_state == "desync") %>%
      pull(peak_times)
  ) %>%
    melt() %>%
    mutate(eeg_state = "desync")
  
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
  
  time_window <- c(50, 100)
  
  ### calculating wavelet
  wave <- analyze.wavelet(EEG_ds_df %>%
                            dplyr::filter(times > time_window[1], times < time_window[2]),
                          "eeg_values",
                          loess.span = 0,
                          dt = 1 / 40, ### 1/sampling rate (number of intervals/time unit)
                          dj = 1 / 20,
                          # lowerPeriod = ,
                          # upperPeriod = 2,
                          make.pval = T,
                          method = "ARIMA",
                          n.sim = 1
  )
  
  ### calculating frequency from period
  # wave$Freq = c(seq(1,1, length.out = length(wave$Period)))/wave$Period
  
  powers <- wave$Power
  freqs <- (1 / wave$Period) %>% round(2)
  period <- (wave$Period) %>% round(2)
  
  
  powers <- powers %>%
    melt() # %>%
  
  
  ### PLOT: wavelet --------------------------------------------------------------
  
  ggplot() +
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
    ### wavelet
    geom_raster(
      data = powers,
      mapping = aes(
        x = Var2 / samp_rate_ds + time_window[1],
        y = Var1,
        fill = scale(value)
      )
    ) +
    xlab("Time (s)") +
    scale_y_reverse(
      name = "Frequency (Hz)",
      breaks = c(1, 25, 50, 75, 100, 125, 150),
      labels = c(
        freqs[1],
        freqs[25],
        freqs[50],
        freqs[75],
        freqs[100],
        freqs[125],
        freqs[150]
      )
    ) +
    scale_fill_distiller(palette = "RdGy", name = "Power") +
    
    ### EEG
    geom_line(
      data = EEG_ds_df %>%
        dplyr::filter(times > time_window[1], times < time_window[2]),
      mapping = aes(x = times, y = -2 * eeg_values + 20), color = "white"
    ) +
    
    ### APs
    geom_point(
      data = ap_peaks %>%
        dplyr::filter(peak_times > time_window[1], peak_times < time_window[2]),
      mapping = aes(x = peak_times, y = 100),
      shape = "|",
      color = "white", size = 5
    ) +
    
    ### first spikes of clusters defined by PSD
    # (if (any(ap_peaks$burst_PSD == 1)) {
    #   geom_point(
    #     data = ap_peaks %>%
    #       dplyr::filter(peak_times > time_window[1], peak_times < time_window[2]) %>%
    #       dplyr::filter(burst_PSD == 1),
    #     mapping = aes(x = peak_times, y = 96),
    #     shape = "ˇ",
    #     color = "green", size = 7
    #   )
  # }) +
  
  
  ### first spikes of clusters defined by ISI
  geom_point(
    data = ap_peaks %>%
      dplyr::filter(peak_times > time_window[1], peak_times < time_window[2]) %>%
      dplyr::filter(burst == 1),
    mapping = aes(x = peak_times, y = 100),
    shape = "|",
    color = "red", size = 7
  ) +
    
    ### eeg state (sync: high, desync: low)
    geom_line(
      data = EEG_ds_df %>%
        dplyr::filter(times > time_window[1], times < time_window[2]),
      mapping = aes(x = times, y = -3 * levels + 5), color = "white"
    )
  
  
  # ggsave(file.path("output_data", paste0(filename, c("_wavelet.png"))),
  #        width = 24,
  #        height = 18,
  #        units = "cm",
  #        dpi = 300
  # )
  
  
  ### write data to csv
  file_exist_test <- file.exists(file.path("output_data", "sync_desync_data.csv"))
  
  write_csv(as.tibble(left_join(EEG_STATE_ap, EEG_STATE_cluster) %>%
                        left_join(EEG_STATE_cluster_PSD) %>%
                        left_join(EEG_STATE_length) %>%
                        mutate(ID = file_to_load) %>%
                        mutate(sd_threshold = SD_THRESHOLD)) %>%
              mutate(burst_threshold_ISI = burst_threshold_isi[1]) %>%
              mutate(burst_threshold_PSD = burst_threshold_power) %>%
              mutate(clustered = burst_threshold_detect$clustered) %>% 
              mutate(clustered_d_log = burst_threshold_detect$clustered_d_log) %>% 
              mutate(MFR = No_APs / state_length),
            file.path("output_data", "sync_desync_data.csv"),
            append = T,
            col_names = !file_exist_test
  )
}


### run on multiple files
lapply(file_list, SyncDesyncAnalysis, sd_threshold = -1)


### LOADING DATA (sync_desync tibble) ###############
### Analyzing results
SYNC_DESYNC_RESULT <- read_csv(file.path("output_data", "sync_desync_data.csv"))

SYNC_DESYNC_RESULT <- SYNC_DESYNC_RESULT %>%
  mutate(No_clusters_norm = No_clusters / state_length) %>%
  mutate(No_clusters_PSD_norm = No_clusters_PSD / state_length)



y_axis <- "No_clusters_norm"

### PLOT : Sync vs Desync (MFR/Clusters) ----------------------
ggplot(data = SYNC_DESYNC_RESULT %>% 
         dplyr::filter(clustered_d_log == T 
                       #clustered_d_log == F
         ),
       mapping = aes(
         x = forcats::fct_relevel(eeg_state, "sync", "desync"),
         y = eval(parse(text = y_axis)))
) +
  theme_minimal() +
  theme(axis.text = element_text(size = 20),
        text = element_text(size = 20),
        legend.text = element_text(size = 20)) +
  xlab("LFP period")+
  scale_x_discrete(labels = c('Synchronous','Deynchronous')) +
  ylab("MFR") +
  
  geom_point(
    aes(fill = forcats::fct_relevel(as.factor(clustered_d_log),"TRUE", "FALSE")), 
    shape = 21,
    size = 4
  ) +
  geom_line(
    mapping = aes(group = ID,
                  color = forcats::fct_relevel(as.factor(clustered_d_log),"TRUE", "FALSE")
    ),
    lwd = 1
  ) +
  geom_signif(comparisons = list(c("sync", "desync")),
              map_signif_level = FALSE,
              test = "wilcox.exact",
              test.args = c("paired" = T)
  ) +
  
  scale_fill_discrete(name = "",
                      breaks = c("FALSE", "TRUE"),
                      labels = c("Non clustered", "Clustered")
                      #guide = F
  ) +
  scale_color_discrete(name = "",
                       breaks = c("FALSE", "RTUE"),
                       labels = c("Non clustered", "Clustered") 
                       #guide = F
  )


ggsave(file.path("output_data","MFR_clust_vs_nonclust.png"),
       width = 8,
       height = 12,
       dpi = 300)

