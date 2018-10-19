library(R.matlab)
library(tidyverse)
library(reshape2)
load(file.path("f:","_R_WD","useful_to_load","colorMatrix.RData"))

raw.rec <- readMat(file.path("data","08_right_PnO_t5_4083_baseline1_wake.mat"))

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

points_to_peak <- which(raw.rec$ap[,,1]$values[150,] == max(raw.rec$ap[,,1]$values[150,])) %>% 
  as.numeric()
raw.rec$ap[,,1]$interval*points_to_peak ### time of the peak of the APs after the first point 

ap <- raw.rec$ap[,,1]$times %>% as.double()
ap_peaks <- tibble(peak_times = (ap + c(raw.rec$ap[,,1]$interval*points_to_peak)))


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
EEG <- raw.rec$EEG[,,1]$values

samp_rate <- 1000/(raw.rec$EEG[,,1]$interval*1000) %>% as.double()
rec_length <- (raw.rec$EEG[,,1]$length/samp_rate) %>% as.double()

EEG_scaled <- (EEG*as.double(raw.rec$EEG[,,1]$scale)) + as.double(raw.rec$EEG[,,1]$offset)

### downsampling and standardizing in R
source(file.path("downSamp.R"))
EEG_ds_scaled <- downSamp(data = EEG_scaled, ds_factor = 512, samp_rate = 20000) %>% scale()
samp_rate_ds <- 11676/298.9031
rec_length <- 298.9031
interval_ds <- 1/samp_rate_ds

### Filtering 

### n:order, W:freq (digital filters: between 0 and 1, 1 is the Nyquist freq. eg: 20/10000 from 20 Hz)
ButterHP = butter(n = 3, 
                  W = c(1/(samp_rate_ds/2), 10/(samp_rate_ds/2) ),
                  type = 'pass',
                  plane = "z") 
EEG_ds_scaled <- signal::filter(ButterHP,EEG_ds_scaled) %>% as.double() %>% scale()



### times of data points
times <- seq(from = 0, to = rec_length-interval_ds, by = interval_ds)


### calculate moving sd
source("slideFunct.R")
slide_sd <- slideFunct(data = EEG_ds_scaled, 
                       window = 1*round(samp_rate_ds), 
                       step = 1*round(samp_rate_ds/2), 
                       type = "sd") %>% scale()


slide_mean <- slideFunct(data = EEG_ds_scaled,
                         window = 1*round(samp_rate_ds),
                         step = 1*round(samp_rate_ds/2), 
                         type = "mean") #%>% scale()


### No. datapoints the slide_sd/slide_mean was calculated from
(length(EEG_ds_scaled) / length(slide_sd)) %>% round()
slide_sd_rep <- rep(slide_sd, each = (length(EEG_ds_scaled) / length(slide_sd)) %>% round())

slide_mean_rep <- rep(slide_mean, each = (length(EEG_ds_scaled) / length(slide_mean)) %>% round())

### length difference of the vectors (due to samp_rate_ds)
length(EEG_ds_scaled) - length(slide_sd_rep)

### filling length difference
slide_sd_rep <- c(slide_sd_rep, rep(slide_sd_rep[length(slide_sd_rep)], 
                                    length(EEG_ds_scaled) - length(slide_sd_rep))
                  )
slide_mean_rep <- c(slide_mean_rep, rep(slide_mean_rep[length(slide_mean_rep)], 
                                     length(EEG_ds_scaled) - length(slide_mean_rep))
)

### constructing data frame with eeg values, moving SDs and means
EEG_ds_df <- as.matrix(EEG_ds_scaled) %>% ### only works with wavelet if it is a matrix!
  ts(start = 0, 
     end = rec_length-interval_ds, 
     deltat = 1/samp_rate_ds) %>% 
  data.frame() %>% 
  rename(eeg_values = Series.1) %>% 
  mutate(times) %>%
  mutate(sd = slide_sd_rep) %>% 
  mutate(mean = slide_mean_rep) %>% 
  mutate(levels = 1) %>% 
  mutate(levels = replace(levels, sd < -1, 0)) %>% 
  mutate(ID = "gv120")

### plot eeg with sync desync periods
ggplot(data = EEG_ds_df, mapping = aes(x = times, y = eeg_values)) +
  geom_line() + 
  xlim(rec_length-20,rec_length) + 
  geom_line(data = EEG_ds_df, mapping = aes(x = times, y = sd), color = "red") +
  geom_line(data = EEG_ds_df, mapping = aes(x = times, y = mean), color = "blue") +
  geom_line(data = EEG_ds_df, mapping = aes(x = times, y = levels+5), color = "purple") +
  geom_point(data = ap_peaks %>% 
               dplyr::filter(peak_times > rec_length-20, peak_times < rec_length), 
             mapping = aes(x = peak_times, y = 7),
             shape = "|",
             color = "black", size = 4)



### Analyzing firing based on EEG periods --------------------------------------------------

### finding starts and ends of sync/desync eeg periods
starts <- which(EEG_ds_df %>% pull(levels) - (EEG_ds_df %>% pull(levels) %>% lag()) != 0)
starts <- c(1,starts)

ends <- which(EEG_ds_df %>% pull(levels) - (EEG_ds_df %>% pull(levels) %>% lag()) != 0) - 1
ends <- c(ends, length(EEG_ds_df$levels)) 

if (EEG_ds_df$levels[1] == 1){
  starts_ends <- cbind(starts, ends) %>% 
    as.tibble() %>% 
    mutate(eeg_period = rep(c("sync","desync"), 
                            length(starts))[1:length(starts)])
} else {
  starts_ends <- cbind(starts, ends) %>% 
    as.tibble() %>% 
    mutate(eeg_period = rep(c("desync","sync"), 
                            length(starts))[1:length(starts)])
}

### constructing eeg_periods tibble to store start/end times of sync/desync periods
eeg_periods <- tibble(
  sync_start = EEG_ds_df$times[starts_ends$starts[seq(1,length(starts_ends$starts),2)]],
  sync_end = EEG_ds_df$times[starts_ends$ends[seq(1,length(starts_ends$starts),2)]],
  desync_start = c(EEG_ds_df$times[starts_ends$starts[seq(2,length(starts_ends$starts),2)]],NA),
  desync_end = c(EEG_ds_df$times[starts_ends$ends[seq(2,length(starts_ends$starts),2)]],NA)
) %>% 
  mutate(
    sync_length  = sync_end - sync_start,
    desync_length = desync_end - desync_start)

# level_lengths <- rle(EEG_ds_df$levels)


### list of ap times during sync and desync eeg periods
# sync_ap <- list()
# for (i in 1:length(eeg_periods$sync_start)){
#   sync_ap[[i]] = ap_peaks$peak_times[ap_peaks$peak_times > eeg_periods$sync_start[i] & 
#                                   ap_peaks$peak_times < eeg_periods$sync_end[i]]
# }
# 
# desync_ap <- list()
# for (i in 1:length(eeg_periods$sync_start)){
#   desync_ap[[i]] = ap_peaks$peak_times[ap_peaks$peak_times > eeg_periods$desync_start[i] & 
#                                          ap_peaks$peak_times < eeg_periods$desync_end[i]]
# }

### calculating and plotting inter spike intervals during sync and desync eeg periods
# desync_isi <- lapply(desync_ap, diff) %>% 
#   unlist() %>% 
#   as.data.frame() %>% 
#   rename(desync_ISI = ".") %>% 
#   melt() %>% 
#   na.omit()
#   
# sync_isi <- lapply(sync_ap, diff) %>% 
#   unlist() %>% 
#   as.data.frame() %>% 
#   rename(sync_ISI = ".") %>% 
#   melt() %>% 
#   na.omit()
# 
# 
# ISI <- rbind(sync_isi,desync_isi) %>% ### DOES NOT WORK WITH bind_rows 
#   rename(time = value, eeg_state = variable) %>% 
#   mutate(ID = "gv25")
# 
# ggplot(data = ISI, mapping = aes(x = time, fill = eeg_state)) +
#   geom_histogram(bins = 150) +
#   scale_fill_brewer(type = "div", palette = "Paired") 
#   scale_y_continuous(trans='log10')


### detecting AP clusters ("bursts") based on ISI (distance between local maxima / 2)
BurstThresholdDetect <- function(hist_data, histbreaks){
  isihist <- hist(hist_data, breaks = histbreaks, plot = F)
  max1 = which(
    isihist$counts[1:(length(isihist$counts)/2)] == 
      isihist$counts[1:(length(isihist$counts)/2)] %>% max()
  ) %>% as.numeric()
  
  max2 = which(
    isihist$counts[(length(isihist$counts)/2+1):length(isihist$counts)] == 
      isihist$counts[(length(isihist$counts)/2+1):length(isihist$counts)] %>% max()
  ) + length(isihist$counts)/2
  
  
  hist(hist_data, breaks = histbreaks)
  abline(v = isihist$mids[max1], col = "red")
  abline(v = isihist$mids[max2], col = "red")
  abline(v = isihist$mids[((max2-max1)/2)+max1], col = "red")
  burst_threshold <- isihist$mids[((max2-max1)/2)+max1]
  return(burst_threshold)
}

burst_threshold <-  BurstThresholdDetect(hist_data = 
                                           diff(ap_peaks$peak_times),
                     histbreaks = 150)

### adding columns to ap_peaks 
ap_peaks <- ap_peaks %>% 
  dplyr::mutate(eeg_state = "desync") %>% 
  dplyr::mutate(isi = c(diff(peak_times), NA)) %>% 
  dplyr::mutate(burst = 0) 
ap_peaks <- ap_peaks %>% 
  dplyr::mutate(burst = replace(burst, isi > burst_threshold[1], 1))

for (i in 1:length(eeg_periods$sync_start)){
  ap_peaks <- ap_peaks %>% 
    dplyr::mutate(eeg_state = replace(
      eeg_state,
      peak_times > eeg_periods$sync_start[i] & peak_times < eeg_periods$sync_end[i],
      values = "sync"))
}    

ap_peaks %>% 
  group_by(eeg_state) %>% 
  summarise(length(isi)) 

ap_peaks %>% 
  group_by(eeg_state) %>% 
  summarise(sum(burst)) 


eeg_periods$desync_length %>% sum(na.rm = T)
eeg_periods$sync_length %>% sum(na.rm = T)

hist(ap_peaks$isi[ap_peaks$eeg_state == "sync"],breaks = 150)
hist(ap_peaks$isi[ap_peaks$eeg_state == "desync"],breaks = 150, col = "red", add = T)
hist(ap_peaks$isi,breaks = 150)

ggplot(data = ap_peaks, aes(isi, fill = eeg_state)) +
  geom_histogram(bins = 150, na.rm = T, position="identity")+
  guides(fill = guide_legend(reverse = TRUE))
 



# isihist$counts %>% sort(decreasing = T) %>% `[[`(2)


### calculating and plotting autocorrelation during sync and desync eeg periods 

### rel_time matrix creator for lapply 
RTM_creator <- function(x) {
  mat_spike <- matrix(x, 
                      nrow = length(x), 
                      ncol = length(x), byrow = T
  )
  mat_reltimes <- mat_spike - t(mat_spike)
  #mat_reltimes <- data.frame(mat_reltimes)
  mat_reltimes[mat_reltimes == 0] = NA 
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

AC <- rbind(sync_AC,desync_AC) %>% ### DOES NOT WORK WITH bind_rows 
  transmute(time = value, eeg_state = eeg_state) 

AC %>% 
  group_by(eeg_state) %>% 
  summarise(length(time))

hist(AC %>% 
       dplyr::filter(eeg_state == "sync") %>% 
       dplyr::filter(time > -1, time <1) %>% 
       pull(time), breaks = 500)

hist(AC %>% 
       dplyr::filter(eeg_state == "desync") %>% 
       dplyr::filter(time > -1, time <1) %>% 
       pull(time), breaks = 500, 
     border  = "red", 
     add = T)


ggplot(data = AC %>% 
         dplyr::filter(time > -1, time <1),
       mapping = aes(x = time, fill = eeg_state)) + 
  geom_histogram(bins = 500, position="identity") +
  scale_fill_brewer(type = "div", palette = "Paired")

# sync_AC <- lapply(sync_ap, RTM_creator) %>% 
#   unlist() %>% 
#   as.data.frame() %>% 
#   rename(sync_AC = ".") %>% 
#   melt()
#   
# 
# desync_AC <- lapply(desync_ap, RTM_creator) %>% 
#   unlist() %>% 
#   as.data.frame() %>% 
#   rename(desync_AC = ".") %>%   
#   melt()
# 
# AC <- rbind(sync_AC,desync_AC) %>% ### DOES NOT WORK WITH bind_rows 
#   rename(time = value) %>% 
#   mutate(ID = "gv03")



### Wavelet -------------------------------------------------------------------

library(WaveletComp)

time_window = c(50,100)
wave <- analyze.wavelet(EEG_ds_df %>% 
                          dplyr::filter(times > time_window[1], times < time_window[2]), 
                        "eeg_values", 
                        loess.span = 0,
                        dt = 1/40,
                        ###1/sampling rate (number of intervals/time unit)  
                        dj = 1/20,
                        #lowerPeriod = ,
                        #upperPeriod = 2,
                        make.pval = T, 
                        method = 'ARIMA',
                        n.sim = 1)

### calculate frequency from period
#wave$Freq = c(seq(1,1, length.out = length(wave$Period)))/wave$Period

powers <- wave$Power
powers %>% dim()
freqs <- (1/wave$Period) %>% round(2)
period <- (wave$Period) %>%  round(2)


library(reshape2)
library(ggplot2)

powers <- powers %>% 
  melt() #%>% 
#mutate(freq = rep(freqs, 801))


wt.image(wave, color.key = "interval", n.levels = 150,
         legend.params = list(lab = "wavelet power levels", mar = 4.7),
         plot.coi = F,
         plot.contour = F,
         plot.ridge = F
)


### Plotting wavelet --------------------------------------------------------------


ggplot() + 
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(margin = margin(t = 0, r = -30, b = 0, l = 20)),
        axis.text.x = element_text(margin = margin(t = -20, r = 0, b = 20, l = 0)),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.ticks.length = unit(5, "mm")) +
  geom_raster(data = powers,
              mapping = aes(x = Var2/samp_rate_ds + time_window[1],
                            y = Var1,
                            fill = scale(value))) +
  xlab("Time (s)") +
  scale_y_reverse(name = "Frequency (Hz)",
                     breaks = c(1,25,50,75,100,125,150),
                     labels = c(freqs[1],
                                freqs[25],
                                freqs[50],
                                freqs[75],
                                freqs[100],
                                freqs[125],
                                freqs[150])) +
  #scale_y_reverse() +
  scale_fill_distiller(palette = "RdGy", name = "Power") +
  geom_line(data = EEG_ds_df %>% 
              dplyr::filter(times > time_window[1], times < time_window[2]), 
            mapping = aes(x = times, y = -2*eeg_values+20), color = "white") + 
  geom_point(data = ap_peaks %>% 
                dplyr::filter(peak_times > time_window[1], peak_times < time_window[2]), 
             mapping = aes(x = peak_times, y = 100),
             shape = "|",
             color = "white", size = 4) +
   geom_line(data = EEG_ds_df %>% 
               dplyr::filter(times > time_window[1], times < time_window[2]), 
             mapping = aes(x = times, y = -3*levels+5), color = "white") 
  

