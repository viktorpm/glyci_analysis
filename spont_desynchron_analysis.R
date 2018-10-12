library(R.matlab)
library(tidyverse)


raw.rec <- readMat(file.path("x:","extra32","_Projects","Glyci_II","DATA","ELPHYS","juxta","GII-25","08_right_PnO_t5_4083_baseline1_wake.mat"))


raw.rec$ap[,,1]$interval
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


which(raw.rec$ap[,,1]$values[150,] == max(raw.rec$ap[,,1]$values[150,]))
raw.rec$ap[,,1]$interval*6*1000 ### time of the peak of the APs after the first point 

ap <- raw.rec$ap[,,1]$times %>% as.double()
ap_peaks <- ap + (raw.rec$ap[,,1]$interval*6)
ap_peaks <- ap_peaks %>% as.tibble()

### EEG ----------------------------------------------------------
EEG <- raw.rec$EEG[,,1]$values
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

samp_rate <- 1000/(raw.rec$EEG[,,1]$interval*1000) %>% as.double()
rec_length <- (raw.rec$EEG[,,1]$length/samp_rate) %>% as.double()
EEG_scaled <- (EEG*as.double(raw.rec$EEG[,,1]$scale)) + as.double(raw.rec$EEG[,,1]$offset)

### EEG_ds (downsampled) --------------------------------------------------------
EEG_ds <- raw.rec$eeg.ds[,,1]$values %>% as.double()
raw.rec$eeg.ds[,,1]$length/40

samp_rate_ds <- 1000/(raw.rec$eeg.ds[,,1]$interval*1000) %>% as.double()
rec_length <- (raw.rec$eeg.ds[,,1]$length/samp_rate_ds) %>% as.double()
interval_ds <- raw.rec$eeg.ds[,,1]$interval %>% as.double()

EEG_ds_scaled <- (EEG_ds*as.double(raw.rec$eeg.ds[,,1]$scale)) + 
  as.double(raw.rec$eeg.ds[,,1]$offset)

# -----------------------------------------------------------------

### calculates time window
plot_from <- 250*samp_rate_ds %>% as.double()
plot_to <- 290*samp_rate_ds %>% as.double()
plot(EEG_ds_scaled[plot_from:plot_to], type = "l")

### times of data points
times <- seq(from = 0, to = rec_length-interval_ds, by = interval_ds)

### calculate moving sd
source("slideFunct.R")
slide_sd <- slideFunct(data = EEG_ds_scaled, window = 40, step = 20, type = "sd")
slide_mean <- slideFunct(data = EEG_ds_scaled, window = 40, step = 20, type = "mean")

### No. datapoints the sd/mean was calculated from
(length(EEG_ds_scaled) / length(slide_sd)) %>% round()
slide_sd_rep <- rep(slide_sd, each = (length(EEG_ds_scaled) / length(slide_sd)) %>% round())

slide_mean_rep <- rep(slide_mean, each = (length(EEG_ds_scaled) / length(slide_mean)) %>% round())


### length difference of the vectors
length(EEG_ds_scaled) - length(slide_sd_rep)


### filling length difference
slide_sd_rep <- c(slide_sd_rep, rep(slide_sd_rep[length(slide_sd_rep)], 
                                    length(EEG_ds_scaled) - length(slide_sd_rep))
                  )
slide_mean_rep <- c(slide_mean_rep, rep(slide_mean_rep[length(slide_mean_rep)], 
                                     length(EEG_ds_scaled) - length(slide_mean_rep))
)


EEG_ds_df <- EEG_ds_scaled %>% 
  ts(start = 0, 
     end = rec_length-interval_ds, 
     deltat = 1/samp_rate_ds) %>% 
  as.data.frame() %>% 
  mutate(times) %>%
  mutate(sd = slide_sd_rep) %>% 
  mutate(mean = slide_mean_rep)

ggplot(data = EEG_ds_df, mapping = aes(x = times, y = x)) +
  geom_line() + 
  xlim(100,200) + 
  geom_line(data = EEG_ds_df, mapping = aes(x = times, y = sd), color = "red") +
  geom_line(data = EEG_ds_df, mapping = aes(x = times, y = mean), color = "blue")


  
identical(
EEG_ds_df %>% 
  dplyr::slice(plot_from:plot_to),

EEG_ds_df %>% 
  dplyr::filter(times > 250, times < 290)
)



# small <- EEG_df %>% 
#   slice(1:40000)
# 
# tmp_df <- data.frame(EEG_scaled %>% ts())
# tmps_small <- tmp_df %>% 
#   slice(1:100000)

### Wavelet -------------------------------------------------------------------

library(WaveletComp)
wave <- analyze.wavelet(EEG_ds_df %>% 
                          dplyr::filter(times > 250, times < 290), "Series.1", 
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
wave$Freq = c(seq(1,1, length.out = length(wave$Period)))/wave$Period


wt.image(wave, color.key = "interval", n.levels = 150,
         legend.params = list(lab = "wavelet power levels", mar = 4.7),
         plot.coi = F,
         plot.contour = F,
         plot.ridge = F
)


powers <- wave$Power
powers %>% dim()
freqs <- (1/wave$Period) %>% round(2)
period <- (wave$Period) %>%  round(2)


library(reshape2)
library(ggplot2)

powers <- powers %>% 
  melt() #%>% 
  #mutate(freq = rep(freqs, 801))

#powers <- as.tibble(powers)
#p <- gather(as.tibble(powers), key = "col", value = "power", scale)

powers$value %>% length()

powers$Var1[1:142] %>% range()

range(freqs)
range(period)


freqs %>% plot()
-log1p(powers$Var1[1:162]) %>% plot(col = "red")

((powers$Var1[1:162] / freqs)/100) %>% 
  rev() %>% 
  points(col = "red")
((powers$Var1[1:162] / freqs)/100) %>% 
  rev()

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
              mapping = aes(x = Var2/samp_rate_ds + 250,
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
              dplyr::slice(plot_from:plot_to), 
            mapping = aes(x = times, y = -2*Series.1+20), color = "white") + 
  geom_point(data = ap_peaks %>% 
                filter(value > 250.1, value < 290), 
             mapping = aes(x = value, y = 100),
             shape = "|",
             color = "white", size = 4) 

  



powers$Var1 %>% length()  /freqs %>% length()


rep(freqs, 1599) %>% length()


#scale_y_continuous(name = "Frequency (Hz)", breaks = freqs)
#ylim(min(freqs),max(freqs))



ggplot(data = EEG_ds_df %>% 
         dplyr::slice(plot_from:plot_to), mapping = aes(x = times, y = Series.1))+
geom_line()
  







  
  

  
  
  
  
  
  
  
  
mean.power = matrix(0,1,801)

for (i in 1:length(wave$Power[1,])){
  mean.power[1,i] = mean(wave$Power[10:20,i])
}

plot(mean.power[1,], type = 'l')
plot(small[1], col = 'red')

### No. of periods x No. of time points 
library(colorRamps)
which(wave$axis.2 < 0.4 & wave$axis.2 > 0.25)
heatmap(wave$Power[151:153,1:5000], Rowv = NA, Colv = NA, col = matlab.like(256))


heatmap(wave$Power, Rowv = NA, Colv = NA, col = matlab.like(256))

wavelet.cols = matlab.like(7)

library(plotly)
plot_ly(z = wave$Power, colors = wavelet.cols, type = 'heatmap')




  