a <- c(1:20)

a <- a[-(seq(2, length(a), 2))] 


source((file.path("downSamp.R")))
eeg.ds.new <- downSamp(data = EEG, ds_factor = 512, 20000)

plot(EEG_ds[1:100], type = "l")
points(eeg.ds.new[1:100], type = "l", col = "red")

rm(eeg.ds.new)
log2(9) %>% round()

log2(10) %% 1 == 0

library(tidyverse)

EEG <- raw.rec$EEG[,,1]$values %>% as.double()
EEGDS <- EEG[-(seq(1, length(EEG), 2))]
EEGDS2 <- EEGDS[-(seq(1, length(EEGDS), 2))]
EEGDS3 <- EEGDS2[-(seq(1, length(EEGDS2), 2))]
EEGDS4 <- EEGDS3[-(seq(1, length(EEGDS3), 2))]
EEGDS5 <- EEGDS4[-(seq(1, length(EEGDS4), 2))]
EEGDS6 <- EEGDS5[-(seq(1, length(EEGDS5), 2))]
EEGDS7 <- EEGDS6[-(seq(1, length(EEGDS6), 2))]
EEGDS8 <- EEGDS7[-(seq(1, length(EEGDS7), 2))]
EEGDS9 <- EEGDS8[-(seq(1, length(EEGDS8), 2))]

plot(EEGDS9[1:100], type = "l")
plot(EEG[1:40000], type = "l", col = "red")

length(EEG)

2^9
2*2*2*2*2*2*2*2*2
20000/512

39.0625*rec_length


log2(512) %% 1 == !0

microbenchmark::microbenchmark(
  sqldf("
        SELECT 
        a.*,
        e.sync_start
        FROM ap_peaks AS a
        LEFT JOIN eeg_periods AS e
        ON a.peak_times > e.sync_start
        AND a.peak_times < e.sync_end
        ;
        ") %>% mutate(
          egg_period2 = ifelse(is.na(sync_start), "desync", "sync")
  ),
  for (i in 1:length(eeg_periods$sync_start)){
    ap_peaks <- ap_peaks %>% 
      dplyr::mutate(egg_period = replace(
        egg_period,
        peak_times > eeg_periods$sync_start[i] & peak_times < eeg_periods$sync_end[i],
        values = "sync"))
  } 
)


library(tidyverse)
freq <- c(rep(1, 4), rep(10,3), rep(27, 2))
stim_number <- c(rep(0,length(freq)))


initial_value <- freq[1]
stim_counter <- 1
index <- 1

stim_number[1] <- stim_counter

repeat{
  if (freq[index+1] == initial_value){
    stim_number[index+1] <-  stim_counter+1
    stim_counter <- stim_counter+1
    index <- index + 1
  } else {
    initial_value <- freq[index+1]
    stim_counter <- 1
    stim_number[index+1] <- stim_counter
    index <- index+1
  }
  
  if (index == length(stim_number)){break}
}


source(file.path("supplementary_functions", "CreateRecTibble.R"))
recordings <- CreateRecTibble(
  AP_times = read_csv(file.path("data", "cortical_stim", "AP_times.csv")),
  stim_times = read_csv(file.path("data", "cortical_stim", "stim_times.csv"))
)

str(recordings)
recordings <- recordings %>% dplyr::filter(animal_id != "not specified")

recordings <- recordings %>% 
  mutate(stim_number = 0)

initial_value <- recordings$stim_freq[1]
stim_counter <- 1
index <- 1

recordings$stim_number[1] = stim_counter

repeat{
  if (recordings$stim_freq[index+1] == initial_value){
    recordings$stim_number[index+1] <-  stim_counter+1
    stim_counter <- stim_counter+1
    index <- index + 1
  } else {
    initial_value <- recordings$stim_freq[index+1]
    stim_counter <- 1
    recordings$stim_number[index+1] <- stim_counter
    index <- index+1
  }
  
  if (index == length(recordings$stim_number)){break}
}

recordings$stim_number[recordings$signal_type == "AP"] = NA
