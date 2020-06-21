readLines(con = "firing_rate_interaction_plot.R") %>% grep(pattern = "library", x = ., value = T)


lines <- map(list.files(pattern = ".R"),
             ~ readLines(con = .x) %>% 
               grep(pattern = "library\\(", x = ., value = T)
             ) %>%
  unlist()
  
lines %>% substr(start = 1,
                 stop = regexpr(text = ., pattern = "\\)") %>%
                   as.vector()
                 ) # regexpr: in case of multiple matches it only takes the first one




map(list.files(pattern = ".R"),
    ~ readLines(con = .x) %>% 
      grep(pattern = "library\\(", x = ., value = T)
) %>%
  unlist() %>%
  substr(start = 1,
         stop = regexpr(text = ., pattern = "\\)") %>% # regexpr: in case of multiple matches it only takes the first one
           as.vector()
  ) %>% `[` (!duplicated(.)) 
  


  

gregexpr(pattern = "\\)") %>% unlist()
  regexpr(pattern = "\\)") 

  gregexpr(pattern = "\\)")
  
  
  
  
  
  

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


###Hilbert ----

library(seewave)

hilbert(EEG_ds_scaled %>%
          ts(
            start = 0,
            end = rec_length_ds - interval_ds,
            deltat = 1 / samp_rate_ds
          ), samp_rate_ds) %>% plot()


EEG_ds_scaled %>%
  ts(
    start = 0,
    end = rec_length_ds - interval_ds,
    deltat = 1 / samp_rate_ds
  ) %>%
  plot(type = "l", col = "gray",xlim = c (0,20)) 


Im(hilbert(EEG_ds_scaled, samp_rate_ds)) %>% 
  ts(start = 0,
     end = rec_length_ds - interval_ds,
     deltat = 1 / samp_rate_ds) %>% 
  points(xlim = c (0,20), type = "l")


##### multiple gausian separatin


library(mixtools)


h_freq <- hist(ap_peaks %>% diff(), breaks = "FD", freq = T)
h_dens <- hist(ap_peaks %>% diff(), breaks = "FD", freq = F)
h_dens$density <- h_dens$density/100

plot(h_dens,freq = F)

d <- h_dens$density


par(mfrow = c(1,2))
hist(ap_peaks %>% diff(), breaks = "FD", freq = T)
hist(ap_peaks %>% diff(), breaks = "FD", freq = F)
par(mfrow = c(1,1))


fit <- normalmixEM(d, k = 2)
summary(fit)

plot(h_dens,freq = F)


#show the respective curves
lines(d,fit$lambda[1]*dnorm(d,fit$mu[1],fit$sigma[1]), col = "green")
lines(d,fit$lambda[2]*dnorm(d,fit$mu[2],fit$sigma[2]), col = "red")


gg_freq <- ggplot(data = tibble(isi = ap_peaks %>% diff()),
                  mapping = aes(x= isi)) +
  geom_histogram(bins = 100) +
  ggtitle("Frequency histogram")

gg_rel_freq_percent <- ggplot(data = tibble(isi = ap_peaks %>% diff()),
                              mapping = aes(x= isi)) +
  geom_histogram(bins = 100, aes(y = ..density..)) +
  ggtitle("Relative frequency histogram (%)")

gg_rel_freq <- ggplot(data = tibble(isi = ap_peaks %>% diff()),
                      mapping = aes(x= isi)) +
  geom_histogram(bins = 100, aes(y = ..count../sum(..count..))) +
  ggtitle("Relative frequency histogram")

gg_kernel_density <- ggplot(data = tibble(isi = ap_peaks %>% diff()),
                            mapping = aes(x= isi)) +
  geom_density() +
  ggtitle("Kernel density estimate (smoothed histogram)")

ggpubr::ggarrange(gg_freq,gg_rel_freq,gg_rel_freq_percent, gg_kernel_density, ncol = 2, nrow = 2)



########################################
### random stim PSTH


ctx_stim_info <- read.csv(file.path("data","cortical_stim","file_info.csv")) %>% as_tibble()

ctx_stim_info$rec_length[1] %>% as.numeric() 



train_starts <- seq(
  from = 1, 
  to = ctx_stim_info$rec_length[1] %>% as.numeric() - 10,
  length.out = 10) * runif(min = 0, max = 1, n = 10) %>%
  sort()


train_freq <- 10
train_length <- 10

train_starts <- runif(min = 1, max = ctx_stim_info$rec_length[1] %>% as.numeric() - 10, n = 10) %>% sort()
train_ends <- train_starts + 1/train_freq * train_length


seq(from = train_starts[1], to = train_ends[1], by = 1/train_freq) 

seq(from = train_starts[1], to = train_ends[1], length.out = 10) %>% diff() %>% round(1)



