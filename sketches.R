a <- c(1:20)

a <- a[-(seq(2, length(a), 2))] 


source((file.path("downSamp.R")))
eeg.ds.new <- downSamp(data = EEG, 512)


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





