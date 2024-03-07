### detecting AP clusters ("bursts") based on ISI (distance between local maxima / 2)
### does not handel multiple local maxima

BurstThresholdDetect <- function(hist_data, histbreaks, isi_range) {
  # browser()
  ### filtering data 
  hist_data_filt <- hist_data[hist_data > isi_range[1] & hist_data < isi_range[2]]
  
  isihist_threshold <- hist(hist_data_filt, 
                            breaks = histbreaks, 
                            plot = F)
  
  isihist_diptest <- hist(hist_data_filt,
                          breaks = histbreaks, 
                          plot = F)
  
  # browser()
  ### indices of maximal values in the first and second half of the histogram (max1 and max2)
  max1 <- which(
    isihist_threshold$counts[1:(length(isihist_threshold$counts) / 2)] ==
      isihist_threshold$counts[1:(length(isihist_threshold$counts) / 2)] %>% max()
  ) %>% as.numeric()
  
  max2 <- which(
    isihist_threshold$counts[(length(isihist_threshold$counts) / 2 + 1):
                               length(isihist_threshold$counts)] ==
      isihist_threshold$counts[(length(isihist_threshold$counts) / 2 + 1):
                                 length(isihist_threshold$counts)] %>% max()
  ) + length(isihist_threshold$counts) / 2
  
  diptest <- diptest::dip.test(isihist_diptest$counts, simulate.p.value = F)
  density <- density(hist_data_filt, na.rm = T)
  diptest_density_log <- diptest::dip.test(log(density$y)[log(density$y) %>% is.finite()], 
                                           simulate.p.value = F)
  
  # browser()
  if(diptest$p.value <= 0.05){
    clustered = T
  } else {
    clustered = F
  }
  
  if(diptest_density_log$p.value <= 0.05){
    clustered_d_log = T
  } else {
    clustered_d_log = F
  }
  

  
  plot(log(density$y), type = "l")
  
  ### half way between the max value of the first and second half of the ISI window (range)
  index_treshold1 <- ((max2 - max1) / 2) + max1
  
  ### minimal value between max1 and max2
  index_treshold2 <- which(isihist_threshold$counts[c(max1:max2)] %>% min() == isihist_threshold$counts[c(max1:max2)]) + (max1-1)
  
  
  ### minimum value between threshold1 and max1. This is the best cluster threshold
  index_threshold3 <- which(isihist_threshold$counts[c(max1:index_treshold1)] %>% min() == isihist_threshold$counts[c(max1:index_treshold1)]) + (max1 - 1)
  
  
  burst_threshold <- isihist_threshold$mids[index_treshold1] 
  burst_threshold2 <- isihist_threshold$mids[index_treshold2] 
  burst_threshold3 <- isihist_threshold$mids[index_threshold3]
  
  hist(hist_data_filt, breaks = histbreaks, xlim = isi_range)
  abline(v = isihist_threshold$mids[max1])
  abline(v = isihist_threshold$mids[max2])
  abline(v = isihist_threshold$mids[index_treshold1], col = "red")
  abline(v = isihist_threshold$mids[index_treshold2], col = "magenta")
  abline(v =isihist_threshold$mids[index_threshold3], col = "blue")
  
  # browser()
  if (isihist_threshold$counts[index_threshold3] < isihist_threshold$counts[max2]) {
    clustered = T
  } else {clustered = F}
  
  # browser()
  return(list(bt = burst_threshold,
              bt2 = burst_threshold2,
              bt3 = burst_threshold3,
              # p_val_multipeak = diptest$p.value,
              # p_val_density_log = diptest_density_log$p.value,
              # clustered = clustered %>% as.character(),
              # clustered_d_log = clustered_d_log %>% as.character(),
              clustered = clustered, 
              bins_n = isihist_threshold$breaks %>% length()))
}