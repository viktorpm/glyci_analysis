### detecting AP clusters ("bursts") based on ISI (distance between local maxima / 2)
### does not handel multiple local maxima

BurstThresholdDetect <- function(hist_data, histbreaks) {
  isihist_threshold <- hist(hist_data[hist_data > 0 & hist_data < 1], 
                            breaks = histbreaks, 
                            plot = F)
  
  isihist_diptest <- hist(hist_data,
                          breaks = histbreaks, 
                          plot = F)
  
  
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
  density <- density(hist_data, na.rm = T)
  diptest_density_log <- diptest::dip.test(log(density$y)[log(density$y) %>% is.finite()], 
                                           simulate.p.value = F)
  
  
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
  
  hist(hist_data, breaks = histbreaks, xlim = c(0,1))
  abline(v = isihist_threshold$mids[max1])
  abline(v = isihist_threshold$mids[max2])
  abline(v = isihist_threshold$mids[((max2 - max1) / 2) + max1], col = "red")
  
  plot(log(density$y), type = "l")
  
  burst_threshold <- isihist_threshold$mids[((max2 - max1) / 2) + max1]
  return(list(bt = burst_threshold, 
              clust = clustered, 
              p_val = diptest$p.value,
              p_val_density_log = diptest_density_log$p.value,
              clustered = clustered,
              clustered_d_log = clustered_d_log))
}