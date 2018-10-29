downSamp <- function(data, ds_factor, samp_rate) {
  counter <- 0
  
  samp_rate_new <- samp_rate/ds_factor
  r_length_decimal <- samp_rate_new %>% log10() %>% as.integer()
  r_length <- (length(data)/samp_rate) %>% round(r_length_decimal)
  
  
  No_points_orig <- length(data) %>% as.numeric()
  No_points_new <- (r_length * samp_rate_new) %>% as.numeric()
  No_points_diff <- No_points_orig - No_points_new %>% as.numeric()
  # data <- data[-seq(from = 1, to = length(data), length.out = No_points_diff)]
  
  if (missing(samp_rate)){
    samp_rate <- 1
    print("No sampling rate was entered. Could not calculate recording length")
  } else {
    
    print(paste0("Length of your recording: ", r_length, "s"))
  }
  
  if (log2(ds_factor) %% 1 == 0) {
    iter <- log2(ds_factor)
  } else {
    iter <- round(log2(ds_factor))
    warning("Downsampling rate is not the power of 2. Rounding was applied")
  }
  
  repeat{
    data <- data[-(seq(2, length(data), 2))]
    counter <- counter + 1
    if (counter == iter) {
      break
    }
  }
  print(paste0("Number of iterations: ", iter))
  print(paste0("Number of datapoints in your output vector: ", length(data)))
  
  #print(paste0("Sampling rate of your new data: ", length(data)/r_length))
  output <- list()
  output$data <- data
  output$iter <- iter
  output$data_points_orig <- No_points_orig
  output$data_points_new <- No_points_new
  output$data_points_diff <- No_points_diff
  output$r_length <- r_length
  return(output)
}
