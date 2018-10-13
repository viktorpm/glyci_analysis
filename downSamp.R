downSamp <- function(data, ds_rate) {
  counter <- 0
  if (log2(ds_rate) %% 1 == 0) {
    iter <- log2(ds_rate)
  } else {
    iter <- round(log2(ds_rate))
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
  return(data)
}
