downSamp2 <- function(data, samp_rate_orig, samp_rate_new) {
  library(surveillance)
  library(dplyr)
  samp_rate_orig <- 20000#samp_rate_orig %>% as.numeric()
  samp_rate_new <- 55#samp_rate_new %>% as.numeric()
  data <- EEG_scaled
  ### to determine the decimal place accuracy of r_length
  r_length_decimal <- samp_rate_new %>% log10() %>% as.integer()

  No_points_orig <- length(data) %>% as.numeric()
  r_length <- (No_points_orig / samp_rate_orig) %>% round(r_length_decimal)
  No_points_new <- (r_length * samp_rate_new) %>% as.numeric()
  No_points_diff <- No_points_orig - No_points_new %>% as.numeric()

  # prime_fact <- surveillance::primeFactors(No_points_diff) %>% round() %>% as.integer()
  # iters <- summary(as.factor(prime_fact)) %>% as.integer()
  # to_remove <- levels(as.factor(prime_fact)) %>% as.integer()

  data <- data[-seq(from = 1, to = length(data), length.out = No_points_diff)]
  
  output <- list()
  output$data <- data
  output$data_points <- length(data)
  output$r_length_ds <- r_length
  return(output)  
  
}

# data <- data[-seq(from = 1, to = length(data), by = 1.002004)]
# for (i in 1:length(iters)) {
#   counter <- 0
#   repeat{
#     data <- data[-(seq(
#       from = 1,
#       to = length(data),
#       by = 1.002004
#     )
#     )                   ]
#     counter <- counter + 1
#     if (counter == iters[i]) {
#       break
#     }
#   }
# }
