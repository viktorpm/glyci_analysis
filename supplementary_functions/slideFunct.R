slideFunct <- function(data, window, step, type) {
  if (missing(type)) {
    type <- "mean"
  }

  total <- length(data)
  spots <- seq(from = 1, to = (total - window), by = step)
  result <- vector(length = length(spots))

  ### calculate sd
  if (type == "sd") {
    for (i in 1:length(spots)) {
      result[i] <- sd(data[spots[i]:(spots[i] + window)])
    }
    return(result)
  }
  
  ### 
  if (type == "mean") {
    for(i in 1:length(spots)){
      result[i] <- mean(data[spots[i]:(spots[i]+window)])
    }
    return(result)
  }
  
  if (type == "cumsum") {
    for(i in 1:length(spots)){
      result[i] <- cumsum(data[spots[i]:(spots[i]+window)])
    }
    return(result)
  }
  
  
}
