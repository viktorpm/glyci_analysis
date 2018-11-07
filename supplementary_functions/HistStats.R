HistStats <- function(gplot_object) {
  # gplot_object <- egxample_plot
  plot_data <- ggplot_build(gplot_object)$data
  
  peak <- plot_data[[1]] %>% 
    group_by(PANEL) %>%
    dplyr::filter(count == max(count)) %>%
    .$x
  
  mean <- gplot_object$data %>%
    group_by(animal_ID) %>%
    summarise(means = mean(first_ap_reltimes)) %>%
    pull(means)
  
  median <- gplot_object$data %>%
    group_by(animal_ID) %>%
    summarise(medians = median(first_ap_reltimes)) %>%
    pull(medians)
  
  Q <- gplot_object$data %>%
    group_by(animal_ID) %>%
    summarise(Qs = list(quantile(first_ap_reltimes))) %>%
    pull(Qs)
  
  Q25 <- lapply(seq_along(Q), function(i) {
    Q[[i]]["25%"]
  }) %>% unlist()
  Q75 <- lapply(seq_along(Q), function(i) {
    Q[[i]]["75%"]
  }) %>% unlist()
  IQR <- (Q75 - Q25) %>% unname()
  
  output <- list(
    peak = peak,
    mean = mean,
    median = median,
    quantile25 = Q25,
    quantile75 = Q75,
    IQR = IQR
  )
  
  return(output)
}





