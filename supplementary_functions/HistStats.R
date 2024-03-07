### calculates basic statistics of the data (from plot)
### requires a ggplot object
### extracts the data the plot was drawn from

HistStats <- function(gplot_object) {
  # browser()
  # gplot_object <- gp_hist
  plot_data <- ggplot_build(gplot_object)$data
  
  animal_id <- gplot_object$data$animal_id %>% unique()

  peak <- plot_data[[1]] %>%
    group_by(PANEL) %>%
    dplyr::filter(count == max(count)) %>%
    dplyr::summarise(x = first(x), n = n()) %>%
    .$x

  mean <- gplot_object$data %>%
    group_by(animal_id) %>%
    dplyr::filter(first_ap_reltimes > 0) %>% 
    dplyr::summarise(means = mean(first_ap_reltimes)) %>%
    pull(means)

  median <- gplot_object$data %>%
    group_by(animal_id) %>%
    dplyr::filter(first_ap_reltimes > 0) %>% 
    dplyr::summarise(medians = median(first_ap_reltimes)) %>%
    pull(medians)

  Q <- gplot_object$data %>%
    group_by(animal_id) %>%
    dplyr::filter(first_ap_reltimes > 0) %>% 
    dplyr::summarise(Qs = list(quantile(first_ap_reltimes))) %>%
    pull(Qs)

  Q25 <- lapply(seq_along(Q), function(i) {
    Q[[i]]["25%"]
  }) %>% unlist()
  Q75 <- lapply(seq_along(Q), function(i) {
    Q[[i]]["75%"]
  }) %>% unlist()
  IQR <- (Q75 - Q25) %>% unname()



  all_count_sum <- plot_data[[1]] %>%
    group_by(PANEL) %>%
    dplyr::summarise(all_count_sum = sum(count)) %>%
    pull(all_count_sum)


  low_range <- replace(
    x = c(Q25 - (1.5 * IQR)),
    list = c(Q25 - (1.5 * IQR)) < 0,
    values = 0
  ) %>% unname()
  high_range <- (Q75 + (1.5 * IQR)) %>% unname()

  
  
  # range_count_sum <- plot_data[[1]] %>%
  #   dplyr::filter(
  #     x < lapply(
  #       seq_along(as.list(high_range)),
  #       function(i) {
  #         (as.list(high_range)[[i]])
  #       }
  #     ),
  #     x > lapply(seq_along(as.list(low_range)), function(i) {
  #       (as.list(low_range)[[i]])
  #     })
  #   ) %>%
  #   group_by(PANEL) %>%
  #   dplyr::summarise(range_count_sum = sum(count)) %>%
  #   pull(range_count_sum)
  
  
  range_count_sum = list()
  
  for (i in seq_along(low_range)){
    range_count_sum[[i]] <- plot_data[[1]] %>%
      group_by(PANEL) %>%
      dplyr::filter(
        x < high_range[i],
        x > low_range[i]
        
      ) %>% 
      dplyr::summarise(range_count_sum = sum(count)) %>%
      pull(range_count_sum) %>% `[[`(i)
      
  }
  range_count_sum <- range_count_sum %>% unlist()


  y_axis_max <- plot_data[[1]] %>%
    group_by(PANEL) %>%
    dplyr::filter(count == max(count)) %>%
    dplyr::summarise(y = first(y), n = n()) %>%
    pull(y)

  No_PANELS <- plot_data[[1]] %>%
    select(PANEL) %>%
    pull() %>%
    as.factor() %>%
    levels() %>%
    length()

  output <- list(
    animal_id = animal_id,
    #stim_freq = frequency,
    peak = peak, ### x positin of the peak of the histogram
    mean = mean, ### x positin of the mean of the data
    median = median, ### x positin of the median of the data
    quantile25 = Q25, ### x positin of the 25th percentile of the data
    quantile75 = Q75, ### x positin of the 75th percentile of the data
    IQR = IQR, ### IQR of the data
    y_max = y_axis_max, ### maximal y axis value (peak value)
    No_PANELS = No_PANELS, ### No. panels of the ggplot object: facet_wrap()
    all_count_sum = all_count_sum, ### sum of all counts in the plot data range
    range_count_sum = range_count_sum ### sum of all counts in the +- 1.5*IQR data range (no outliers)
  )

  return(output)
}
