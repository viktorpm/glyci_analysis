### calculates basic statistics of the data (from plot)
### requires a ggplot object
### extracts the data the plot was drawn from

HistStats <- function(gplot_object) {
  # gplot_object <- gp_hist
  plot_data <- ggplot_build(gplot_object)$data
  variables <- gplot_object$data %>% colnames()
  print(paste0("Variables in your data: ", paste(variables, collapse = ", ")))
  user_input <- readline(prompt = "Chose variable to process: ")
  group_input <- readline(prompt = "Would you like to group your data? (y/n): ")
  
  if (group_input == "y") {
    group_select <- readline(prompt = "Chose variable to group by: ")
  }
  
  
  peak <- plot_data[[1]] %>%
    group_by(PANEL) %>%
    dplyr::filter(count == max(count)) %>%
    summarise(x = first(x), n = n()) %>%
    .$x
  
  mean <- gplot_object$data %>%
  {if (group_input == "y")
    gplot_object$data %>%
      group_by(eval(parse(text = group_select))) %>%
      summarise(means = mean(eval(parse(text = user_input)))) %>%
      pull(means) else gplot_object$data %>% pull(eval(parse(text = user_input))) %>% mean() }
  
  median <- gplot_object$data %>%
  {if (group_input == "y")
    gplot_object$data %>%
      group_by(eval(parse(text = group_select))) %>%
      summarise(medians = median(eval(parse(text = user_input)))) %>%
      pull(medians) else gplot_object$data %>% pull(eval(parse(text = user_input))) %>% median() }
  
  Q <- gplot_object$data %>%
  {if (group_input == "y")
    gplot_object$data %>%
      group_by(eval(parse(text = group_select))) %>%
      summarise(Qs = list(quantile(eval(parse(text = user_input))))) %>%
      pull(Qs) else gplot_object$data %>% 
      pull(eval(parse(text = user_input))) %>%  
      quantile() %>% 
      list() }
  
  
  Q25 <- lapply(seq_along(Q), function(i) {
    Q[[i]]["25%"]
  }) %>% unlist()
  Q75 <- lapply(seq_along(Q), function(i) {
    Q[[i]]["75%"]
  }) %>% unlist()
  IQR <- (Q75 - Q25) %>% unname()
  
  
  
  all_count_sum <- plot_data[[1]] %>%
    group_by(PANEL) %>%
    summarise(all_count_sum = sum(count)) %>%
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
  #   summarise(range_count_sum = sum(count)) %>%
  #   pull(range_count_sum)
  
  
  range_count_sum <- list()
  
  for (i in seq_along(low_range)) {
    range_count_sum[[i]] <- plot_data[[1]] %>%
      group_by(PANEL) %>%
      dplyr::filter(
        x < high_range[i],
        x > low_range[i]
      ) %>%
      summarise(range_count_sum = sum(count)) %>%
      pull(range_count_sum) %>%
      `[[`(i)
  }
  range_count_sum <- range_count_sum %>% unlist()
  
  
  y_axis_max <- plot_data[[1]] %>%
    group_by(PANEL) %>%
    dplyr::filter(count == max(count)) %>%
    summarise(y = first(y), n = n()) %>%
    pull(y)
  
  No_PANELS <- plot_data[[1]] %>%
    select(PANEL) %>%
    pull() %>%
    as.factor() %>%
    levels() %>%
    length()
  
  output <- list(
    peak = peak, ### x positin of the peak of the histogram
    mean = mean, ### x positin of the mean of the data
    median = median, ### x positin of the median of the data
    quantile25 = unname(Q25), ### x positin of the 25th percentile of the data
    quantile75 = unname(Q75), ### x positin of the 75th percentile of the data
    IQR = IQR, ### IQR of the data
    y_max = y_axis_max, ### maximal y axis value (peak value)
    No_PANELS = No_PANELS, ### No. panels of the ggplot object: facet_wrap()
    all_count_sum = all_count_sum, ### sum of all cunts in the plot data range
    range_count_sum = range_count_sum ### sum of all cunts in the +- 1.5*IQR data range (no outliers)
  )
  
  return(output)
}
