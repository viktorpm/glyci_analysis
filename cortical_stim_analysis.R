### LOADING PACKAGES ############

library(R.matlab)
library(tidyverse)
library(reshape2)
library(bspec) ### power spectrum
library(WaveletComp) ### wavelet
library(diptest) ### to test distribution uni/multimodality (ISI)
library(ggsignif)
library(ggrepel)
library(gghalves) # half violin plot
library(see) # half violin half dot

### LOADING RECORDINGS (tibble with AP and stim times) ###############
source(file.path("supplementary_functions", "CreateRecTibble.R"))
RECORDINGS <- CreateRecTibble(
  AP_times = read_csv(file.path("data", "cortical_stim", "AP_times.csv")),
  stim_times = read_csv(file.path("data", "cortical_stim", "stim_times.csv"))
)

### adding stimulus number within train
RECORDINGS <- RECORDINGS %>%
  mutate(stim_number = 0)

### calculating stimulus number within train
initial_value <- RECORDINGS$stim_freq[1] %>%
  `comment<-`("First value of stim_freq variable. When it changes stimulus counting restarts")

stim_counter <- 1 %>%
  `comment<-`("Counts stimuli in a train")
index <- 1 %>%
  `comment<-`("Tracks the position (index) of stim_freq")

RECORDINGS$stim_number[1] <- stim_counter
repeat{
  if (RECORDINGS$stim_freq[index + 1] == initial_value) {
    RECORDINGS$stim_number[index + 1] <- stim_counter + 1
    stim_counter <- stim_counter + 1
    index <- index + 1
  } else {
    initial_value <- RECORDINGS$stim_freq[index + 1]
    stim_counter <- 1
    RECORDINGS$stim_number[index + 1] <- stim_counter
    index <- index + 1
  }

  if (index == length(RECORDINGS$stim_number)) {
    break
  }
}

### Alternative to repeat loop to number stimuli
# RECORDINGS %>%
#   filter(signal_type == "stim") %>%
#   mutate(
#     lag_sf = lag(stim_freq),
#     asd = ifelse(lag_sf == stim_freq & !is.na(lag_sf), 0, 1),
#     cum_asd = cumsum(asd)
#   ) %>%
#   select(file_name, stim_freq, cum_asd, stim_number) %>%
#   group_by(cum_asd) %>%
#   mutate(stim_number_asd = row_number()) %>%
#   ungroup() %>%
#   View()

### replacing stim_number with NA at "AP"
RECORDINGS$stim_number[RECORDINGS$signal_type == "AP"] <- NA

### categorizing stimuli
RECORDINGS <- RECORDINGS %>%
  mutate(stim_freq_categ = stim_freq) %>%
  mutate(stim_freq_categ = replace(
    x = stim_freq_categ,
    list = (stim_freq_categ == 12 | stim_freq_categ == 8),
    values = 10
  )) %>%
  mutate(stim_freq_categ = replace(
    x = stim_freq_categ,
    list = (stim_freq_categ == 18),
    values = 20
  ))

RECORDINGS$stim_freq %>%
  as.factor() %>%
  levels()
RECORDINGS$stim_freq %>% unique()

RECORDINGS$stim_freq_categ %>%
  as.factor() %>%
  levels()

# RECORDINGS <- RECORDINGS %>% dplyr::filter(animal_id != "not specified")

animal_ID_list <- RECORDINGS$animal_id %>%
  as.factor() %>%
  levels() %>%
  as.list()

freq_filter_1 <- "1"
freq_filter_10 <- c("8", "10", "12")
freq_filter_20 <- c("18", "20")


CreatePSTHTibble <- function(animal_id, RECORDINGS, freqs) {

  ### parameters ----------------
  animal_filter <- animal_id
  RECORDINGS <- RECORDINGS
  freq_filter <- freqs
  animal_name <- as.character(animal_id)

  # browser()


  ### supplementary function to create reltimes matrices
  RTM_creator <- function(ap, stim) {
    ap_vect <- ap %>% as.vector()
    stim_vect <- stim %>% as.vector()

    mat_AP <- matrix(
      ap_vect,
      nrow = length(stim_vect),
      ncol = length(ap_vect), byrow = T
    )

    mat_stim <- matrix(
      stim_times$signal_time,
      nrow = length(stim_vect),
      ncol = length(ap_vect)
    )


    mat_reltimes <- mat_AP - mat_stim
    mat_reltimes[mat_reltimes == 0] <- NA
    return(mat_reltimes)
  }


  ### filtering AP_times ---------------------
  AP_times <- RECORDINGS %>%
    dplyr::filter(
      signal_type == "AP",
      animal_id == animal_filter,
      unit_id == 1,
      # stim_freq == freq_filter
      sapply(
        str_split(
          string = RECORDINGS$available_freqs,
          pattern = ", "
        ),
        function(x) any(x %in% freq_filter)
      )
    ) %>%
    dplyr::select(signal_time)



  ### filtering stim_times ---------------------
  stim_times <- RECORDINGS %>%
    dplyr::filter(
      signal_type == "stim",
      animal_id == animal_filter,
      stim_freq %in% freq_filter
    ) %>%
    dplyr::select(signal_time)




  rel_times <- RTM_creator(
    ap = AP_times %>% pull(signal_time),
    stim = stim_times %>% pull(signal_time)
  )


  ### calculating first APs after stimuli ------------------------------
  if (AP_times$signal_time %>% length() != 0) {
    first_spikes <- matrix(0, length(stim_times$signal_time), 1)
    for (i in 1:length(stim_times$signal_time)) {
      first_spikes[i, 1] <- AP_times %>%
        dplyr::filter(signal_time > stim_times$signal_time[i]) %>%
        pull() %>%
        min()
    }
    first_spikes[is.infinite(first_spikes)] <- NA
  } else {
    first_spikes <- 0
    print(paste0(animal_name, ": ", "No APs in this category"))
  }





  # first_spikes[119, 1] <- AP_times %>%
  #   dplyr::filter(signal_time > stim_times$signal_time[119]) %>%
  #   pull() %>% min()
  #
  # %>% is_empty() %>% ifelse(min(NA),min(.))




  first_ap_reltimes <- RTM_creator(
    ap = first_spikes,
    stim = stim_times %>% pull(signal_time)
  )

  ### constructing tibble to store first AP reltimes ------------------------------

  dir.create("output_data")

  file_exist_test <- file.exists(file.path("output_data", "cortical_stim_analysis.csv"))

  # if (file_exist_test == T) {
  #   isTRUE(freqs %in% read.csv(file.path("output_data", "cortical_stim_analysis.csv"))$freq)
  # }

  write_csv(tibble(
    animal_id = animal_filter,
    first_ap_reltimes = first_ap_reltimes %>%
      melt() %>%
      pull(value),
    freq = freq_filter[1]
  ),
  file.path("output_data", "cortical_stim_analysis.csv"),
  append = T,
  col_names = !file_exist_test
  )


  # PSTH_data <- tibble(
  #   animal_id = animal_filter,
  #   first_ap_reltimes = first_ap_reltimes %>%
  #     melt() %>%
  #     pull(value),
  #   freq = freq_filter[1]
  #   )
  # return(PSTH_data)
}


# freq_filter_10
# freq_filter_1
# freq_filter_20

RECORDINGS %>%
  dplyr::group_by(animal_id, stim_freq_categ, signal_type, stim_number) %>%
  summarise(time_stamps = paste(signal_time, collapse = ", ")) -> tmp

tmp %>%
  ungroup() %>%
  slice(1) %>%
  pull(time_stamps) %>%
  strsplit(split = ", ") %>%
  unlist() %>%
  as.numeric()

### Replace freq_filter_*, CAREFUL!!! running it multiple times cretes duplicates!!!!!
lapply(animal_ID_list,
  CreatePSTHTibble,
  RECORDINGS = RECORDINGS,
  freqs = freq_filter_20
)




### LOADING STIM_RESULTS (produced by CreatePSTHTibble.R, data for PSTH)
STIM_RESULTS <- read_csv(file.path("output_data", "cortical_stim_analysis.csv"))
STIM_RESULTS$freq %>%
  as.factor() %>%
  levels()
STIM_RESULTS$animal_id %>%
  as.factor() %>%
  levels()


### Boxplots ----------------------------------------

source(file.path("supplementary_functions", "geom_flat_violin.R"))
PSTH_range <- c(0, 0.05)
gp_box_1 <- ggplot(
  data = STIM_RESULTS %>%
    dplyr::filter(
      first_ap_reltimes > PSTH_range[1],
      first_ap_reltimes < PSTH_range[2],
      freq == 1
    ),
  mapping = aes(y = first_ap_reltimes, x = animal_id)
) +
  coord_flip() +
  gghalves::geom_half_violin(side = "r", trim = T, scale = "width") +
  # geom_dotplot(binaxis = "y",
  #              dotsize = 4,
  #              stackdir = "down",
  #              binwidth = 0.0002,method = "histodot",binpositions = "all",
  #              ) +
  # see::geom_violindot(binwidth = 0.0001,size_dots = 6,) +
  geom_jitter(position = position_jitter(0), shape = "|") +
  # geom_boxplot(width = .1, outlier.colour = NA, position = "dodge") +
  geom_hline(yintercept = 0) +
  ggtitle("1 Hz")

gp_box_10 <- ggplot(
  data = STIM_RESULTS %>%
    dplyr::filter(
      first_ap_reltimes > PSTH_range[1],
      first_ap_reltimes < PSTH_range[2],
      freq == 8
    ),
  mapping = aes(y = first_ap_reltimes, x = animal_id)
) +
  # geom_dotplot(binaxis = "y", dotsize = 0.5, stackdir = "up", binwidth = 0.001) +
  coord_flip() +
  geom_half_violin(side = "r", trim = T, scale = "width") +
  geom_jitter(position = position_jitter(0), shape = "|") +
  # geom_flat_violin() +
  # geom_boxplot(width = .1, outlier.colour = NA, position = "dodge") +
  geom_hline(yintercept = 0) +
  ggtitle("10 Hz")

gp_box_20 <- ggplot(
  data = STIM_RESULTS %>%
    dplyr::filter(
      first_ap_reltimes > PSTH_range[1],
      first_ap_reltimes < PSTH_range[2],
      freq == 18
    ),
  mapping = aes(y = first_ap_reltimes, x = animal_id)
) +
  # geom_dotplot(binaxis = "y", dotsize = 0.5, stackdir = "up", binwidth = 0.001) +
  coord_flip() +
  geom_half_violin(side = "r", trim = T, scale = "width") +
  geom_jitter(position = position_jitter(0), shape = "|") +
  # geom_flat_violin() +
  # geom_boxplot(width = .1, outlier.colour = NA, position = "dodge") +
  geom_hline(yintercept = 0) +
  ggtitle("20 Hz")

source(file.path("supplementary_functions", "multiplot.R"))
multiplot(gp_box_1, gp_box_10, gp_box_20, cols = 3)


PSTH_range <- c(0, 0.05)
ggplot(
  data = STIM_RESULTS %>%
    dplyr::filter(
      first_ap_reltimes > PSTH_range[1],
      first_ap_reltimes < PSTH_range[2],
      freq == 1
    ),
  mapping = aes(x = first_ap_reltimes)
) +
  geom_histogram() +
  facet_wrap(~animal_id)


### Histograms ----------------------------------------


### PLOT: PSTH all cells -----------
PSTH_range <- c(-0.01, 0.04)

### calculating basic statistics from the histogram
source(file.path("supplementary_functions", "HistStats.R"))

STIM_RESULTS %>%
  dplyr::group_by(freq) %>%
  summarise(n_times = length(unique(animal_id))) %>%
  ungroup() %>%
  as.matrix() %>%
  unname() -> rep_matrix


hist_stats <- STIM_RESULTS %>%
  dplyr::filter(
    first_ap_reltimes > PSTH_range[1],
    first_ap_reltimes < PSTH_range[2]
  ) %>%
  dplyr::group_by(freq) %>%
  do(gp_hist = ggplot(
    data = .,
    mapping = aes(x = first_ap_reltimes)
  ) +
    xlim(PSTH_range[1], PSTH_range[2]) +
    geom_histogram(aes(fill = animal_id), binwidth = 0.001) +
    scale_fill_brewer(
      palette = "Set3",
      name = "Animal"
      # guide = FALSE
    ) +
    facet_wrap(~animal_id, ncol = 3) +
    geom_vline(xintercept = 0) +
    xlab("Time (s)") +
    ylab("Count") +
    ggtitle(paste0(.$freq, " Hz"))) %>%
  `[[`(2) %>%
  lapply(. %>% HistStats() %>% as.tibble()) %>%
  bind_rows() %>%
  add_column(
    freq = c(
      rep(rep_matrix[1, 1], rep_matrix[1, 2]),
      rep(rep_matrix[2, 1], rep_matrix[2, 2]),
      rep(rep_matrix[3, 1], rep_matrix[3, 2])
    ),
    .after = "animal_id"
  )




### main plot

### 1, 8, 18
freq_to_plot <- 8 %>% as.character()

gp_hist <- ggplot(
  data = STIM_RESULTS %>%
    dplyr::filter(
      first_ap_reltimes > PSTH_range[1],
      first_ap_reltimes < PSTH_range[2],
      freq == freq_to_plot
      # animal_id == "GII_21"
    ),
  mapping = aes(x = first_ap_reltimes)) +
  xlim(PSTH_range[1], PSTH_range[2]) +
  # theme_minimal() +
  geom_histogram(aes(fill = animal_id), binwidth = 0.001) +
  geom_jitter(
    height = 3,
    aes(y = first_ap_reltimes + hist_stats %>% 
          dplyr::filter(freq == freq_to_plot) %>% 
          dplyr::group_by(animal_id) %>% 
          pull(y_max),
        col = animal_id
        ),inherit.aes = T
    ) +

  # stat_bin(aes(label = ..count..),
  #   binwidth = 0.001,
  #   geom = "text",
  #   vjust = -.5
  # ) +
  # theme(
  #   text = element_text(size = 18),
  #   axis.text.x = element_text(angle = 45, hjust = 1.2, vjust = 1.2, size = 18)
  # ) +
  scale_fill_brewer(
    palette = "Set3",
    name = "Animal"
    # guide = FALSE
  ) +
  scale_color_brewer(
    palette = "Set3",
    name = "Animal"
    # guide = FALSE
  ) +
  facet_wrap(~animal_id, ncol = 3) +
  geom_vline(xintercept = 0) +
  xlab("Time (s)") +
  ylab("Count")




gp_hist +
  ### median
  geom_segment(
    data = hist_stats %>%
      dplyr::filter(freq == freq_to_plot) %>%
      dplyr::group_by(animal_id),
    mapping = aes(
      x = median,
      xend = median,
      y = y_max + 7,
      yend = y_max + 18
    ),
    size = 1
  ) +

  ### data range
  geom_segment(
    data = hist_stats %>%
      dplyr::filter(freq == freq_to_plot) %>%
      dplyr::group_by(animal_id),
    mapping = aes(
      x = quantile25 - 1.5 * IQR,
      xend = quantile25,
      y = y_max + 12.5,
      yend = y_max + 12.5
    )
  ) +

  geom_segment(
    data = hist_stats %>%
      dplyr::filter(freq == freq_to_plot) %>%
      dplyr::group_by(animal_id),
    mapping = aes(x = quantile75, xend = quantile75 + 1.5 * IQR, y = y_max + 12.5, yend = y_max + 12.5)
  ) +

  ### IQR (box)
  geom_rect(
    data = hist_stats %>%
      dplyr::filter(freq == freq_to_plot) %>%
      dplyr::group_by(animal_id),
    mapping = aes(xmin = quantile25, xmax = quantile75, ymin = y_max + 7, ymax = y_max + 18),
    color = "red",
    alpha = 0,
    inherit.aes = F
  ) +

  ### mean
  geom_point(
    data = hist_stats %>%
      dplyr::filter(freq == freq_to_plot) %>%
      dplyr::group_by(animal_id),
    mapping = aes(x = mean, y = y_max + 4), size = 2
  ) +

  ### peak
  geom_point(
    data = hist_stats %>%
      dplyr::filter(freq == freq_to_plot) %>%
      dplyr::group_by(animal_id),
    mapping = aes(x = peak, y = y_max + 4),
    shape = 25,
    size = 2,
    fill = "red"
  ) 









# gp_hist_stats <- HistStats(gplot_object = gp_hist)


### vertical lines
# for (i in seq_along(animal_ID_list)) {
#   gp_hist <- gp_hist +
#     geom_vline(
#       data = STIM_RESULTS %>%
#         dplyr::filter(
#           first_ap_reltimes > PSTH_range[1],
#           first_ap_reltimes < PSTH_range[2],
#           freq == freq_to_plot,
#           animal_id == as.character(animal_ID_list[[i]])
#         ),
#       aes_(xintercept = gp_hist_stats$median[i])
#     )
#
#   gp_hist <- gp_hist +
#     geom_vline(
#       data = STIM_RESULTS %>%
#         dplyr::filter(
#           first_ap_reltimes > PSTH_range[1],
#           first_ap_reltimes < PSTH_range[2],
#           freq == freq_to_plot,
#           animal_id == as.character(animal_ID_list[[i]])
#         ),
#       aes_(xintercept = gp_hist_stats$quantile25[i]),
#       color = "red"
#     )
#
#   gp_hist <- gp_hist +
#     geom_vline(
#       data = STIM_RESULTS %>%
#         dplyr::filter(
#           first_ap_reltimes > PSTH_range[1],
#           first_ap_reltimes < PSTH_range[2],
#           freq == freq_to_plot,
#           animal_id == as.character(animal_ID_list[[i]])
#         ),
#       aes_(xintercept = gp_hist_stats$quantile75[i]),
#       color = "red"
#     )
#
#   gp_hist <- gp_hist +
#     geom_vline(
#       data = STIM_RESULTS %>%
#         dplyr::filter(
#           first_ap_reltimes > PSTH_range[1],
#           first_ap_reltimes < PSTH_range[2],
#           freq == freq_to_plot,
#           animal_id == as.character(animal_ID_list[[i]])
#         ),
#       aes_(xintercept = gp_hist_stats$quantile75[i] + 1.5 * gp_hist_stats$IQR[i]),
#       color = "red",
#       linetype = "dashed"
#     )
#
#   gp_hist <- gp_hist +
#     geom_vline(
#       data = STIM_RESULTS %>%
#         dplyr::filter(
#           first_ap_reltimes > PSTH_range[1],
#           first_ap_reltimes < PSTH_range[2],
#           freq == freq_to_plot,
#           animal_id == as.character(animal_ID_list[[i]])
#         ),
#       aes_(xintercept = gp_hist_stats$quantile25[i] - 1.5 * gp_hist_stats$IQR[i]),
#       color = "red",
#       linetype = "dashed"
#     )
#
#   gp_hist <- gp_hist +
#     geom_point(
#       data = STIM_RESULTS %>%
#         dplyr::filter(
#           first_ap_reltimes > PSTH_range[1],
#           first_ap_reltimes < PSTH_range[2],
#           freq == freq_to_plot,
#           animal_id == as.character(animal_ID_list[[i]])
#         ),
#       aes_(x = gp_hist_stats$mean[i], y = gp_hist_stats$y_max[i] + 4), size = 2
#     )
#
#   gp_hist <- gp_hist +
#     geom_point(
#       data = STIM_RESULTS %>%
#         dplyr::filter(
#           first_ap_reltimes > PSTH_range[1],
#           first_ap_reltimes < PSTH_range[2],
#           freq == freq_to_plot,
#           animal_id == as.character(animal_ID_list[[i]])
#         ),
#       aes_(x = gp_hist_stats$peak[i], y = gp_hist_stats$y_max[i] + 4),
#       shape = 25,
#       size = 2,
#       fill = "red"
#     )
# }
# gp_hist

# ggsave(file.path("output_data", "PSTH_all_cells.png"),
#   width = 12,
#   height = 10,
#   dpi = 300
# )



### Calculatin response probability ------------------------
hist_stats$all_count_sum
hist_stats$range_count_sum

### No. APs after stimuli @ given frequencies
No_APs <- STIM_RESULTS %>%
  dplyr::filter(first_ap_reltimes > 0, first_ap_reltimes < 0.05, freq == 18) %>%
  group_by(animal_id) %>%
  summarise(No_APs = length(first_ap_reltimes)) %>%
  pull(No_APs)


### No. stim @ given frequency
No_stim <- RECORDINGS %>%
  dplyr::filter(
    signal_type == "stim", stim_freq_categ == 1
  ) %>%
  group_by(animal_id) %>%
  summarise(No_stim = length(signal_time)) %>%
  pull(No_stim)




resp_probability <- RECORDINGS %>%
  dplyr::filter(
    signal_type == "stim", stim_freq_categ == 1
    #| stim_freq == 20 
    #| stim_freq == 12
  ) %>%
  group_by(animal_id) %>%
  summarise(No_stim = length(signal_time)) %>%
  add_column(freq = 1, .after = "animal_id") %>%
  add_column(No_APs = hist_stats %>% dplyr::filter(freq == 1) %>% pull(all_count_sum)) %>%
  add_column(No_APs_range = hist_stats %>% dplyr::filter(freq == 1) %>% pull(range_count_sum)) %>%
  add_column(peak_latency = hist_stats %>% dplyr::filter(freq == 1) %>% pull(peak)) %>%
  add_column(mean_latency = hist_stats %>% dplyr::filter(freq == 1) %>% pull(mean))


  RECORDINGS %>%
  dplyr::filter(
    signal_type == "stim"
    #| stim_freq == 20 
    #| stim_freq == 12
  ) %>%
  group_by(animal_id, stim_freq_categ) %>%
  summarise(No_stim = length(signal_time)) %>% 
  do(No_APs = hist_stats %>% pull(all_count_sum)) %>% 
  pull(No_APs)

  add_column(freq = 1, .after = "animal_id") %>%
  add_column(No_APs = hist_stats %>% dplyr::filter(freq == 1) %>% pull(all_count_sum)) %>%
  add_column(No_APs_range = hist_stats %>% dplyr::filter(freq == 1) %>% pull(range_count_sum)) %>%
  add_column(peak_latency = hist_stats %>% dplyr::filter(freq == 1) %>% pull(peak)) %>%
  add_column(mean_latency = hist_stats %>% dplyr::filter(freq == 1) %>% pull(mean))



resp_probability10 <- RECORDINGS %>%
  dplyr::filter(
    signal_type == "stim",
    stim_freq_categ == 10
    # | stim_freq == 10
    # | stim_freq == 12
  ) %>%
  group_by(animal_id) %>%
  summarise(No_stim = length(signal_time)) %>%
  add_column(freq = 10, .after = "animal_id") %>%
  add_column(No_APs = gp_hist_stats$all_count_sum) %>%
  add_column(No_APs_range = gp_hist_stats$range_count_sum) %>%
  add_column(peak_latency = gp_hist_stats$peak) %>%
  add_column(mean_latency = gp_hist_stats$mean)

resp_probability20 <- RECORDINGS %>%
  dplyr::filter(
    signal_type == "stim",
    stim_freq_categ == 20
    #| stim_freq == 12
  ) %>%
  group_by(animal_id) %>%
  summarise(No_stim = length(signal_time)) %>%
  add_column(freq = 20, .after = "animal_id") %>%
  add_column(No_APs = gp_hist_stats$all_count_sum) %>%
  add_column(No_APs_range = gp_hist_stats$range_count_sum) %>%
  add_column(peak_latency = gp_hist_stats$peak) %>%
  add_column(mean_latency = gp_hist_stats$mean)



RESP_PROB <- bind_rows(resp_probability, resp_probability10, resp_probability20)

RESP_PROB <- RESP_PROB %>%
  mutate(resp_prob_range = (No_APs_range / No_stim) * 100)


# save(RESP_PROB, file = file.path("output_data","RData", "RESP_PROB.RData"))

?RESload()

### PLOT: resp brob ----------
ggplot(
  data = RESP_PROB %>%
    dplyr::filter(animal_id != "GII_20"),
  mapping = aes(
    x = as.factor(freq),
    y = resp_prob_range
  )
) +
  ### theme, design, labels
  theme_bw() +
  xlab("Stimulus frequency") +
  ylab("Response probability") +
  theme(
    text = element_text(size = 35),
    axis.text = element_text(size = 35),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()
  ) +
  scale_fill_brewer(palette = "Reds") +
  guides(fill = FALSE) +

  ### data visualization layers
  geom_point(
    aes(fill = animal_id),
    color = "black",
    size = 5,
    shape = 21,
    position = position_jitter(0.1)
  ) +
  # geom_line(aes(group = animal_id)) +
  geom_boxplot(alpha = 0, width = 0.2, lwd = 1, fatten = 1.2)

### additional texts
### https://cran.r-project.org/web/packages/ggrepel/vignettes/ggrepel.html
# geom_text_repel(aes(label = animal_id),
#   nudge_x = 0.01,
#   direction = "y",
#   hjust = -0.5,
#   segment.size = 0.5
# ) +

### statistics
# geom_signif(comparisons = list(c("10", "20")),
#             map_signif_level = TRUE,
#             test = "mood.test")


ggsave(file.path("output_data", "resp_prob.eps"),
  width = 4,
  height = 8,
  dpi = 300
)


### PLOT: Latencies ----------

ggplot(
  data = bind_rows(
    RESP_PROB %>% select(-peak_latency) %>%
      rename(latency_value = mean_latency) %>%
      mutate(latency_type = "mean") %>%
      dplyr::filter(animal_id != "GII_20"),
    RESP_PROB %>% select(-mean_latency) %>%
      rename(latency_value = peak_latency) %>%
      mutate(latency_type = "peak") %>%
      dplyr::filter(animal_id != "GII_20")
  ),
  mapping = aes(
    x = as.factor(freq),
    y = latency_value,
    fill = latency_type,
    color = latency_type
  )
) +

  ### theme, design, labels
  theme_bw() +
  theme(
    text = element_text(size = 35),
    axis.text = element_text(size = 35),
    axis.ticks.length = unit(5, "mm"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()
  ) +
  # scale_color_brewer(palette = "Blacks") +
  xlab("Stimulus frequency") +
  ylab("Peak latency") +
  scale_y_continuous(
    name = "Latency (ms)",
    breaks = c(0.000, 0.005, 0.010, 0.015, 0.020),
    labels = c(0, 5, 10, 15, 20),
    limits = c(0, 0.022)
  ) +
  labs(fill = "", color = "") +
  # guides(color = T) +
  scale_fill_discrete(
    name = "Latency",
    breaks = c("mean", "peak"),
    labels = c("Mean", "Peak")
  ) +
  scale_color_discrete(
    name = "Latency",
    breaks = c("mean", "peak"),
    labels = c("Mean", "Peak")
  ) +

  ### data visualization layers
  geom_point(
    color = "black",
    size = 5,
    shape = 21,
    position = position_dodge(0.4)
  ) +
  geom_boxplot(
    position = position_dodge(0.4),
    width = 0.3,
    alpha = 0,
    lwd = 1,
    fatten = 1.2
  )

### https://cran.r-project.org/web/packages/ggrepel/vignettes/ggrepel.html
# geom_text_repel(aes(label = animal_id),
#                 nudge_x = 0.15,
#                 direction = "y",
#                 hjust = -0.5,
#                 segment.size = 0.2
# ) +

### Statistics
# facet_wrap("latency_type") +
# stat_compare_means(comparisons = l
#                    )
#
# geom_signif(
#   comparisons = list(c("10", "20")),
#   map_signif_level = TRUE,
#   test = "mood.test")


ggsave(file.path("output_data", "latency.eps"),
  width = 8,
  height = 10,
  dpi = 300
)



### AP number within stimulus train -----------------------------

### counting no APs after each stimuli in a train
CountAPperStim <- function(animal, stim_frequency) {
  numbered_stim <- list()
  RTM_creator <- function(ap, stim) {
    ap_vect <- ap %>% as.vector()
    stim_vect <- stim %>% as.vector()

    mat_AP <- matrix(
      ap_vect,
      nrow = length(stim_vect),
      ncol = length(ap_vect), byrow = T
    )

    mat_stim <- matrix(
      stim_vect,
      nrow = length(stim_vect),
      ncol = length(ap_vect)
    )


    mat_reltimes <- mat_AP - mat_stim
    mat_reltimes[mat_reltimes == 0] <- NA
    return(mat_reltimes)
  }
  for (i in 1:10) {
    stim <- RECORDINGS %>%
      dplyr::filter(
        stim_freq_categ == stim_frequency,
        stim_number == i,
        animal_id == animal
      ) %>%
      pull(signal_time)

    ap <- RECORDINGS %>%
      dplyr::filter(
        signal_type == "AP",
        animal_id == animal
      ) %>%
      pull(signal_time)

    numbered_stim[[i]] <- RTM_creator(ap = ap, stim = stim)
  }


  AP_counts <- list()

  for (i in c(1:10)) {
    AP_counts[[i]] <- length(numbered_stim[[i]][numbered_stim[[i]] > 0 &
      numbered_stim[[i]] < 0.05 ])
  }

  AP_counts <- AP_counts %>%
    unlist() %>%
    tibble(No_AP = .) %>%
    mutate(stim_number = c(1:10)) %>%
    mutate(animal_id = animal)

  do.call(rbind, AP_counts)
  return(AP_counts)
}

AP_PER_STIM <- do.call(
  rbind,
  lapply(RECORDINGS$animal_id %>% as.factor() %>% levels(),
    CountAPperStim,
    stim_frequency = 20
  )
)


hz1 <- AP_PER_STIM %>%
  # dplyr::filter(animal_id != "GII_26", animal_id != "GII_20") %>%
  group_by(stim_number) %>%
  summarise(
    No_AP_mean = mean(No_AP),
    No_AP_sd = sd(No_AP),
    SEM = sd(No_AP) / sqrt(length(No_AP))
  ) %>%
  add_column(freq = 1)

hz10 <- AP_PER_STIM %>%
  # dplyr::filter(animal_id != "GII_26", animal_id != "GII_20") %>%
  group_by(stim_number) %>%
  summarise(
    No_AP_mean = mean(No_AP),
    No_AP_sd = sd(No_AP),
    SEM = sd(No_AP) / sqrt(length(No_AP))
  ) %>%
  add_column(freq = 10)

hz20 <- AP_PER_STIM %>%
  # dplyr::filter(animal_id != "GII_26", animal_id != "GII_20") %>%
  group_by(stim_number) %>%
  summarise(
    No_AP_mean = mean(No_AP),
    No_AP_sd = sd(No_AP),
    SEM = sd(No_AP) / sqrt(length(No_AP))
  ) %>%
  add_column(freq = 20)

AP_COUNT <- bind_rows(hz1, hz10, hz20) ### GII26 and GII20 were removed
AP_COUNT_ALL <- bind_rows(hz1, hz10, hz20) ### with GII26 and GII20
save(AP_COUNT, file = file.path("output_data", "RData", "AP_COUNT.RData"))
save(AP_COUNT_ALL, file = file.path("output_data", "RData", "AP_COUNT_ALL.RData"))



### PLOT: ap numbers ----------
ggplot(
  data = AP_COUNT_ALL,
  mapping = aes(
    x = stim_number,
    y = No_AP_mean,
    fill = as.factor(freq),
    group = as.factor(freq)
  )
) +

  ### theme, design, labels
  theme_minimal() +
  theme(
    text = element_text(size = 20),
    axis.text = element_text(size = 20),
    legend.text = element_text(size = 20),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()
  ) +
  xlab("No. stimulus in train") +
  ylab("Mean no. APs") +
  scale_x_discrete(limits = as.character(c(1:10))) +
  labs(fill = "Frequency", color = "Frequency") +
  # scale_fill_brewer(type = "div", palette = "Reds") +
  # scale_color_brewer(type = "div", palette = "Reds") +

  ### data visualization layers
  geom_point(aes(col = as.factor(freq)),
    size = 3,
    # shape = 21
  ) +
  geom_line(aes(color = as.factor(freq)), lwd = 1) +
  geom_errorbar(aes(
    ymin = No_AP_mean - SEM,
    ymax = No_AP_mean + SEM,
    color = as.factor(freq)
  ),
  width = .5,
  position = position_dodge(0.3),
  lty = 5
  )


ggsave(file.path("output_data", "ap_numbers.eps"),
  width = 12,
  height = 8,
  device = "eps"
  # dpi = 300
)


ggplot(
  data = AP_PER_STIM,
  mapping = aes(
    x = stim_number,
    y = No_AP,
    # shape = as.factor(animal_id),
    color = as.factor(animal_id),
    linetype = as.factor(animal_id)
  )
) +
  geom_line() +
  geom_point(size = 3) +
  scale_x_discrete(limits = as.character(c(1:10))) +
  # scale_shape_manual(values = 1:nlevels(as.factor(AP_PER_STIM$animal_id))) +
  theme_minimal() +
  geom_label_repel(
    data = AP_PER_STIM %>%
      dplyr::filter(stim_number == 10),
    aes(label = animal_id),
    direction = "y",
    hjust = -1
  )
