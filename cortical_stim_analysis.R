library(R.matlab)
library(tidyverse)
library(reshape2)
library(bspec) ### power spectrum
library(WaveletComp) ### wavelet
library(diptest) ### to test distribution uni/multimodality (ISI)


source(file.path("supplementary_functions", "CreateRecTibble.R"))
recordings <- CreateRecTibble(
  AP_times = read_csv(file.path("data", "cortical_stim", "AP_times.csv")),
  stim_times = read_csv(file.path("data", "cortical_stim", "stim_times.csv"))
)

str(recordings)
# recordings <- recordings %>% dplyr::filter(animal_id != "not specified")



recordings$stim_freq %>% as.factor() %>% levels()
animal_ID_list <- recordings$animal_id %>% as.factor() %>% levels() %>% as.list()


freq_filter_20 <- c("18", "20")
freq_filter_10 <- c("8", "10", "12")
freq_filter_1 <- "1"
CreatePSTHTibble <- function(animal_ID, recordings, freqs) {

  ### parameters ----------------
  animal_filter <- animal_ID
  recordings <- recordings
  freq_filter <- freqs
  animal_name <- as.character(animal_ID)




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
  AP_times <- recordings %>%
    dplyr::filter(
      signal_type == "AP",
      animal_id == animal_filter,
      unit_id == 1,
      # stim_freq == freq_filter
      sapply(
        str_split(
          string = recordings$available_freqs,
          pattern = ", "
        ),
        function(x) any(x %in% freq_filter)
      )
    ) %>%
    dplyr::select(signal_time)



  ### filtering stim_times ---------------------
  stim_times <- recordings %>%
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

  file_exist_test <- file.exists(file.path("output_data", "cortical_stim_analysis.csv"))
  write_csv(tibble(
    animal_ID = animal_filter,
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
  #   animal_ID = animal_filter,
  #   first_ap_reltimes = first_ap_reltimes %>%
  #     melt() %>%
  #     pull(value),
  #   freq = freq_filter[1]
  #   )
  # return(PSTH_data)
}
lapply(animal_ID_list, CreatePSTHTibble, recordings = recordings, freq_filter_20)




stim_result_tibble <- read_csv(file.path("output_data", "cortical_stim_analysis.csv"))
stim_result_tibble$freq %>% as.factor() %>% levels()
stim_result_tibble$animal_ID %>% as.factor() %>% levels()
animal_ID_list <- stim_result_tibble$animal_ID %>% as.factor() %>% levels() %>% as.list()


### Boxplots ----------------------------------------

source(file.path("supplementary_functions", "geom_flat_violin.R"))

gp_box_1 <- ggplot(
  data = stim_result_tibble %>%
    dplyr::filter(
      first_ap_reltimes > PSTH_range[1],
      first_ap_reltimes < PSTH_range[2],
      freq == 1
    ),
  mapping = aes(y = first_ap_reltimes, x = animal_ID)
) +
  # geom_dotplot(binaxis = "y", dotsize = 0.5, stackdir = "up", binwidth = 0.001) +
  coord_flip() +
  # geom_flat_violin() +
  geom_boxplot(width = .1, outlier.colour = NA, position = "dodge") +
  geom_hline(yintercept = 0)

gp_box_10 <- ggplot(
  data = stim_result_tibble %>%
    dplyr::filter(
      first_ap_reltimes > PSTH_range[1],
      first_ap_reltimes < PSTH_range[2],
      freq == 8
    ),
  mapping = aes(y = first_ap_reltimes, x = animal_ID)
) +
  # geom_dotplot(binaxis = "y", dotsize = 0.5, stackdir = "up", binwidth = 0.001) +
  coord_flip() +
  # geom_flat_violin() +
  geom_boxplot(width = .1, outlier.colour = NA, position = "dodge") +
  geom_hline(yintercept = 0)

gp_box_20 <- ggplot(
  data = stim_result_tibble %>%
    dplyr::filter(
      first_ap_reltimes > PSTH_range[1],
      first_ap_reltimes < PSTH_range[2],
      freq == 18
    ),
  mapping = aes(y = first_ap_reltimes, x = animal_ID)
) +
  # geom_dotplot(binaxis = "y", dotsize = 0.5, stackdir = "up", binwidth = 0.001) +
  coord_flip() +
  # geom_flat_violin() +
  geom_boxplot(width = .1, outlier.colour = NA, position = "dodge") +
  geom_hline(yintercept = 0)

source(file.path("supplementary_functions", "multiplot.R"))
multiplot(gp_box_1, gp_box_10, gp_box_20, cols = 3)




### Histograms ----------------------------------------

PSTH_range <- c(-0.01, 0.04)
freq_to_plot <- 8 %>% as.character()

### main plot
gp_hist <- ggplot(
  data = stim_result_tibble %>%
    dplyr::filter(
      first_ap_reltimes > PSTH_range[1],
      first_ap_reltimes < PSTH_range[2],
      freq == freq_to_plot
      # animal_ID == "GII_21"
    ),
  mapping = aes(x = first_ap_reltimes)
) +
  xlim(PSTH_range[1], PSTH_range[2]) +
  theme_minimal() +
  geom_histogram(aes(fill = animal_ID), binwidth = 0.001) +
  # stat_bin(aes(label = ..count..),
  #   binwidth = 0.001,
  #   geom = "text",
  #   vjust = -.5
  # ) +
  theme(
    text = element_text(size = 15),
    axis.text.x = element_text(angle = 45, hjust = 1.2, vjust = 1.2)
  ) +
  scale_fill_brewer(
    palette = "Set3",
    name = "Animal"
    # guide = FALSE
  ) +
  facet_wrap(~animal_ID, ncol = 3) +
  geom_vline(xintercept = 0) +
  xlab("Time (s)") +
  ylab("Count")

### calculating basic statistics from the histogram
source(file.path("supplementary_functions", "HistStats.R"))
gp_hist_stats <- HistStats(gplot_object = gp_hist)


### vertical lines
for (i in seq_along(animal_ID_list)) {
  gp_hist <- gp_hist +
    geom_vline(
      data = stim_result_tibble %>%
        dplyr::filter(
          first_ap_reltimes > PSTH_range[1],
          first_ap_reltimes < PSTH_range[2],
          freq == freq_to_plot,
          animal_ID == as.character(animal_ID_list[[i]])
        ),
      aes_(xintercept = gp_hist_stats$median[i])
    )

  gp_hist <- gp_hist +
    geom_vline(
      data = stim_result_tibble %>%
        dplyr::filter(
          first_ap_reltimes > PSTH_range[1],
          first_ap_reltimes < PSTH_range[2],
          freq == freq_to_plot,
          animal_ID == as.character(animal_ID_list[[i]])
        ),
      aes_(xintercept = gp_hist_stats$quantile25[i]),
      color = "red"
    )

  gp_hist <- gp_hist +
    geom_vline(
      data = stim_result_tibble %>%
        dplyr::filter(
          first_ap_reltimes > PSTH_range[1],
          first_ap_reltimes < PSTH_range[2],
          freq == freq_to_plot,
          animal_ID == as.character(animal_ID_list[[i]])
        ),
      aes_(xintercept = gp_hist_stats$quantile75[i]),
      color = "red"
    )

  gp_hist <- gp_hist +
    geom_vline(
      data = stim_result_tibble %>%
        dplyr::filter(
          first_ap_reltimes > PSTH_range[1],
          first_ap_reltimes < PSTH_range[2],
          freq == freq_to_plot,
          animal_ID == as.character(animal_ID_list[[i]])
        ),
      aes_(xintercept = gp_hist_stats$quantile75[i] + 1.5 * gp_hist_stats$IQR[i]),
      color = "red",
      linetype = "dashed"
    )

  gp_hist <- gp_hist +
    geom_vline(
      data = stim_result_tibble %>%
        dplyr::filter(
          first_ap_reltimes > PSTH_range[1],
          first_ap_reltimes < PSTH_range[2],
          freq == freq_to_plot,
          animal_ID == as.character(animal_ID_list[[i]])
        ),
      aes_(xintercept = gp_hist_stats$quantile25[i] - 1.5 * gp_hist_stats$IQR[i]),
      color = "red",
      linetype = "dashed"
    )

  gp_hist <- gp_hist +
    geom_point(
      data = stim_result_tibble %>%
        dplyr::filter(
          first_ap_reltimes > PSTH_range[1],
          first_ap_reltimes < PSTH_range[2],
          freq == freq_to_plot,
          animal_ID == as.character(animal_ID_list[[i]])
        ),
      aes_(x = gp_hist_stats$mean[i], y = gp_hist_stats$y_max[i] + 4), size = 2
    )

  gp_hist <- gp_hist +
    geom_point(
      data = stim_result_tibble %>%
        dplyr::filter(
          first_ap_reltimes > PSTH_range[1],
          first_ap_reltimes < PSTH_range[2],
          freq == freq_to_plot,
          animal_ID == as.character(animal_ID_list[[i]])
        ),
      aes_(x = gp_hist_stats$peak[i], y = gp_hist_stats$y_max[i] + 4),
      shape = 25,
      size = 2,
      fill = "red"
    )
}
gp_hist


gp_hist_stats$all_count_sum
gp_hist_stats$range_count_sum

### No. APs after stimuli @ given frequencies
No_APs <- stim_result_tibble %>%
  dplyr::filter(first_ap_reltimes > 0, first_ap_reltimes < 0.05, freq == 18) %>%
  group_by(animal_ID) %>%
  summarise(No_APs = length(first_ap_reltimes)) %>%
  pull(No_APs)


### No. stim @ given frequency
No_stim <- recordings %>%
  dplyr::filter(
    signal_type == "stim", stim_freq == 18
    | stim_freq == 20
    #| stim_freq == 12
  ) %>%
  group_by(animal_id) %>%
  summarise(No_stim = length(signal_time)) %>%
  pull(No_stim)




resp_probability <- recordings %>%
  dplyr::filter(
    signal_type == "stim", stim_freq == 1
    #| stim_freq == 20 
    #| stim_freq == 12
  ) %>%
  group_by(animal_id) %>%
  summarise(No_stim = length(signal_time)) %>%
  add_column(freq = 1, .after = "animal_id") %>%
  add_column(No_APs = gp_hist_stats$all_count_sum) %>%
  add_column(No_APs_range = gp_hist_stats$range_count_sum)

resp_probability10 <- recordings %>%
  dplyr::filter(
    signal_type == "stim",
    stim_freq == 8
    | stim_freq == 10
    | stim_freq == 12
  ) %>%
  group_by(animal_id) %>%
  summarise(No_stim = length(signal_time)) %>%
  add_column(freq = 10, .after = "animal_id") %>%
  add_column(No_APs = gp_hist_stats$all_count_sum) %>%
  add_column(No_APs_range = gp_hist_stats$range_count_sum)

resp_probability20 <- recordings %>%
  dplyr::filter(
    signal_type == "stim", stim_freq == 18
    | stim_freq == 20
    #| stim_freq == 12
  ) %>%
  group_by(animal_id) %>%
  summarise(No_stim = length(signal_time)) %>%
  add_column(freq = 20, .after = "animal_id") %>%
  add_column(No_APs = gp_hist_stats$all_count_sum) %>%
  add_column(No_APs_range = gp_hist_stats$range_count_sum)



resp_prob <- bind_rows(bind_rows(resp_probability, resp_probability10), resp_probability20)

resp_prob <- resp_prob %>%
  mutate(resp_prob_range = (No_APs_range / No_stim) * 100)

ggplot(
  data = resp_prob,
  mapping = aes(
    x = as.factor(freq),
    y = resp_prob_range
  )
) +
  theme_minimal() +
  xlab("Response probability") +
  ylab("Stimulus frequency") +
  theme(text = element_text(size = 15)) +
  geom_point() +
  geom_boxplot(alpha = 0, width = 0.2)


### Latency and probability within stimulus train ################################


source(file.path("supplementary_functions", "CreateRecTibble.R"))
recordings <- CreateRecTibble(
  AP_times = read_csv(file.path("data", "cortical_stim", "AP_times.csv")),
  stim_times = read_csv(file.path("data", "cortical_stim", "stim_times.csv"))
)

str(recordings)
# recordings <- recordings %>% dplyr::filter(animal_id != "not specified")

recordings <- recordings %>%
  mutate(stim_number = 0)


### numbering stimuli within train ---------------------
initial_value <- recordings$stim_freq[1]
stim_counter <- 1
index <- 1

recordings$stim_number[1] <- stim_counter

repeat{
  if (recordings$stim_freq[index + 1] == initial_value) {
    recordings$stim_number[index + 1] <- stim_counter + 1
    stim_counter <- stim_counter + 1
    index <- index + 1
  } else {
    initial_value <- recordings$stim_freq[index + 1]
    stim_counter <- 1
    recordings$stim_number[index + 1] <- stim_counter
    index <- index + 1
  }

  if (index == length(recordings$stim_number)) {
    break
  }
}
recordings$stim_number[recordings$signal_type == "AP"] <- NA

recordings$animal_id %>% as.factor() %>% levels()


CountAPperStim <- function(animal) {
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
    stim <- recordings %>%
      dplyr::filter(
        stim_freq == 1,
        # stim_freq == 20,
        stim_number == i,
        animal_id == animal
      ) %>%
      pull(signal_time)

    ap <- recordings %>%
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
    mutate(anima_id = animal)

  do.call(rbind, AP_counts)
  return(AP_counts)
}

tmp <- lapply(recordings$animal_id %>% as.factor() %>% levels(), CountAPperStim)



ggplot(
  data =
    do.call(rbind, tmp),
  mapping = aes(
    x = stim_number,
    y = No_AP,
    shape = as.factor(anima_id),
    color = as.factor(anima_id),
    linetype = as.factor(anima_id)
  )
) +
  geom_line() +
  geom_point() +
  scale_x_discrete(limits = as.character(c(1:10))) +
  scale_shape_manual(values = 1:nlevels(as.factor(do.call(rbind, tmp)$anima_id)))







# PSTH_range <- c(0, 0.05)
#
# par(mfrow = c(5,2))
# par(mar=c(1,1,1,1))
#
# for (i in c(1:10)) {
#   hist(numbered_stim[[i]][numbered_stim[[i]] > 0 & numbered_stim[[i]] < 0.05 ],
#        xlim = PSTH_range,
#        breaks = 20,
#        main = NULL,
#        xlab = NULL,
#        ylab = NULL)
# }




ggplot(
  data = as.tibble(RTM_1st_ap %>% melt()),
  mapping = aes(x = value)
) +
  xlim(PSTH_range[1], PSTH_range[2]) +
  geom_histogram(binwidth = 0.001) +
  # stat_bin(aes(label = ..count..),
  #   binwidth = 0.001,
  #   geom = "text",
  #   vjust = -.5
  # ) +
  # scale_fill_brewer(palette = "Set3") +
  # facet_wrap(~animal_ID, ncol = 2) +
  geom_vline(xintercept = 0)
