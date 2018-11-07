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


animal_ID_list <- recordings$animal_id %>% as.factor() %>% levels() %>% as.list()
recordings$stim_freq %>% as.factor() %>% levels()



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

lapply(animal_ID_list, CreatePSTHTibble, recordings = recordings, freq_filter_1)

stim_result_tibble <- read_csv(file.path("output_data", "cortical_stim_analysis.csv"))

stim_result_tibble$freq %>% as.factor() %>% levels()




### Boxplots ----------------------------------------

source(file.path("supplementary_functions", "geom_flat_violin.R"))

gp_box_1 <- ggplot(
  data = stim_result_tibble %>%
    dplyr::filter(first_ap_reltimes > PSTH_range[1], 
                  first_ap_reltimes < PSTH_range[2], 
                  freq == 1),
  mapping = aes(y = first_ap_reltimes, x = animal_ID)
) +
  # geom_dotplot(binaxis = "y", dotsize = 0.5, stackdir = "up", binwidth = 0.001) +
  coord_flip() +
  # geom_flat_violin() +
  geom_boxplot(width = .1, outlier.colour = NA, position = "dodge") +
  geom_hline(yintercept = 0)

gp_box_10 <- ggplot(
  data = stim_result_tibble %>%
    dplyr::filter(first_ap_reltimes > PSTH_range[1], 
                  first_ap_reltimes < PSTH_range[2], 
                  freq == 8),
  mapping = aes(y = first_ap_reltimes, x = animal_ID)
) +
  # geom_dotplot(binaxis = "y", dotsize = 0.5, stackdir = "up", binwidth = 0.001) +
  coord_flip() +
  # geom_flat_violin() +
  geom_boxplot(width = .1, outlier.colour = NA, position = "dodge") +
  geom_hline(yintercept = 0)

gp_box_20 <- ggplot(
  data = stim_result_tibble %>%
    dplyr::filter(first_ap_reltimes > PSTH_range[1], 
                  first_ap_reltimes < PSTH_range[2], 
                  freq == 18),
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
    dplyr::filter(first_ap_reltimes > PSTH_range[1],
                  first_ap_reltimes < PSTH_range[2], 
                  freq == freq_to_plot),
  mapping = aes(x = first_ap_reltimes)
) +
  xlim(PSTH_range[1], PSTH_range[2]) +
  # geom_histogram(aes(fill = animal_ID), bins = 50) +
  stat_bin(aes(fill = animal_ID), binwidth = 0.001) +
  scale_fill_brewer(palette = "Set3") +
  facet_wrap(~animal_ID, ncol = 1) +
  geom_vline(xintercept = 0)


### calculating basic statistics from the histogram
source(file.path("supplementary_functions", "HistStats.R"))
gp_hist_1_stats <- HistStats(gp_hist)



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
      aes_(xintercept = gp_hist_1_stats$median[i])
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
      aes_(xintercept = gp_hist_1_stats$quantile25[i]),
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
      aes_(xintercept = gp_hist_1_stats$quantile75[i]),
      color = "red"
    )
}

gp_hist




stim_result_tibble %>%
  dplyr::filter(first_ap_reltimes >= 0, first_ap_reltimes <= 0.05, freq == 18) %>%
  group_by(animal_ID) %>%
  summarise(length(first_ap_reltimes))


### No. stim
recordings %>%
  dplyr::filter(signal_type == "stim", animal_id == "GV_136", stim_freq == 18) %>%
  select(signal_time) %>%
  pull() %>%
  length()



