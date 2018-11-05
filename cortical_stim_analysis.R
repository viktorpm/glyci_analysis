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
recordings <- recordings %>% dplyr::filter(animal_id != "not specified")


animal_ID_list <- recordings$animal_id %>% as.factor() %>% levels() %>% as.list()
# available_freqs_list <- recordings$stim_freq %>% as.factor() %>% levels()


freq_filter_20 <- c("17", "18", "19", "20")
freq_filter_10 <- c("9", "10")
freq_filter_1 <- c("1")





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
  if (AP_times$signal_time %>% length() != 0){
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







### plotting -----------------------------
hist.range.lower <- -0.01
hist.range.upper <- 0.04

stim_result_tibble %>% head()



stat_summary(fun.y = mean, geom = "point", shape = "|", size = 2)

library(gridExtra)
install.packages("ggridges")
library(ggridges)




BP <- stim_result_tibble %>%
  dplyr::filter(first_ap_reltimes > -0.01, first_ap_reltimes < 0.04) %>% 
  group_by(animal_ID) %>% 
  summarise(plot = list(boxplot(first_ap_reltimes, plot = F))) %>% data.frame()
BP$plot[[1]]$stats[3,1]
BP$plot

stim_result_tibble %>% 
  filter(animal_ID == "GII-21") %>% 
  mutate(lat_mean = BP$plot[[1]]$stats[3,1]) %>% 
  mutate(lat_lower_range = BP$plot[[1]]$stats[1,1]) %>% 
  mutate(lat_uper_range = BP$plot[[1]]$stats[5,1])
  
  



ggplot(
  data = stim_result_tibble %>%
    dplyr::filter(first_ap_reltimes > -0.01, first_ap_reltimes < 0.04) %>%
    dplyr::filter(freq)
    group_by(animal_ID),
  mapping = aes(y = first_ap_reltimes, x = animal_ID)
) +
  geom_boxplot() +
  coord_flip() +
  geom_hline(yintercept = 0)



scale = c(-0.01, 0.04)


ggplot(
  data = stim_result_tibble %>%
    dplyr::filter(first_ap_reltimes > scale[1], first_ap_reltimes < scale[2], freq == 9),
  mapping = aes(x = first_ap_reltimes)
) +
  xlim(scale[1],scale[2]) +
  #geom_histogram(aes(fill = animal_ID), bins = 50) + 
  stat_bin(aes(fill = animal_ID), bins = 100) +
  scale_fill_brewer(palette = "Set3") +
  # geom_vline(xintercept = BP$plot[[1]]$stats[1,1]) +
  # geom_vline(xintercept = BP$plot[[1]]$stats[2,1]) +
  # geom_vline(xintercept = BP$plot[[1]]$stats[5,1]) +
  facet_wrap(~animal_ID, ncol = 1) +
  geom_vline(xintercept = 0) 

  geom_vline(xintercept = stim_result_tibble %>%
                                 dplyr::filter(first_ap_reltimes > -0.01, 
                                               first_ap_reltimes < 0.04) %>% 
                                 group_by(animal_ID) %>% 
                                 summarise(M = mean(first_ap_reltimes)) %>% 
                                 pull(M)
                               ) +
  geom_rug() 






  

