### LOADING PACKAGES ############

library(R.matlab)
library(tidyverse)
kelibrary(reshape2)
library(bspec) ### power spectrum
library(WaveletComp) ### wavelet
library(diptest) ### to test distribution uni/multimodality (ISI)

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
initial_value <- RECORDINGS$stim_freq[1]
stim_counter <- 1
index <- 1

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

RECORDINGS$stim_freq %>% as.factor() %>% levels()
RECORDINGS$stim_freq_categ %>% as.factor() %>% levels()

# RECORDINGS <- RECORDINGS %>% dplyr::filter(animal_id != "not specified")

animal_ID_list <- RECORDINGS$animal_id %>% as.factor() %>% levels() %>% as.list()

freq_filter_20 <- c("18", "20")
freq_filter_10 <- c("8", "10", "12")
freq_filter_1 <- "1"
CreatePSTHTibble <- function(animal_id, RECORDINGS, freqs) {

  ### parameters ----------------
  animal_filter <- animal_id
  RECORDINGS <- RECORDINGS
  freq_filter <- freqs
  animal_name <- as.character(animal_id)




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

  file_exist_test <- file.exists(file.path("output_data", "cortical_stim_analysis.csv"))
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
lapply(animal_ID_list, CreatePSTHTibble, RECORDINGS = RECORDINGS, freq_filter_1)



### LOADING STIM_RESULTS (produced by CreatePSTHTibble.R, data for PSTH)
STIM_RESULTS <- read_csv(file.path("output_data", "cortical_stim_analysis.csv"))
STIM_RESULTS$freq %>% as.factor() %>% levels()
STIM_RESULTS$animal_id %>% as.factor() %>% levels()


### Boxplots ----------------------------------------

source(file.path("supplementary_functions", "geom_flat_violin.R"))

gp_box_1 <- ggplot(
  data = STIM_RESULTS %>%
    dplyr::filter(
      first_ap_reltimes > PSTH_range[1],
      first_ap_reltimes < PSTH_range[2],
      freq == 1
    ),
  mapping = aes(y = first_ap_reltimes, x = animal_id)
) +
  # geom_dotplot(binaxis = "y", dotsize = 0.5, stackdir = "up", binwidth = 0.001) +
  coord_flip() +
  # geom_flat_violin() +
  geom_boxplot(width = .1, outlier.colour = NA, position = "dodge") +
  geom_hline(yintercept = 0)

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
  # geom_flat_violin() +
  geom_boxplot(width = .1, outlier.colour = NA, position = "dodge") +
  geom_hline(yintercept = 0)

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
  # geom_flat_violin() +
  geom_boxplot(width = .1, outlier.colour = NA, position = "dodge") +
  geom_hline(yintercept = 0)

source(file.path("supplementary_functions", "multiplot.R"))
multiplot(gp_box_1, gp_box_10, gp_box_20, cols = 3)




### Histograms ----------------------------------------


### PLOT: PSTH all cells -----------
PSTH_range <- c(0, 0.05)
freq_to_plot <- 18 %>% as.character()

### main plot
gp_hist <- ggplot(
  data = STIM_RESULTS %>%
    dplyr::filter(
      first_ap_reltimes > PSTH_range[1],
      first_ap_reltimes < PSTH_range[2],
      freq == freq_to_plot
      # animal_id == "GII_21"
    ),
  mapping = aes(x = first_ap_reltimes)
) +
  xlim(PSTH_range[1], PSTH_range[2]) +
  theme_minimal() +
  geom_histogram(aes(fill = animal_id), binwidth = 0.001) +
  # stat_bin(aes(label = ..count..),
  #   binwidth = 0.001,
  #   geom = "text",
  #   vjust = -.5
  # ) +
  theme(
    text = element_text(size = 18),
    axis.text.x = element_text(angle = 45, hjust = 1.2, vjust = 1.2, size = 18)
  ) +
  scale_fill_brewer(
    palette = "Set3",
    name = "Animal"
    # guide = FALSE
  ) +
  facet_wrap(~animal_id, ncol = 3) +
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
      data = STIM_RESULTS %>%
        dplyr::filter(
          first_ap_reltimes > PSTH_range[1],
          first_ap_reltimes < PSTH_range[2],
          freq == freq_to_plot,
          animal_id == as.character(animal_ID_list[[i]])
        ),
      aes_(xintercept = gp_hist_stats$median[i])
    )

  gp_hist <- gp_hist +
    geom_vline(
      data = STIM_RESULTS %>%
        dplyr::filter(
          first_ap_reltimes > PSTH_range[1],
          first_ap_reltimes < PSTH_range[2],
          freq == freq_to_plot,
          animal_id == as.character(animal_ID_list[[i]])
        ),
      aes_(xintercept = gp_hist_stats$quantile25[i]),
      color = "red"
    )

  gp_hist <- gp_hist +
    geom_vline(
      data = STIM_RESULTS %>%
        dplyr::filter(
          first_ap_reltimes > PSTH_range[1],
          first_ap_reltimes < PSTH_range[2],
          freq == freq_to_plot,
          animal_id == as.character(animal_ID_list[[i]])
        ),
      aes_(xintercept = gp_hist_stats$quantile75[i]),
      color = "red"
    )

  gp_hist <- gp_hist +
    geom_vline(
      data = STIM_RESULTS %>%
        dplyr::filter(
          first_ap_reltimes > PSTH_range[1],
          first_ap_reltimes < PSTH_range[2],
          freq == freq_to_plot,
          animal_id == as.character(animal_ID_list[[i]])
        ),
      aes_(xintercept = gp_hist_stats$quantile75[i] + 1.5 * gp_hist_stats$IQR[i]),
      color = "red",
      linetype = "dashed"
    )

  gp_hist <- gp_hist +
    geom_vline(
      data = STIM_RESULTS %>%
        dplyr::filter(
          first_ap_reltimes > PSTH_range[1],
          first_ap_reltimes < PSTH_range[2],
          freq == freq_to_plot,
          animal_id == as.character(animal_ID_list[[i]])
        ),
      aes_(xintercept = gp_hist_stats$quantile25[i] - 1.5 * gp_hist_stats$IQR[i]),
      color = "red",
      linetype = "dashed"
    )

  gp_hist <- gp_hist +
    geom_point(
      data = STIM_RESULTS %>%
        dplyr::filter(
          first_ap_reltimes > PSTH_range[1],
          first_ap_reltimes < PSTH_range[2],
          freq == freq_to_plot,
          animal_id == as.character(animal_ID_list[[i]])
        ),
      aes_(x = gp_hist_stats$mean[i], y = gp_hist_stats$y_max[i] + 4), size = 2
    )

  gp_hist <- gp_hist +
    geom_point(
      data = STIM_RESULTS %>%
        dplyr::filter(
          first_ap_reltimes > PSTH_range[1],
          first_ap_reltimes < PSTH_range[2],
          freq == freq_to_plot,
          animal_id == as.character(animal_ID_list[[i]])
        ),
      aes_(x = gp_hist_stats$peak[i], y = gp_hist_stats$y_max[i] + 4),
      shape = 25,
      size = 2,
      fill = "red"
    )
}
gp_hist
ggsave(file.path("output_data","PSTH_all_cells.png"),
       width = 12,
       height = 10,
       dpi = 300)



### Calculatin response probability ------------------------
gp_hist_stats$all_count_sum
gp_hist_stats$range_count_sum

### No. APs after stimuli @ given frequencies
No_APs <- STIM_RESULTS %>%
  dplyr::filter(first_ap_reltimes > 0, first_ap_reltimes < 0.05, freq == 18) %>%
  group_by(animal_id) %>%
  summarise(No_APs = length(first_ap_reltimes)) %>%
  pull(No_APs)


### No. stim @ given frequency
No_stim <- RECORDINGS %>%
  dplyr::filter(
    signal_type == "stim", stim_freq_categ == 20
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
  add_column(No_APs = gp_hist_stats$all_count_sum) %>%
  add_column(No_APs_range = gp_hist_stats$range_count_sum)

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
  add_column(No_APs_range = gp_hist_stats$range_count_sum)

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
  add_column(No_APs_range = gp_hist_stats$range_count_sum)



RESP_PROB <- bind_rows(bind_rows(resp_probability, resp_probability10), resp_probability20)

RESP_PROB <- RESP_PROB %>%
  mutate(resp_prob_range = (No_APs_range / No_stim) * 100)

library(ggrepel)


### PLOT: resp brob ----------
ggplot(
  data = RESP_PROB,
  mapping = aes(
    x = as.factor(freq),
    y = resp_prob_range
  )
) +
  theme_minimal() +
  xlab("Stimulus frequency") +
  ylab("Response probability") +
  theme(text = element_text(size = 18),
        axis.text = element_text(size = 18)) +
  geom_boxplot(alpha = 0, width = 0.2, lwd = 1, fatten = 1.2) +
  geom_point(
    fill = "dark gray",
    color = "black",
    size = 2.5,
    shape = 21
  ) +
  ### https://cran.r-project.org/web/packages/ggrepel/vignettes/ggrepel.html
  geom_text_repel(aes(label = animal_id),
    nudge_x = 0.15,
    direction = "y",
    hjust = -0.5,
    segment.size = 0.2
  )
ggsave(file.path("output_data","resp_prob.png"),
       width = 6,
       height = 8,
       dpi = 300)



### Latency and probability within stimulus train -----------------------------

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
    stim_frequency = 1
  )
)


hz1 <- AP_PER_STIM %>%
  filter(animal_id != "GII_26") %>%
  group_by(stim_number) %>%
  summarise(
    No_AP_mean = mean(No_AP),
    No_AP_sd = sd(No_AP),
    SEM = sd(No_AP) / sqrt(length(No_AP))
  ) %>%
  add_column(freq = 1)

hz10 <- AP_PER_STIM %>%
  filter(animal_id != "GII_26") %>%
  group_by(stim_number) %>%
  summarise(
    No_AP_mean = mean(No_AP),
    No_AP_sd = sd(No_AP),
    SEM = sd(No_AP) / sqrt(length(No_AP))
  ) %>%
  add_column(freq = 10)

hz20 <- AP_PER_STIM %>%
  filter(animal_id != "GII_26") %>%
  group_by(stim_number) %>%
  summarise(
    No_AP_mean = mean(No_AP),
    No_AP_sd = sd(No_AP),
    SEM = sd(No_AP) / sqrt(length(No_AP))
  ) %>%
  add_column(freq = 20)

bind_rows(bind_rows(hz1, hz10), hz20)


### PLOT: ap numbers ----------
ggplot(
  data = bind_rows(bind_rows(hz1, hz10), hz20),
  mapping = aes(
    x = stim_number,
    y = No_AP_mean,
    fill = as.factor(freq),
    group = as.factor(freq)
  )
) +
  theme_minimal() +
  theme(text = element_text(size = 18),
        axis.text = element_text(size = 18), 
        legend.text = element_text(size = 18)
        ) +
  xlab("Mean no. APs") +
  ylab("No. stimulus in train") + 
  geom_point(size = 3, shape = 21) +
  geom_line(aes(color = as.factor(freq)), lwd = 1) +
  geom_errorbar(aes(
    ymin = No_AP_mean - SEM,
    ymax = No_AP_mean + SEM,
    color = as.factor(freq),
  ),
  width = .5,
  position = position_dodge(0.1)
  ) +
  scale_x_discrete(limits = as.character(c(1:10))) +
  labs(fill = "Frequency", color = "Frequency" ) 
ggsave(file.path("output_data","ap_numbers.png"),width = 14,height = 8,dpi = 300)


ggplot(
  data = AP_PER_STIM,
  mapping = aes(
    x = stim_number,
    y = No_AP,
    shape = as.factor(animal_id),
    color = as.factor(animal_id),
    linetype = as.factor(animal_id)
  )
) +
  geom_line() +
  geom_point() +
  scale_x_discrete(limits = as.character(c(1:10))) +
  scale_shape_manual(values = 1:nlevels(as.factor(AP_PER_STIM$animal_id)))



