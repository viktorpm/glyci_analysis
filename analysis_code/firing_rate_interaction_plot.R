### recording site: IL
### stimulation: glycinergic fibers in the IL around the recorded cells
### to visualize firing rate change before, during and after the stimulation @30Hz
### firing rate: no. of APs/time (length of the stimulus)
### control firing rates (before and after the stimulus): no. of APs in a time window with the same length as the stimulus

########################## Reading and organizing data ##########################

### function: separated by CapitalLetters
### data from external source: separated.by '.'
### variable: separated_by '_'


### NEW VERSION DEVELOPEMENT ------------------------

library(tidyverse)
library(reshape2)
library(ggrepel)
source(file.path("supplementary_functions", "CreateRecTibble.R"))


### LOADING DATA ------------

### stimulus data
IL_stim_firing <- CreateRecTibble(
  AP_times = read_csv(file.path("data", "IL_MFR", "stimulus", "AP_times.csv")),
  stim_times = read_csv(file.path("data", "IL_MFR", "stimulus", "stim_times.csv"))
)
info_stim <- read_csv(file.path("data", "IL_MFR", "stimulus", "file_info.csv"))

### baseline data
IL_baseline_firing <- CreateRecTibble(
  AP_times = read_csv(file.path("data", "IL_MFR", "baseline", "AP_times.csv")),
  stim_times = read_csv(file.path("data", "IL_MFR", "baseline", "stim_times.csv"))
)
info_baseline <- read_csv(file.path("data", "IL_MFR", "baseline", "file_info.csv"))



### CELL INFO to store categories
CELL_INFO <- info_baseline %>%
  select(file_name) %>%
  add_column(cell_id = substr(info_baseline$file_name, 1, 6), .before = "file_name") %>%
  mutate(bl_activity = T) %>%
  mutate(bl_activity = replace(
    bl_activity,
    .$cell_id == "cell05" |
      .$cell_id == "cell10" |
      .$cell_id == "cell12" |
      .$cell_id == "cell13" |
      .$cell_id == "cell14" |
      .$cell_id == "cell15" |
      .$cell_id == "cell21" |
      .$cell_id == "cell24",
    F
  )) %>%
  mutate(pinch = T) %>%
  mutate(pinch = replace(
    pinch,
    .$cell_id == "cell01" |
      .$cell_id == "cell06" |
      .$cell_id == "cell07" |
      .$cell_id == "cell08" |
      .$cell_id == "cell18" |
      .$cell_id == "cell19" |
      .$cell_id == "cell20" |
      .$cell_id == "cell21" |
      .$cell_id == "cell23" |
      .$cell_id == "cell24" |
      .$cell_id == "cell25" |
      .$cell_id == "cell26" |
      .$cell_id == "cell27" |
      .$cell_id == "cell28",
    F
  )) %>%
  mutate(ident = T) %>%
  mutate(ident = replace(
    ident,
    .$cell_id == "cell06" |
      .$cell_id == "cell07" |
      .$cell_id == "cell11" |
      .$cell_id == "cell12" |
      .$cell_id == "cell13" |
      .$cell_id == "cell14" |
      .$cell_id == "cell21" |
      .$cell_id == "cell22" |
      .$cell_id == "cell23" |
      .$cell_id == "cell24" |
      .$cell_id == "cell25" |
      .$cell_id == "cell26" |
      .$cell_id == "cell27",
    F
  )) %>%
  mutate(control = F) %>%
  mutate(control = replace(
    control,
    .$cell_id == "cell25" |
      .$cell_id == "cell26" |
      .$cell_id == "cell27" |
      .$cell_id == "cell28" |
      .$cell_id == "cell29",
    T
  ))


IL_stim_firing$file_name %>%
  unique() %>%
  length()
info_stim$file_name %>% length()

### CALCULATING --------

IL_stim_firing <- left_join(IL_stim_firing,
  info_stim %>%
    select(
      file_name,
      rec_length,
      train_starts,
      train_ends,
      train_length
    ),
  by = "file_name"
)

### calculating stim start times
start_times <- IL_stim_firing %>%
  group_by(file_name) %>%
  summarise(train_starts = unique(train_starts), train_ends = unique(train_ends)) %>%
  mutate(cell_id = substr(info_stim$file_name, 1, 6)) %>%
  group_by(cell_id) %>%
  summarise(train_starts = unique(train_starts)) %>%
  pull(train_starts) %>%
  strsplit(",") %>%
  set_names(substr(info_stim$file_name, 1, 6)) %>%
  lapply(as.numeric)

### calculating stim end times
end_times <- IL_stim_firing %>%
  group_by(file_name) %>%
  summarise(train_starts = unique(train_starts), train_ends = unique(train_ends)) %>%
  mutate(cell_id = substr(info_stim$file_name, 1, 6)) %>%
  group_by(cell_id) %>%
  summarise(train_ends = unique(train_ends)) %>%
  pull(train_ends) %>%
  strsplit(",") %>%
  set_names(substr(info_stim$file_name, 1, 6)) %>%
  lapply(as.numeric)


# cell_list <- substr(info_stim$file_name, 1, 6)
CELL_INFO$cell_id

### function to calculate #AP b/d/a stim
BDACalculator <- function(data, list) {
  cell_list <- list

  train_length <- data %>%
    filter(substr(file_name, 1, 6) == cell_list) %>%
    select(train_length) %>%
    pull() %>%
    `[[`(1)



  if (
    (data %>%
      filter(substr(file_name, 1, 6) == cell_list, signal_type == "AP") %>%
      select(unit_id) %>%
      unique() %>%
      pull() %>%
      length()) == 1
  ) {
    ### In case of one unit in the file:
    AP_times <- data %>%
      filter(substr(file_name, 1, 6) == cell_list, signal_type == "AP") %>%
      select(signal_time) %>%
      pull()

    AP_numbers <- matrix(0, nrow = length(start_times[[cell_list]]), ncol = 3)
    colnames(AP_numbers) <- c("b", "d", "a")

    for (train_num in 1:length(start_times[[cell_list]])) {
      ### before
      AP_numbers[train_num, 1] <- length(AP_times[start_times[[cell_list]][train_num] - train_length
      < AP_times & AP_times < start_times[[cell_list]][train_num]])

      ### during
      AP_numbers[train_num, 2] <- length(AP_times[AP_times > start_times[[cell_list]][train_num] & AP_times < end_times[[cell_list]][train_num]])

      ### after
      AP_numbers[train_num, 3] <- length(AP_times[AP_times > end_times[[cell_list]][train_num] &
        AP_times < end_times[[cell_list]][train_num] + train_length])
    }

    melt(AP_numbers, varnames = c("train", "stim_cond"), value.name = "No_AP") %>%
      as.tibble() %>%
      add_column(cell_id = cell_list) %>%
      add_column(train_length = train_length)
  } else {
    ### in case of multiple units in the file:
    AP_times_1 <- data %>%
      filter(substr(file_name, 1, 6) == cell_list, signal_type == "AP", unit_id == 1) %>%
      select(signal_time) %>%
      pull()

    AP_numbers_1 <- matrix(0, nrow = length(start_times[[cell_list]]), ncol = 3)
    colnames(AP_numbers_1) <- c("b", "d", "a")

    for (train_num in 1:length(start_times[[cell_list]])) {
      ### before
      AP_numbers_1[train_num, 1] <- length(AP_times_1[start_times[[cell_list]][train_num] - train_length
      < AP_times_1 & AP_times_1 < start_times[[cell_list]][train_num]])

      ### during
      AP_numbers_1[train_num, 2] <- length(AP_times_1[AP_times_1 > start_times[[cell_list]][train_num] &
        AP_times_1 < end_times[[cell_list]][train_num]])

      ### after
      AP_numbers_1[train_num, 3] <- length(AP_times_1[AP_times_1 > end_times[[cell_list]][train_num] &
        AP_times_1 < end_times[[cell_list]][train_num] + train_length])
    }


    AP_times_2 <- data %>%
      filter(substr(file_name, 1, 6) == cell_list, signal_type == "AP", unit_id == 2) %>%
      select(signal_time) %>%
      pull()

    AP_numbers_2 <- matrix(0, nrow = length(start_times[[cell_list]]), ncol = 3)
    colnames(AP_numbers_2) <- c("b", "d", "a")

    for (train_num in 1:length(start_times[[cell_list]])) {
      ### before
      AP_numbers_2[train_num, 1] <- length(AP_times_2[start_times[[cell_list]][train_num] - train_length
      < AP_times_2 & AP_times_2 < start_times[[cell_list]][train_num]])

      ### during
      AP_numbers_2[train_num, 2] <- length(AP_times_2[AP_times_2 > start_times[[cell_list]][train_num] &
        AP_times_2 < end_times[[cell_list]][train_num]])

      ### after
      AP_numbers_2[train_num, 3] <- length(AP_times_2[AP_times_2 > end_times[[cell_list]][train_num] &
        AP_times_2 < end_times[[cell_list]][train_num] + train_length])
    }



    if (cell_list %>% substr(5, 6) %>% as.numeric() %>% nchar() == 1) {
      new_name <- paste0(
        "cell0",
        cell_list %>% substr(5, 6) %>% as.numeric() + 1
      )
    }

    if (cell_list %>% substr(5, 6) %>% as.numeric() %>% nchar() == 2) {
      new_name <- paste0(
        "cell",
        cell_list %>% substr(5, 6) %>% as.numeric() + 1
      )
    }

    bind_rows(
      melt(AP_numbers_1, varnames = c("train", "stim_cond"), value.name = "No_AP") %>%
        as.tibble() %>%
        add_column(cell_id = cell_list) %>%
        add_column(train_length = train_length),

      melt(AP_numbers_2, varnames = c("train", "stim_cond"), value.name = "No_AP") %>%
        as.tibble() %>%
        add_column(cell_id = new_name) %>%
        add_column(train_length = train_length)
    )
  } # else
} # foo


b_d_a_MFR <- lapply(CELL_INFO$cell_id, BDACalculator, data = IL_stim_firing) %>%
  bind_rows() %>%
  mutate(FR = No_AP / train_length) %>%
  dplyr::group_by(stim_cond, cell_id) %>%
  summarise(MFR = mean(FR))




IL_baseline_firing <- left_join(IL_baseline_firing,
  info_baseline %>%
    select(file_name, rec_length),
  by = "file_name"
) %>%
  mutate(stim_cond = "baseline")

### function to calculate ISI mean and SD during baseline activity of the cells
SDMeanISI <- function(f_name) {


  # mean_isi <- IL_baseline_firing %>%
  #   filter(file_name == f_name, signal_type == "AP") %>%
  #   select(signal_time) %>%
  #   pull() %>%
  #   diff() %>%
  #   mean() * 1000 %>%
  #   `names<-`(f_name)
  #
  # sd_isi <- IL_baseline_firing %>%
  #   filter(file_name == f_name, signal_type == "AP") %>%
  #   select(signal_time) %>%
  #   pull() %>%
  #   diff() %>%
  #   sd() * 1000 %>%
  #   `names<-`(f_name)
  #
  # return(list(MFR = mean_isi,
  #             SDFR = sd_isi))

  tmp2 <- IL_baseline_firing %>%
    filter(file_name == f_name, signal_type == "AP") %>%
    select(signal_time) %>%
    pull() %>%
    diff() %>%
    mean() * 1000 %>%
      tibble(mean_isi = .)
  tmp2 %>%
    mutate(sd_isi = IL_baseline_firing %>%
      filter(file_name == f_name, signal_type == "AP") %>%
      select(signal_time) %>%
      pull() %>%
      diff() %>%
      sd() * 1000) %>%
    mutate(file_name = f_name)
}

sd_mean_isi <- lapply(info_baseline$file_name, SDMeanISI) %>%
  bind_rows() %>%
  as.tibble()

sd_mean_isi <- sd_mean_isi %>%
  mutate(cell_id = substr(sd_mean_isi$file_name, 1, 6)) %>%
  mutate(MFR = 1000 / mean_isi) %>%
  mutate(SD_FR = 1000 / sd_isi) %>%
  mutate(MFR_inc = MFR + SD_FR) %>%
  mutate(MFR_dec = MFR - SD_FR) %>%
  # mutate(MFR_inc = 1000/(mean_isi-sd_isi)) %>%
  # mutate(MFR_dec = 1000/(mean_isi+sd_isi)) %>%
  mutate(stim_cond = "baseline") %>%
  mutate(MFR_dec = replace(MFR_dec, sd_mean_isi$MFR_dec < 0, 0))

### calculating ranks based on activity change from "baseline" to "during stimulus"

cellranks <- b_d_a_MFR %>%
  group_by(stim_cond) %>%
  mutate(base_MFR = sd_mean_isi$MFR) %>%
  mutate(activity_change = ((MFR - base_MFR) / base_MFR) * 100) %>%
  dplyr::filter(stim_cond == "d") %>%
  mutate(change_rank = ifelse(activity_change > 0,
    rank(activity_change),
    -rank(-activity_change)
  )) %>%
  ungroup()


### binding b_d_a_MFR (firing rate b/d/a stim) with sd_mean_isi table (baseline firing rate of the same 29 neurons), joining with CELL_INFO containing important information of the cells (baseline activity, identified, pinched, control)

TO_PLOT <- bind_rows(
  sd_mean_isi %>%
    select(MFR, cell_id, stim_cond),
  b_d_a_MFR
) %>%
  left_join(CELL_INFO %>%
    select(-file_name),
  by = "cell_id"
  ) %>%
  left_join(cellranks %>%
    select(cell_id, change_rank, activity_change),
  by = "cell_id"
  )



ggplot(
  data = TO_PLOT %>%
    filter(
      control == F,
      stim_cond == "baseline" | stim_cond == "d"
    ),
  mapping = aes(
    x = forcats::fct_relevel(stim_cond, "baseline", "d", "a"),
    y = MFR
  )
) +
  geom_boxplot(width = 0.2, alpha = 0.5) +
  geom_point(aes(fill = change_rank),
    shape = 21,
    # fill = "#EB8104",
    # color = "white",
    size = 4
  ) +
  # stroke = 2) +
  geom_line(
    aes(group = cell_id, col = change_rank)
  ) +
  geom_label_repel(
    data = TO_PLOT %>%
      filter(
        control == F,
        stim_cond == "baseline"
      ),
    mapping = aes(label = cell_id, col = change_rank),
    nudge_x = -.9,
    direction = "y",
    # hjust = 2,
    segment.size = .9
  )





ggplot(data = cellranks,
       mapping = aes(x = change_rank, y = activity_change)) +
  geom_point(size = 3) +
  geom_text_repel(
    aes(label = cell_id),
    direction = "y",
    angle = 45,
    hjust = .5,
    vjust = 1
    ) +
  geom_hline(yintercept = 0, col = "grey", lty = "dashed") +
  theme_minimal()





ggplot(
  data = TO_PLOT %>%
    filter(
      control == F,
      stim_cond == "b" | stim_cond == "d" | stim_cond == "a"
    ),
  mapping = aes(
    x = forcats::fct_relevel(stim_cond, "b", "d", "a"),
    y = MFR
  )
) +
  theme(panel.background = element_rect(fill = 0)) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 20),
    text = element_text(size = 20)
  ) +
  ggtitle("Pinch", ) +
  xlab("Stimulus") +
  ylim(0, 35) +

  geom_boxplot(width = 0.2, alpha = 0.5) +
  geom_point(
    shape = 21,
    fill = "#EB8104",
    # color = "white",
    size = 4
  ) +
  # stroke = 2) +
  geom_line(
    aes(group = cell_id),
    color = "#EB8104"
  ) +
  geom_label_repel(
    data = TO_PLOT %>%
      filter(
        control == F,
        stim_cond == "a"
      ),
    mapping = aes(label = cell_id),
    nudge_x = 0.15,
    direction = "y",
    hjust = -0.5,
    segment.size = .9
  ) +
  facet_wrap(vars(pinch))
