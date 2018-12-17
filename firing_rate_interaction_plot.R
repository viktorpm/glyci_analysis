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
RECORDINGS <- CreateRecTibble(
  AP_times = read_csv(file.path("data", "IL_MFR", "stimulus", "AP_times.csv")),
  stim_times = read_csv(file.path("data", "IL_MFR","stimulus", "stim_times.csv"))
)
info <- read_csv(file.path("data", "IL_MFR", "stimulus", "file_info.csv"))


RECORDINGS$file_name %>% unique() %>% length()
info$file_name %>% length()



RECORDINGS <- left_join(RECORDINGS, info, by = "file_name")


start_times <- RECORDINGS %>%
  group_by(file_name) %>%
  summarise(train_starts = unique(train_starts), train_ends = unique(train_ends)) %>%
  mutate(cell_id = substr(info$file_name, 1, 6)) %>%
  group_by(cell_id) %>%
  summarise(train_starts = unique(train_starts)) %>%
  pull(train_starts) %>%
  strsplit(",") %>%
  set_names(substr(info$file_name, 1, 6)) %>%
  lapply(as.numeric)


end_times <- RECORDINGS %>%
  group_by(file_name) %>%
  summarise(train_starts = unique(train_starts), train_ends = unique(train_ends)) %>%
  mutate(cell_id = substr(info$file_name, 1, 6)) %>%
  group_by(cell_id) %>%
  summarise(train_ends = unique(train_ends)) %>%
  pull(train_ends) %>%
  strsplit(",") %>%
  set_names(substr(info$file_name, 1, 6)) %>%
  lapply(as.numeric)




# split(.,.$cell_id)

cell_list <- substr(info$file_name, 1, 6)

foo <- function(data, list) {
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
    
    
    
    if (cell_list %>% substr(5,6) %>% as.numeric() %>% nchar() == 1) {
      new_name <- paste0(
        "cell0",
        cell_list %>% substr(5,6) %>% as.numeric() + 1
      )
    } 
      
    if (cell_list %>% substr(5,6) %>% as.numeric() %>% nchar() == 2) {
      new_name <- paste0(
        "cell",
        cell_list %>% substr(5,6) %>% as.numeric() + 1
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
   }#else

}#foo


# allCells %>% filter(pinch == F) %>% select(cell) %>% unique()


IL_stim_firing <- lapply(cell_list, foo, data = RECORDINGS) %>% 
  bind_rows() %>% 
  mutate(pinch = T) %>% 
  mutate(pinch = replace(pinch,
                         pinch, .$cell_id == "cell01"|
                           .$cell_id == "cell06"|
                           .$cell_id == "cell18"|
                           .$cell_id == "cell19"|
                           .$cell_id == "cell20"|
                           .$cell_id == "cell21"|
                           .$cell_id == "cell23"|
                           .$cell_id == "cell24"|
                           .$cell_id == "cell25"|
                           .$cell_id == "cell26"|
                           .$cell_id == "cell27"|
                           .$cell_id == "cell28",
                         F)) 




filter_condititon <- "pinch == T" 

ALL_CELLS_PLOT <- IL_stim_firing %>%
  dplyr::filter(eval(parse(text = filter_condititon))) %>% 
  #mutate(No_AP = as.numeric(levels(No_AP))[No_AP]) %>%
  mutate(No_AP_p_stimLength = No_AP / train_length) %>%
  dplyr::group_by(stim_cond, cell_id) %>%
  summarise(MFR = mean(No_AP_p_stimLength))

cells_to_highlight <- IL_stim_firing %>%
  filter(eval(parse(text = filter_condititon))) %>% 
  #mutate(No.AP = as.numeric(levels(No.AP))[No.AP]) %>%
  mutate(No_AP_p_stimLength = No_AP / train_length) %>%
  dplyr::group_by(stim_cond, cell_id) %>%
  summarise(MFR = mean(No_AP_p_stimLength)) %>%
  filter(cell_id == "cell01" )
           # cell == "cell02" |
           # cell == "cell03" |
           # cell == "cell08" |
           # cell == "cell09" |
           # cell == "cell10" |
           # cell == "cell15" |
           # cell == "cell16" |
           # cell == "cell17" |
           # cell == "cell18" |
           # cell == "cell19" |
           # cell == "cell20" )



ggplot(data = ALL_CELLS_PLOT,
       mapping = aes(x = forcats::fct_relevel(stim_cond, "b", "d", "a"),
                     y = MFR)) +
  #theme(panel.background = element_rect(fill = 0)) +
  theme_minimal() +
  theme(axis.text = element_text(size = 20),
        text = element_text(size = 20)) +
  ylim(0,35)+
  
  geom_boxplot(width = 0.2, alpha = 0.5) +
  geom_point(shape = 21, 
             fill = "#EB8104", 
             #color = "white",
             size = 4) +
  #stroke = 2) +
  geom_line(aes(group = cell_id), color = "#EB8104") +
  geom_label_repel(mapping = aes(label = cell_id),
                                    nudge_x = 0.15,
                                    direction = "y",
                                    hjust = -0.5,
                                    segment.size = .9) +
  
  ### Highlight
  # geom_point(data = cells_to_highlight,
  #            #shape = 21,
  #            color = "#1D4871",
  #            #color = "white",
  #            size = 4,
  #            stroke = 2) +
  # geom_line(data = cells_to_highlight, color = "#1D4871", aes(group = cell_id)) +
  # geom_text(data = cells_to_highlight,
  #           aes(label = cell_id),
  #           color = "#1D4871",
  #           hjust = -0.2, vjust = -0.2) +
  # geom_label_repel(data = cells_to_highlight,
  #                  mapping = aes(label = cell_id),
  #                  nudge_x = 0.15,
  #                  direction = "y",
  #                  hjust = -0.5,
  #                  segment.size = .9) +


  scale_x_discrete(name = "Stimulus",labels = c("Before", "During", "After")) +
  labs(y = "Mean firing rate") 



### BASELINE ---------------------

source(file.path("supplementary_functions", "CreateRecTibble.R"))
IL_baseline_firing <- CreateRecTibble(
  AP_times = read_csv(file.path("data", "IL_MFR","baseline", "AP_times.csv")),
  stim_times = read_csv(file.path("data", "IL_MFR","baseline", "stim_times.csv"))
)
info_baseline <- read_csv(file.path("data", "IL_MFR", "baseline","file_info.csv"))

CELL_INFO <- info_baseline %>% 
  select(file_name) %>% 
  add_column(cell_id = substr(info_baseline$file_name, 1,6), .before = "file_name") %>% 
  mutate(bl_activity = T) %>%
  mutate(bl_activity = replace(bl_activity, .$cell_id == "cell05"|
                                 .$cell_id == "cell10"|
                                 .$cell_id == "cell12"|
                                 .$cell_id == "cell13"|
                                 .$cell_id == "cell14"|
                                 .$cell_id == "cell15"|
                                 .$cell_id == "cell21"|
                                 .$cell_id == "cell22",
                               F)) %>% 
  mutate(pinch = T) %>% 
  mutate(pinch = replace(pinch, .$cell_id == "cell01"|
                           .$cell_id == "cell06"|
                           .$cell_id == "cell18"|
                           .$cell_id == "cell19"|
                           .$cell_id == "cell20"|
                           .$cell_id == "cell21"|
                           .$cell_id == "cell23"|
                           .$cell_id == "cell24"|
                           .$cell_id == "cell25"|
                           .$cell_id == "cell26"|
                           .$cell_id == "cell27"|
                           .$cell_id == "cell28",
                         F)) %>% 
  mutate(ident = F) %>% 
  mutate(ident = replace(ident, .$cell_id == "cell01"|
                           .$cell_id == "cell02"|
                           .$cell_id == "cell03"|
                           .$cell_id == "cell08"|
                           .$cell_id == "cell09"|
                           .$cell_id == "cell10"|
                           .$cell_id == "cell15"|
                           .$cell_id == "cell16"|
                           .$cell_id == "cell17"|
                           .$cell_id == "cell18"|
                           .$cell_id == "cell19"|
                           .$cell_id == "cell20"|
                           .$cell_id == "cell28"|
                           .$cell_id == "cell29",
                         T))


IL_baseline_firing <- left_join(IL_baseline_firing, 
                      info_baseline %>%
                        select(file_name,rec_length),
                      by = "file_name") %>% 
  mutate(stim_cond = "baseline")


SDMeanISI <- function(f_name){
  
  
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

sd_mean_isi <- lapply(info_baseline$file_name, SDMeanISI) %>% bind_rows()
sd_mean_isi <- sd_mean_isi %>% 
  mutate(cell_id = substr(sd_mean_isi$file_name,1,6)) %>% 
  mutate(MFR = 1000/mean_isi) %>% 
  mutate(SD_FR = 1000/sd_isi ) %>% 
  mutate(MFR_inc = MFR + SD_FR) %>% 
  mutate(MFR_dec = MFR - SD_FR) %>% 
  # mutate(MFR_inc = 1000/(mean_isi-sd_isi)) %>% 
  # mutate(MFR_dec = 1000/(mean_isi+sd_isi)) %>% 
  mutate(stim_cond = "baseline") %>% 
  mutate(MFR_dec = replace(MFR_dec, sd_mean_isi$MFR_dec < 0, 0))


# info_baseline$No_AP_unit/info_baseline$rec_length

BASE_TO_PLOT <- bind_rows(sd_mean_isi %>%
                            select(MFR, cell_id, stim_cond),
                          ALL_CELLS_PLOT) %>%
  mutate(pinch = T) %>% 
  mutate(pinch = replace(pinch,
                         .$cell_id == "cell01"|
                           .$cell_id == "cell06"|
                           .$cell_id == "cell07"|
                           .$cell_id == "cell08"|
                           .$cell_id == "cell18"|
                           .$cell_id == "cell19"|
                           .$cell_id == "cell20"|
                           .$cell_id == "cell21"|
                           .$cell_id == "cell23"|
                           .$cell_id == "cell24"|
                           .$cell_id == "cell25"|
                           .$cell_id == "cell26"|
                           .$cell_id == "cell27"|
                           .$cell_id == "cell28",
                         F)) %>% 
  mutate(bl_activity = T) %>% 
  mutate(bl_activity = replace(bl_activity,
                               .$cell_id == "cell05"|
                                 .$cell_id == "cell10"|
                                 .$cell_id == "cell12"|
                                 .$cell_id == "cell13"|
                                 .$cell_id == "cell14"|
                                 .$cell_id == "cell15"|
                                 .$cell_id == "cell21"|
                                 .$cell_id == "cell24",
                               F))
                         
  



ggplot(data = BASE_TO_PLOT %>% filter(bl_activity == T, pinch == F, stim_cond == "baseline" | stim_cond == "d"),
       mapping = aes(x = forcats::fct_relevel(stim_cond, "baseline", "d"),
                     y = MFR)) +
  #theme(panel.background = element_rect(fill = 0)) +
  theme_minimal() +
  theme(axis.text = element_text(size = 20),
        text = element_text(size = 20)) +
  ylim(0,35)+
  
  geom_boxplot(width = 0.2, alpha = 0.5) +
  geom_point(shape = 21, 
             fill = "#EB8104", 
             #color = "white",
             size = 4) +
  #stroke = 2) +
  geom_line(aes(group = cell_id), color = "#EB8104") +
  geom_label_repel(mapping = aes(label = cell_id),
                   nudge_x = 0.15,
                   direction = "y",
                   hjust = -0.5,
                   segment.size = .9) 
