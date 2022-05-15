library(tidyverse)
RawDataFilter <- function(df) {
  #browser()
  df %>%
    dplyr::filter(
      stringr::str_detect(animal, pattern = "b")
      
    ) %>% 
    as_tibble() %>% return()
}

stims <- read.csv(file.path("data", "PnO_stim_rotation", "rotation_stims.csv")) %>%
  as_tibble() %>% 
  RawDataFilter() %>% 
  dplyr::mutate(
    intensity = paste0(intensity, " mW"),
  ) 
  
  
  



coords <- read.csv(file.path("data", "PnO_stim_rotation", "rotation_all.csv")) %>%
  as_tibble() %>%
  dplyr::filter(
    stringr::str_detect(animal, pattern = "b")
  ) %>% 
  dplyr::mutate(
    animal = gsub(pattern = ".*rotation", replacement = "", animal) %>%
      str_remove_all(pattern = "\\\\") %>% str_remove(pattern = "fix") %>% str_remove(pattern = "control")
  )


identical(
stims$animal %>% unique() %>% sort(),
coords$animal %>% unique() %>% sort()
)  

### sessions
stims %>% 
  group_by(animal) %>% 
  summarise(unique(session) %>% length())

stims %>% 
  dplyr::mutate(
    stim_length = end - start,
    s2 = lag(session,default = "1"),
    s3 = ifelse(s2 == session, 0,1),
  ) %>% 
  group_by(animal) %>%
  mutate(session = cumsum(s3)) %>% 
  summarise(unique(session) %>% length())



coords %>% 
  group_by(animal) %>% 
  summarise(unique(session) %>% length())
  
coords %>% 
  dplyr::mutate(
    s2 = lag(session,default = "1"),
    s3 = ifelse(s2 == session, 0,1),
  ) %>% 
  group_by(animal) %>%
  mutate(session = cumsum(s3)) %>% 
  summarise(unique(session) %>% length())


###intensities

stims %>% 
  group_by(animal) %>% 
  summarise(unique(intensity) %>% length())

stims %>% 
  group_by(animal, session) %>% 
  summarise(unique(intensity))

stims %>% 
  group_by(animal, session) %>% 
  summarise(intensity) %>% group_by(animal, session) %>% count(intensity) %>% View()





coords %>% 
  dplyr::mutate(
    s2 = lag(session,default = "1"),
    s3 = ifelse(s2 == session, 0,1),
  ) %>% 
  group_by(animal) %>%
  mutate(session = cumsum(s3)) %>% 
  group_by(animal) %>%
  mutate(session = cumsum(s3)) %>%
  select(-s2, -s3) %>% 
  ungroup() %>% 
  
  full_join(
    stims %>% 
      dplyr::mutate(
        stim_length = end - start,
        s2 = lag(session,default = "1"),
        s3 = ifelse(s2 == session, 0,1),
      ) %>% 
      group_by(animal) %>%
      mutate(session = cumsum(s3)) %>% 
      select(-stim_duration, -s2, -s3) %>%
      rename("stim_number" = stim.number) %>%
      gather(key = "timing", value = "frame", start, end)
  ) 
  # %>% 
  # group_by(animal) %>% 
  # summarise(unique(intensity) %>% na.omit() %>% length())
  