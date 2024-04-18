# Run this only once to load raw data, clean it from junk and save it again to reduce file size for GitHub. 
# Original raw data is not version controlled. Cleaned data is loaded for further analysis


stims_tmp <- read.csv(file.path("data", "raw_data", "PnO_stim_rotation", "rotation_stims.csv")) %>%
  as_tibble() %>%
  ### filter junk data
  dplyr::filter(
    stringr::str_detect(animal, pattern = "b")
  ) %>%
  write.csv(file.path("data", "PnO_stim_rotation", "rotation_stims_cleaned.csv"), row.names = F)


coords_raw_tmp <- read.csv(file.path("data", "raw_data", "PnO_stim_rotation", "rotation_all.csv")) %>%
  as_tibble() %>%
  ### filter junk data
  dplyr::filter(
    stringr::str_detect(animal, pattern = "b")
  ) %>%
  write.csv(file.path("data", "PnO_stim_rotation", "rotation_all_cleaned.csv"), row.names = F)