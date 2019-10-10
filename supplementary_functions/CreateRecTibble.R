CreateRecTibble <- function(AP_times, stim_times) {
  recordings <- bind_rows(AP_times, stim_times)

  recordings <- recordings %>%
    mutate(stim_freq = 0)


  reciprocal <- function(x) {
    y <- 1 / x
    return(y)
  }

  remove_zero <- function(x) {
    x[which(x < 1)] <- x[which(x < 1) - 1]
  }


  recordings$stim_freq <-
    recordings %>%
    # dplyr::filter(signal_type == 'stim') %>%
    dplyr::select(signal_time) %>%
    as.matrix() %>%
    diff() %>%
    append(2) %>%
    reciprocal() %>%
    round(2) ### to differentiate test pulses from last elements of the stim train

  ### last elements of the stim trains replaced by the previous element
  ### only test pulses' freqs are smaller than 1
  recordings$stim_freq[which(recordings$stim_freq < 1)] <- recordings$stim_freq[which(recordings$stim_freq < 1) - 1] 
  
  ### elemnts smaller than 1 == test pulse
  ### round the rest as numeric (Warning: NAs introduced by coercion)
  ### replaces NAs with "test pulse"
  recordings$stim_freq[which(recordings$stim_freq < 1)] <- "test pulse"
  recordings$stim_freq <- round(as.numeric(recordings$stim_freq))
  recordings$stim_freq[is.na(recordings$stim_freq) == T] <- "test pulse"
  
  
  ### stim freq at AP's has no value
  recordings$stim_freq[recordings$signal_type == "AP"] <- "no value"
 

  recordings <- recordings %>%
    mutate(available_freqs = NA)




  stim_freqs_by_file <- recordings %>%
    dplyr::filter(signal_type == "stim") %>%
    dplyr::group_by(file_name) %>%
    distinct(stim_freq)

  ### I have NO IDEA of what happens here but it works
  ### https://stackoverflow.com/questions/40093595/dplyr-group-by-and-convert-multiple-columns-to-a-vector
  library(purrrlyr)
  stim_freqs_by_file <- stim_freqs_by_file %>%
    slice_rows("file_name") %>%
    by_slice(function(x) unlist(x), .to = "freqs_in_file")

  ### to filter by available freqs in the selected file

  for (i in 1:length(stim_freqs_by_file$file_name))
  {
    name_string_test <- recordings$file_name %>%
      str_detect(pattern = stim_freqs_by_file$file_name[i])
    recordings$available_freqs[which(name_string_test == T)] <- paste(
      stim_freqs_by_file$freqs_in_file[[i]], collapse = ", "
    )
  }


  return(recordings)
}