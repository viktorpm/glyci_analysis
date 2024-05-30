# Load necessary libraries
library(tidyverse)

# Define custom operators and functions
`%!in%` <- Negate(`%in%`)
sem <- function(x) sd(x) / sqrt(length(x)) # Standard error of the mean function

# Check and install required packages
required_packages <- c("imputeTS")
missing_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}
# Load packages
library(imputeTS)

# Find and load all necessary libraries
library_paths <- list.files(pattern = ".Rmd|.R", full.names = TRUE) %>%
  map(~ readLines(con = .x) %>%
        grep(pattern = "library\\(", x = ., value = TRUE)) %>%
  unlist() %>%
  str_extract(pattern = "(?<=library\\()[^)]+") %>%
  unique()

library_paths <- library_paths[library_paths != "STARS"]

lapply(library_paths, library, character.only = TRUE)


# Define global theme for ggplot
global_facet_theme <- theme(
  strip.text = element_text(size = 15),
  axis.text = element_text(size = 12),
  axis.title = element_text(size = 15)
)

# Custom function to remove junk animals from the data
RawDataFilter <- function(df) {
  df %>%
    filter(str_detect(animal, pattern = "b")) %>%
    as_tibble()
}

# Custom function to recalculate distance and speed after removing frames 
RecalculateDistance <- function(df) {
  df %>%
    dplyr::mutate(
      X_diff = c(diff(X_interp), NA) %>% round(4),
      Y_diff = c(diff(Y_interp), NA) %>% round(4),
      X_diff = ifelse(
        test = X_diff < quantile(X_diff, na.rm = T)["25%"] - 3 * IQR(X_diff, na.rm = T) |
          X_diff > quantile(X_diff, na.rm = T)["75%"] + 3 * IQR(X_diff, na.rm = T),
        yes = 0,
        no = X_diff
      ),
      Y_diff = ifelse(
        Y_diff < quantile(Y_diff, na.rm = T)["25%"] - 3 * IQR(Y_diff, na.rm = T) |
          Y_diff > quantile(Y_diff, na.rm = T)["75%"] + 3 * IQR(Y_diff, na.rm = T),
        yes = 0,
        no = Y_diff
      ),
      d = sqrt(X_diff^2 + Y_diff^2),
      d_mm = d / scale,
      speed_instant = d_mm / time_res,
      mmps = speed_instant * (1 / time_res),
      mps = (speed_instant * (1 / time_res)) / 1000
    )
}