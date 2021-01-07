coords <- read.csv(file.path("turning", "coords.csv")) %>%
  as_tibble() %>%
  mutate(rec_time = (frame / 30) %>% round(digits = 2)) %>%
  unite(col = "head", head1, head2, sep = "_") %>%
  unite(col = "center", center1, center2, sep = "_") %>%
  unite(col = "tail", tail1, tail2, sep = "_") %>%
  gather(
    key = "body_part",
    value = "coordinate",
    head, center, tail
  ) %>%
  extract(col = coordinate, into = c("X", "Y"), regex = "([[:alnum:]]+)_([[:alnum:]]+)") %>%
  mutate(X = as.integer(X), Y = as.integer(Y)) %>%
  dplyr::group_by(body_part) %>%
  mutate(
    X_diff = diff(X) %>% c(NA),
    Y_diff = diff(Y) %>% c(NA)
  ) %>%
  mutate(d = sqrt(X_diff^2 + Y_diff^2)) %>%
  mutate(speed = d / frame) %>%
  ungroup()

coords[!duplicated(coords$body_part), ]

coords %>%
  dplyr::group_by(body_part) %>%
  slice(1, 2)

### to animate the plot these packages need to be installed
library(gganimate)
library(transformr)
library(gifski)

# coords %>%
#   dplyr::filter(sec%%1==0) %>%
p <- coords %>%
  # dplyr::group_by(body_part) %>%
  # slice(1:50) %>%

  ggplot(
    mapping = aes(x = X, y = Y)
  ) +
  geom_point(aes(shape = body_part), size = 2, alpha = .2) +
  geom_line(aes(group = frame), alpha = 0.2)
# scale_color_continuous(guide = F)
p

p +
  gganimate::transition_time(rec_time) + 
  labs(title = "Time: {frame_time %>% round(digits = 2)}s") +
  shadow_mark(alpha = 0.1, size = 0.5)

coords %>% dplyr::filter(body_part == "head") %>% pull(rec_time)

# geom_segment(
#   aes(group = frame
#     # xend = c(tail(X, n = -1), NA),
#     # yend = c(tail(Y, n = -1), NA)
#   )
# )
