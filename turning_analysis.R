### reading data, adding new variables ----

coords <- read.csv(file.path("turning", "coords.csv")) %>%
  as_tibble() %>%
  mutate(rec_time = (frame / 30) %>% round(digits = 2)) %>%
  mutate(stim = ifelse(
    frame >= 1159 & frame <= 1307 | frame >= 1798 & frame <= 2094,
    yes = T, no = F
  )) %>%
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
  mutate(speed = d / frame,
         sum_d = ifelse(rec_time %% 1 == 0, cumsum(d), NA),
         cumsumd = cumsum(d),
         sec = cumsum (ifelse(rec_time %% 1 == 0, 1, 0))
         ) %>% 
  group_by(sec, body_part) %>% 
  mutate(sum_d2 = sum(d))
  # ungroup() %>% 
  # mutate(sum_d = NA) %>% 
  # dplyr::group_by(body_part) %>% 
  # dplyr::mutate(sum_d = replace(sum_d, list = (rec_time %% 1 == 0), values = rowSums(d, 30)))
  


library(RcppRoll)
 
coords[!duplicated(coords$body_part), ]

coords %>%
  dplyr::group_by(body_part) %>%
  slice(1)


### plotting and animating data ----
# to animate the plot these packages need to be installed
library(gganimate)
library(transformr)
library(gifski)


# generating ggplot default colors
gg_color_hue <- function(n) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

n <- 2
cols <- gg_color_hue(n)

p <- coords %>%
  ggplot(
    mapping = aes(x = X, y = Y)
  ) +
  geom_point(aes(shape = body_part), size = 2) +
  geom_line(aes(group = frame, col = stim)) +
  scale_color_manual(
    values = c("black", cols[2]),
    name = "Stimulus",
    labels = c("Off", "On")
  ) +
  scale_shape(
    name = "Body part",
    labels = c("Center", "Head", "Tail")
  )
p

anim <- p +
  gganimate::transition_time(rec_time) +
  labs(title = "Time: {frame_time %>% round(digits = 2)}s") +
  shadow_mark(alpha = 0.1, size = 0.5)


### saving animation
gganimate::anim_save(
  "mouse_trajec",
  animation = gganimate::animate(anim, height = 800, width = 800, fps = 5)
)


### calculating angles ----


### angles

coords %>%
  dplyr::group_by(body_part) %>%
  slice(1:6) %>%
  dplyr::filter(body_part != "tail") %>%
  ggplot(
    mapping = aes(x = X, y = Y)
  ) +
  geom_point(aes(shape = body_part), size = 2) +
  geom_line(aes(group = frame, color = as.character(frame))) +
  geom_text_repel(aes(label = paste0(X, ",", Y))) +
  geom_label_repel(
    data = coords %>%
      dplyr::group_by(body_part) %>%
      slice(1:6) %>%
      dplyr::filter(body_part != "tail") %>%
      dplyr::group_by(frame) %>%
      summarise(vect_X = diff(X), vect_Y = diff(Y)),
    aes(x = 750, y = 525, label = paste0(vect_X, ",", vect_Y), col = as.character(frame))
  )


coords %>%
  dplyr::filter(body_part != "tail") %>%
  #dplyr::group_by(body_part) %>%
  #slice(1:3) %>%
  ### each frame is present twice (for each body part)
  dplyr::group_by(frame) %>%
  summarise(
    # vectxabs = diff(X) %>% abs(),
    # vectyabs = diff(Y) %>% abs(),
    vectx = diff(X),
    vecty = diff(Y)
  ) %>%
  mutate(
    # leadxabs = lead(vectxabs),
    # leadyabs = lead(vectyabs),
    leadx = lead(vectx),
    leady = lead(vecty)
  ) %>%
  dplyr::group_by(frame) %>%
  mutate(
    angle = matlib::angle(
      x = c(vectx, vecty),
      y = c(leadx, leady),
      degree = T
    ) %>% as.vector(),
    signed_angle = (atan2(leadx, leady) - atan2(vectx, vecty)) * 180 / (pi)   
    # angleabs = matlib::angle(
    #   x = c(vectxabs, vectyabs),
    #   y = c(leadxabs, leadyabs),
    #   degree = T) %>% as.vector()
  ) %>% 
  add_column(
    rec_time = coords %>% dplyr::filter(body_part == "head") %>% pull(rec_time),
    stim = coords %>% dplyr::filter(body_part == "head") %>% pull(stim)
  ) %>% 
  ungroup() %>% 
  dplyr::group_by(stim) %>% 
  mutate(sum_angle = cumsum(signed_angle)) 
  dplyr::group_by(stim) %>% 
  summarise(sum(signed_angle, na.rm = T))



atan2()

library(matlib)
matA <- matrix(c(3, 1), nrow = 2) ## column vectors
matB <- matrix(c(5, 5), nrow = 2)
angle(as.vector(matA), as.vector(matB))

coords %>%
  dplyr::group_by(body_part) %>%
  slice(1, 2) %>%
  dplyr::filter(body_part != "tail") %>%
  group_by(frame) %>%
  summarise(as.vector(X), as.vector(Y))

center0 <- c(750, 500)
head0 <- c(760, 500)
frame0 <- c(center0, head0)

center1 <- c(755, 525)
head1 <- c(755, 550)
frame1 <- c(center1, head1)

matlib::angle(frame0, frame1, degree = T)

vectors <- tibble(
  body = c("cent", "h", "cent", "h"),
  frame = c(0, 0, 1, 1),
  X = c(750, 760, 755, 755),
  Y = c(500, 500, 525, 550)
)

vectors %>%
  group_by(frame) %>%
  summarise(dx = diff(X), dy = diff(Y)) %>%
  dplyr::filter(frame == "1") %>%
  select(-1) %>%
  as.matrix() %>%
  as.vector()

matlib::angle(
  x = vectors %>%
    group_by(frame) %>%
    summarise(dx = diff(X), dy = diff(Y)) %>%
    dplyr::filter(frame == "0") %>%
    select(-1) %>% as.matrix() %>% as.vector(),
  y = vectors %>%
    group_by(frame) %>%
    summarise(dx = diff(X), dy = diff(Y)) %>%
    dplyr::filter(frame == "1") %>%
    select(-1) %>% as.matrix() %>% as.vector()
)

vectors %>%
  ggplot(mapping = aes(x = X, y = Y)) +
  geom_point(aes(col = body)) +
  geom_line(aes(group = frame))
