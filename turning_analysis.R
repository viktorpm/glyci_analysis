### reading data, adding new variables ----
scale <- 1.38 ### pixel/mm
fps <- 30 ### frame rate
time_res <- 1 / fps ### time between frames

stims <- read.csv(file.path("data", "turning", "GIIFM19",  "GIIFM19Frames_20201116.csv")) %>% 
  as_tibble() %>%
  mutate(session = str_remove(session, "_")) %>% 
  gather(key = "timing", value = "frame", start, end)


coords <- read.csv(file.path("data", "turning", "GIIFM19.csv")) %>%
  as_tibble() %>%
  mutate(session = str_remove(session, pattern = ".dat")) %>%
  full_join(stims) %>% 
  dplyr::group_by(session) %>% 
  group_by(grp = cumsum(!is.na(stim.number))) %>% 
  mutate(stim = F, 
         stim = replace(stim, first(timing) == 'start', T)) %>% 
  ungroup() %>% 
  select(-grp, -timing, -animal) %>% 
  group_by(session) %>% 
  dplyr::mutate(rec_time = (frame / fps) %>% round(digits = 2)) %>%
  unite(col = "head", head1, head2, sep = "_") %>%
  unite(col = "center", center1, center2, sep = "_") %>%
  unite(col = "tail", tail1, tail2, sep = "_") %>%
  gather(
    key = "body_part",
    value = "coordinate",
    head, center, tail
  ) %>%
  tidyr::extract(col = coordinate, into = c("X", "Y"), regex = "([[:alnum:]]+)_([[:alnum:]]+)") %>%
  dplyr::mutate(
    X = ifelse(
      as.integer(X) == 0 & as.integer(Y) == 0, 
      yes = NA,
      no = as.integer(X)),
    Y = ifelse(
      as.integer(X) == 0 & as.integer(Y) == 0, 
      yes = NA,
      no = as.integer(Y))) %>%
  dplyr::group_by(body_part) %>%
  mutate(
    X_diff = c(diff(X), NA),
    Y_diff = c(diff(Y), NA)
  ) %>%##
  # dplyr::mutate(
  #   X_diff = ifelse(abs(X_diff) > 30, NA, X_diff),
  #   Y_diff = ifelse(abs(Y_diff) > 30, NA, Y_diff)
  # ) %>%
  dplyr::mutate(
    d = sqrt(X_diff^2 + Y_diff^2),
    d_mm = d / scale
  ) %>%
  dplyr::mutate(
    speed_instant = d_mm / time_res,
    mmps = speed_instant * (1 / time_res),
    mps = (speed_instant * (1 / time_res)) / 1000,
    # sum_d_mm = ifelse(rec_time %% 1 == 0, cumsum(d_mm), NA),
    # cumsumd = cumsum(d_mm),
    sec = cumsum(ifelse(rec_time %% 1 == 0, 1, 0))
  ) %>%
  group_by(sec, body_part) %>%
  dplyr::mutate(speed_sec = sum(d_mm, na.rm = T)) %>%
  ungroup() 
  
  # group_by(stim) %>% 
  # dplyr::mutate(mean_mps_roll = zoo::rollmean(x = mps, k = 60,fill = NA))
# ungroup() %>%
# mutate(sum_d = NA) %>%
# dplyr::group_by(body_part) %>%
# dplyr::mutate(sum_d = replace(sum_d, list = (rec_time %% 1 == 0), values = rowSums(d, 30)))



# generating ggplot default colors
gg_color_hue <- function(n) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

n <- 2
cols <- gg_color_hue(n)

coords %>%
  dplyr::filter(session == "video0coords") %>% 
  ggplot(
    mapping = aes(x = X, y = Y, label = frame)
  ) +
  geom_point(aes(shape = body_part), size = 2) +
  geom_line(aes(group = frame, col = stim_number)) +
  scale_color_manual(
    values = c(cols[1], cols[2], cols[3], cols[4], cols[5], "black"),
    name = "Stimulus"
  ) +
  scale_shape(
    name = "Body part",
    labels = c("Center", "Head", "Tail")
  ) +
  xlim(c(0,1280)) +
  ylim(c(0,720)) +
  #facet_wrap(~ session) +
  gganimate::transition_time(rec_time,range = c(0,20))




coords %>%
  dplyr::filter(body_part == "center") %>%
  ggplot(mapping = aes(x = X_diff)) +
  geom_histogram()

coords %>%
  dplyr::filter(body_part == "head") %>% 
  ggplot() +
  geom_line(mapping = aes(x = rec_time, y = d_mm, color = stim, group = 1)) + 
  ylim(c(0,100)) + 
  facet_wrap(~ session)






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



p <- coords %>%
  dplyr::filter(session == "video4coords") %>% 
  ggplot(
    mapping = aes(x = X_interp, y = Y_interp)
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
  ) + 
  xlim(c(0,1280)) +
  ylim(c(0,720)) 
#  facet_wrap(~session)
p

p +
  gganimate::transition_time(rec_time)
  
  
anim <- p +
  gganimate::transition_time(rec_time) +
  labs(title = "Time: {frame_time %>% round(digits = 2)}s") +
  shadow_mark(alpha = 0.1, size = 0.5)


### saving animation
gganimate::anim_save(
  "mouse_trajec_4.gif",
  animation = gganimate::animate(anim, height = 800, width = 800, fps = 5)
)


### calculating angles ----


### angles

coords %>%
  dplyr::group_by(body_part) %>%
  slice(1:3) %>%
  dplyr::filter(body_part != "tail") %>%
  ggplot(
    mapping = aes(x = X, y = Y)
  ) +
  geom_point(aes(shape = body_part), size = 2) +
  geom_line(aes(group = frame, color = as.character(frame))) +
  geom_text_repel(aes(label = paste0(X, ",", Y),col = as.character(frame))) 
  # geom_label_repel(
  #   data = coords %>%
  #     dplyr::group_by(body_part) %>%
  #     slice(1:3) %>%
  #     dplyr::filter(body_part != "tail") %>%
  #     dplyr::group_by(frame) %>%
  #     summarise(vect_X = diff(X), vect_Y = diff(Y)),
  #   aes(x = 750, y = 525, label = paste0(vect_X, ",", vect_Y), col = as.character(frame))
  # )


coords %>%
  dplyr::filter(body_part != "tail") %>%
  # dplyr::group_by(body_part) %>%
  # slice(1:3) %>%
  ### each frame is present twice (for each body part)
  dplyr::group_by(session,frame) %>%
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
  # dplyr::group_by(frame) %>%
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
  dplyr::group_by(session, stim) %>%
  mutate(sum_angle = cumsum(signed_angle)) %>% View() 
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
  geom_point(aes(shape = body)) +
  geom_line(aes(group = frame, color = frame %>% as.character())) +
  geom_text_repel(aes(label = paste0(X, ",", Y),col = as.character(frame)))

