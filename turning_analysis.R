### reading data, adding new variables ----

coords <- read.csv(file.path("turning", "coords.csv")) %>%
  as_tibble() %>%
  mutate(rec_time = (frame / 30) %>% round(digits = 2)) %>%
  mutate(stim = ifelse(
    frame >= 1159 & frame <= 1307 | frame >= 1798 & frame <= 2094,
    yes = T, no = F)) %>% 
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
  slice(1)


### plotting and animating data ----
# to animate the plot these packages need to be installed
library(gganimate)
library(transformr)
library(gifski)


# generating ggplot default colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

n = 2
cols = gg_color_hue(n)

p <- coords %>%
  ggplot(
    mapping = aes(x = X, y = Y)
  ) +
  geom_point(aes(shape = body_part), size = 2) +
  geom_line(aes(group = frame, col = stim)) +
  scale_color_manual(
    values = c("black" ,cols[2]),
    name = "Stimulus",
    labels = c("Off", "On")) +    
  scale_shape(
    name = "Body part",
    labels = c("Center", "Head", "Tail"))
p

anim <- p +
  gganimate::transition_time(rec_time) + 
  labs(title = "Time: {frame_time %>% round(digits = 2)}s") +
  shadow_mark(alpha = 0.1, size = 0.5)

gganimate::animate(anim, height = 800, width = 800)

### saving animation
gganimate::anim_save("mouse_trajec",animation = anim)


### calculating angles ----

coords %>% dplyr::filter(body_part == "head") %>% pull(rec_time)

# geom_segment(
#   aes(group = frame
#     # xend = c(tail(X, n = -1), NA),
#     # yend = c(tail(Y, n = -1), NA)
#   )
# )

### angles

coords %>%
  dplyr::group_by(body_part) %>%
  slice(1,2) %>%
  dplyr::filter(body_part != "tail") %>% 
  
  ggplot(
    mapping = aes(x = X, y = Y)
  ) +
  geom_point(aes(shape = body_part), size = 2) +
  geom_line(aes(group = frame, color = as.character(frame))) +
  geom_text_repel(aes(label = paste0(X,",",Y))) +
  geom_label_repel(
    data = coords %>%
      dplyr::group_by(body_part) %>%
      slice(1,2) %>%
      dplyr::filter(body_part != "tail") %>% 
      dplyr::group_by(frame) %>% 
      summarise(vect_X = diff(X), vect_Y = diff(Y)),
    aes(x = 750, y = 525, label = paste0(vect_X,",",vect_Y))
    )
  

coords %>%
  dplyr::group_by(body_part) %>%
  #slice(1:3) %>%
  dplyr::filter(body_part != "tail") %>% 
  dplyr::group_by(frame) %>% 
  summarise(vectx = diff(X) %>% abs(), vecty = diff(Y) %>% abs) %>%
  mutate(leadx = lead(vectx), leady = lead(vecty)) %>% 
  dplyr::group_by(frame) %>%
  mutate(
    # vectx,
    # vecty,
    # leadx,
    # leady,
    angle = matlib::angle(
      x = c(vectx, vecty), 
      y = c(leadx,leady), 
      degree = T) %>% as.vector()) 
  


library(matlib)
matA <- matrix(c(3, 1), nrow = 2)  ##column vectors
matB <- matrix(c(5, 5), nrow = 2)
angle(as.vector(matA), as.vector(matB)) 

coords %>%
  dplyr::group_by(body_part) %>%
  slice(1,2) %>%
  dplyr::filter(body_part != "tail") %>% 
  group_by(frame) %>% 
  summarise(as.vector(X), as.vector(Y))

center0 <- c(750,500)
head0 <- c(760,500)
frame0 <- c(center0,head0) 

center1 <- c(755,525)
head1 <- c(755,550)
frame1 <- c(center1,head1)

matlib::angle(frame0,frame1,degree = T)

vectors <- tibble(
  body = c("cent","h","cent","h"),
  frame = c(0,0,1,1),
  X = c(750,760,755,755),
  Y = c(500,500,525,550)
  )

vectors %>% 
  group_by(frame) %>% 
  summarise(dx=diff(X),dy=diff(Y)) %>% 
  dplyr::filter(frame=="1") %>% 
  select(-1) %>% as.matrix() %>% as.vector()

matlib::angle(
  x = vectors %>% 
    group_by(frame) %>% 
    summarise(dx=diff(X),dy=diff(Y)) %>% 
    dplyr::filter(frame=="0") %>% 
    select(-1) %>% as.matrix() %>% as.vector(),
  y = vectors %>% 
    group_by(frame) %>% 
    summarise(dx=diff(X),dy=diff(Y)) %>% 
    dplyr::filter(frame=="1") %>% 
    select(-1) %>% as.matrix() %>% as.vector()
  
    )

vectors %>%   
ggplot(mapping = aes(x = X, y = Y)) +
  geom_point(aes(col=body)) + 
  geom_line(aes(group = frame))

