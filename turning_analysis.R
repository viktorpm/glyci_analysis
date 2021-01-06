coords <- read.csv(file.path("turning", "coords.csv")) %>%
  as_tibble() %>%
  mutate(sec = frame / 30) %>% 
  mutate(diff_head1 = diff(head1) %>% c(NA),
         diff_head2 = diff(head2) %>% c(NA)) %>% 
  mutate(d = sqrt(diff_head1^2+diff_head2^2)) %>% 
  mutate(speed = d/frame)
  
  





coords %>%
  #dplyr::filter(sec%%1==0) %>% 
  unite(col = "head", head1, head2, sep = "_") %>% 
  unite(col = "center", center1, center2, sep = "_") %>% 
  unite(col = "tail", tail1, tail2, sep = "_") %>% 
  gather(key = "body_part", 
         value = "coordinate", 
         head, center, tail) %>% 
  extract(col = coordinate, into = c("X", "Y"), regex = "([[:alnum:]]+)_([[:alnum:]]+)") %>% 
  mutate(X = as.integer(X), Y = as.integer(Y)) %>% 
  dplyr::filter(body_part == "head") %>% 

  ggplot(
    mapping = aes(x = X, y = Y)
  ) +
  geom_point() 
  geom_segment(
    aes(
      xend = c(tail(body_part, n = -1), NA),
      yend = c(tail(body_part, n = -1), NA)
    )
  )
?tail
