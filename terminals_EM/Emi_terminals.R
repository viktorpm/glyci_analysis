wd_path <- file.path('f:',
                     '_R_WD',
                     '_data_frames'
)



library(R.matlab)
library(STAR)
library(plyr)
library(ggplot2)
library(see) #for half violinplot
library(gghalves) #for half violinplot
library(beanplot)
library(tidyverse)

### cunstructing the path to the file and saving it as a character variable
file_path <- file.path('data',
                       'Emi_terminals',
                       'terminal_measures.csv' )  



TERMINALS <- read.table(file = file_path, sep = ';', header = T, na.string = 'na')
head(TERMINALS)




ggplot(data = TERMINALS, mapping = aes(x = postion, y = b_area)) +
  geom_half_violin(position = position_nudge(x=-.15),
                  trim = F,
                  side = "l",
                  scale = "area",
                  draw_quantiles = c(0.5),
                  aes(fill = postion)) +
  geom_jitter(position = position_jitter(0.1))



  
ggplot(data = TERMINALS,
       mapping = aes(x = origin, y = b_area)) +
  geom_half_violin(position = position_nudge(x=-.15),
                   trim = F,
                   side = "l",
                   scale = "count",
                   draw_quantiles = c(0.5)
                   #fill = alpha(0.1),
                   #aes(col = postion)
                   ) +
  geom_jitter(position = position_jitter(0.1),
              aes(col = postion)) +
  annotate("text",
           x = 0.7, 
           y = 1, 
           label = paste("n = ", 
                         TERMINALS %>%
                           group_by(origin) %>% 
                           summarise(n = length(b_area)) %>%
                           dplyr::filter(origin == "m2") %>% 
                           select(n) %>%
                           as.character())
           ) +
  annotate("text",
           x = 1.6, 
           y = 1, 
           label = paste("n = ", 
                         TERMINALS %>%
                           group_by(origin) %>% 
                           summarise(n = length(b_area)) %>%
                           dplyr::filter(origin == "unknown") %>% 
                           select(n) %>%
                           as.character())
  )



########################################################

with(TERMINALS, hist(b_area[postion == 'non_labeled']))
with(TERMINALS, 
     plot(density(b_area[postion == 'non_labeled'])
          )
     )


with(TERMINALS, hist(b_area[postion == 'glyt2']))
with(TERMINALS, 
     plot(density(b_area[postion == 'glyt2'])
     )
)


with(TERMINALS, hist(b_area[postion == 'glyt2_cell']))
with(TERMINALS, 
     plot(density(b_area[postion == 'glyt2_cell'])
     )
)




with(TERMINALS, hist(b_area))
with(TERMINALS, 
     plot(density(b_area)
     )
)


#################### terminal area distribution ###########################


### bouton area on glyt2 cell
with(TERMINALS, 
     plot(density(b_area[postion == 'glyt2_cell'])
          ,ylim = c(0,1.5)
          ,xaxs = 'i'
          ,yaxs = 'i'
          #,axes = F
          #,xlim = c(0,2)
          )
) 


### bouton area on proximal glycinergic dendrite
with(TERMINALS, 
     points(density(b_area[postion == 'glyt2_prox'])
            ,type = 'l'
            ,col = 2
          )
)




### bouton areas combined (plot #1)


# with(TERMINALS, 
#      plot(density(b_area[origin == 'unknown'])
#           ,ylim = c(0,1.8)
#           ,xaxs = 'i'
#           ,yaxs = 'i'          
#           )
# ) 


with(TERMINALS, 
     hist(b_area[origin == 'unknown']
          #,ylim = c(0,1.8)
          ,xaxs = 'r'
          ,yaxs = 'r'
          ,breaks = 25
          ,col = rgb(0, 146, 146, maxColorValue = 255, 255)
          ,lty = 'blank'
          )
     )

abline(v = with(TERMINALS,
                mean(b_area[origin == 'unknown'])
                )
       ,col = rgb(0, 146, 146, maxColorValue = 255, 255)
       )



### bouton area with m2 origin and glyt2 target (plot #2)


# with(TERMINALS, 
#      points(density(b_area[postion == 'glyt2'])
#             ,type = 'l'
#             ,col = 2
#           )
# )

with(TERMINALS, 
     hist(b_area[postion == 'glyt2']
          #,ylim = c(0,1.8)
          ,xaxs = 'r'
          ,yaxs = 'r'
          ,breaks = 15
          ,add = T
          ,col = rgb(182, 109, 255, maxColorValue = 255, 100)
          ,lty = 'blank'
     )
)

abline(v = with(TERMINALS,
                mean(b_area[postion == 'glyt2'])
                )
       ,col = rgb(182, 109, 255, maxColorValue = 255, 255)
       )


### bouton area with m2 origin targeting unidentified dendrite (plot #3)


# with(TERMINALS, 
#      points(density(b_area[postion == 'non_labeled' & origin == 'm2'])
#             ,type = 'l'
#             ,col = 3
#             )
# ) 


with(TERMINALS, 
     hist(b_area[postion == 'non_labeled' & origin == 'm2']
          #,ylim = c(0,1.8)
          ,xaxs = 'r'
          ,yaxs = 'r'
          ,breaks = 15
          ,add = T
          ,col = rgb(219, 209, 0, maxColorValue = 255, 100)
          ,lty = 'blank'
          )
     )

abline(v = with(TERMINALS,
                mean(b_area[postion == 'non_labeled' & origin == 'm2'])
                )
       ,col = rgb(219, 209, 0, maxColorValue = 255, 255)
       )


stripchart(TERMINALS$b_area[TERMINALS$origin == 'unknown']
           ,pch = 19
           ,vertical = T
           ,method = 'jitter'
           ,xlim = c(0,4)
           ,col = rgb(0, 146, 146, maxColorValue = 255, 255)
           ,ylim = c(-.5,3)
           )

beanplot(TERMINALS$b_area[TERMINALS$origin == 'unknown']
        ,add = T
        ,at = .8
        #,boxwex = .3
        ,col = rgb(0, 146, 146, maxColorValue = 255, 255)
        ,side = 'first'
        ,what = c(0,1,1,1)
        ,overallline = 'mean'
        #,kernel = 'rectangular'
        ,ll = 0.2
        )




stripchart(TERMINALS$b_area[TERMINALS$postion == 'glyt2']
           ,pch = 19
           ,vertical = T
           ,method = 'jitter'
           ,add = T
           ,col = rgb(182, 109, 255, maxColorValue = 255, 255)
           ,at = 2
           )

beanplot(TERMINALS$b_area[TERMINALS$postion == 'glyt2']
        ,add = T
        ,at = 1.8
        #,boxwex = .3
        ,col = rgb(182, 109, 255, maxColorValue = 255, 255)
        ,side = 'first'
        ,what = c(0,1,1,0)
        ,overallline = 'mean'
        #,kernel = 'rectangular'
)


stripchart(TERMINALS$b_area[TERMINALS$postion == 'non_labeled' & TERMINALS$origin == 'm2']
           ,pch = 19
           ,vertical = T
           ,method = 'jitter'
           ,add = T
           ,col = rgb(219, 209, 0, maxColorValue = 255, 255)
           ,at = 3
           )

beanplot(TERMINALS$b_area[TERMINALS$postion == 'non_labeled' & TERMINALS$origin == 'm2']
        ,add = T
        ,at = 2.8
        #,boxwex = .3
        ,col = rgb(219, 209, 0, maxColorValue = 255, 255)
        ,side = 'first'
        ,what = c(0,1,1,0)
        ,overallline = 'mean'
        #,kernel = 'rectangular'
        )



####################### dendrite diameter distribution ##########################


par(mfrow = c(2,2))
  with(TERMINALS, 
       plot(density(d_dendrite[postion == 'glyt2'])
            #,ylim = c(0,1.5)
            )
  )
  
  with(TERMINALS, 
       plot(density(d_dendrite[postion == 'non_labeled']
                    ,na.rm = T
                    )
            #,ylim = c(0,1.5)
            )
  )
  
  
  with(TERMINALS, hist(d_dendrite[postion == 'glyt2'] ,breaks = 10))
  with(TERMINALS, hist(d_dendrite[postion == 'non_labeled'] ,breaks = 15))
par(mfrow = c(1,1))




m2_origin_target <- hist(TERMINALS$d_dendrite[TERMINALS$origin == 'm2'] ,breaks = 15)
hist(TERMINALS$d_dendrite[TERMINALS$origin == 'm2'] ,breaks = 15)
max(m2_origin_target$counts)

m2_origin_target_density <- density(TERMINALS$d_dendrite[TERMINALS$origin == 'm2'], na.rm = T)
plot(density(TERMINALS$d_dendrite[TERMINALS$origin == 'm2'], na.rm = T))
abline(v = m2_origin_target_density$x[135], col = 2)
abline(v = m2_origin_target_density$x[245], col = 4)




# stripchart(dendrites$d_dendrite ~ droplevels(dendrites$postion),
#            vertical = T,
#            pch = 19,
#            xlim = c(0,3),
#            ylim = c(0,2),
#            ylab = "Dendrite diameter")
# 
# beanplot(TERMINALS$d_dendrite[TERMINALS$postion == 'non_labeled']
#          ,add = T
#          ,at = 1.8
#          ,cut = 1
#          ,boxwex = 1.2
#          ,col = rgb(219, 209, 0, maxColorValue = 255, 255)
#          ,side = 'first'
#          ,what = c(0,1,1,0)
#          #,overallline = 'median'
#          #,kernel = 'rectangular'
# )
# 
# beanplot(TERMINALS$d_dendrite[TERMINALS$origin == 'm2']
#          ,add = T
#          ,at = 0.8
#          ,cut = 1
#          ,boxwex = 1.2
#          ,col = rgb(219, 209, 0, maxColorValue = 255, 255)
#          ,side = 'first'
#          ,what = c(0,1,1,0)
#          #,overallline = 'median'
#          #,kernel = 'rectangular'
# )













#which(m2_origin_target_density$y == max(m2_origin_target_density$y))
which.max(m2_origin_target_density$y)

local_max_2 <- max(m2_origin_target_density$y[m2_origin_target_density$x > 0.5])
which(m2_origin_target_density$y == local_max_2)







length(TERMINALS$d_dendrite[TERMINALS$origin == 'm2'])

levels(TERMINALS$origin)


dendrites <- subset(TERMINALS, TERMINALS$postion == 'glyt2' | TERMINALS$postion == 'non_labeled')
levels(dendrites$postion)
length(levels(droplevels(dendrites$postion)))


with(dendrites,
     stripchart(d_dendrite ~ droplevels(postion)
                ,vertical = T
                ,pch = 19
                ,xlim = c(-1,4)
                )
     )
abline(h = m2_origin_target_density$x[245], col = 2)
abline(h = m2_origin_target_density$x[135], col = 4)



with(dendrites,
     boxplot(d_dendrite ~ droplevels(postion),
             add = T,
             TERMINALS$b_area[TERMINALS$origin == 'unknown'] = 0.1,
             col = rgb(0,0,0, maxColorValue = 255, alpha=0),
             at = c(1.2, 2.2)
             )
     )



with(TERMINALS, length(d_dendrite[postion == 'non_labeled']))
with(TERMINALS, length(d_dendrite[postion == 'glyt2']))


