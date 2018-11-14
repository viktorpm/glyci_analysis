### recording site: IL 
### stimulation: glycinergic fibers in the IL around the recorded cells
### to visualize firing rate change before, during and after the stimulation @30Hz
### firing rate: no. of APs/time (length of the stimulus)
### control firing rates (before and after the stimulus): no. of APs in a time window with the same length as the stimulus

########################## Reading and organizing data ##########################

### function: separated by CapitalLetters
### data from external source: separated.by '.'
### variable: separated_by '_'


### notebook wd path (RUN ON NOTEBOOK):

wd_path <- file.path('c:',
                     'Users',
                     'Viktor',
                     'Documents',
                     'R_WD'
)

setwd(wd_path_NB)


file_path <- file.path('c:',
                       'Users',
                       'Viktor',
                       'OneDrive - MTA KOKI',
                       'professional',
                       '_R_WD',
                       '_data_frames',
                       'IL_juxta_GlyT2_fiber_stim'
                       )


### lab PC wd path (RUN ON lab PC):

wd_path <- file.path('f:',
                     '_R_WD'
)

setwd(wd_path)


file_path <- file.path('f:',
                       '_R_WD',
                       '_data_frames',
                       'folder',
                       '*.mat'
)


### loading packages
library(R.matlab)


### reading data from file
raw.rec <- readMat(paste0(file_path,
                          '/cell21_13_right_IL_t3_3030_stim_control.mat'
                          )
                   )



names(raw.rec)

names(raw.rec) = c('unit','level','stim','eeg')

level.stim = raw.rec$level[[5]]

end = length(level.stim)

level.start = raw.rec$level[[5]][seq(1,end,2),1]
level.stop = raw.rec$level[[5]][seq(1,end,2)+1,1]
level.stimLength = level.stop[1] - level.start[1] #calculates the length of the stim


AP.times = raw.rec$unit[[12]]

# # #APs during the first stim
# length(AP.times[level.start[1]<AP.times & AP.times<level.stop[1]])
# 
# # #APs 33.3s before first stim (less than 33s, first 13.7s)
# length(AP.times[level.start[1]-level.stimLength < AP.times & AP.times<level.start[1]])
# 
# 
# # #APs 33.3s after first stim
# length(AP.times[AP.times>level.stop[1] & AP.times<level.stop[1] + level.stimLength])


# number of APs before during and after stim

AP.numbers = matrix(0, nrow = length(level.start), ncol = 3)

for (i in 1:length(level.start)) {
  AP.numbers[i,1] = length(AP.times[level.start[i] - level.stimLength < AP.times & AP.times < level.start[i]])
  AP.numbers[i,2] = length(AP.times[level.start[i] < AP.times & AP.times < level.stop[i]])
  AP.numbers[i,3] = length(AP.times[AP.times > level.stop[i] & AP.times < level.stop[i] + level.stimLength])
}

# AP.numbers = AP.numbers[-6:-7,]

rownames(AP.numbers) = c('1st','2nd')#,'3d','4th','5th','6th','7th','8th','9th')#,'10th')
colnames(AP.numbers) = c('before','during','after')

# MFR = AP.numbers/33.3

cell24 = AP.numbers


# cellList = list()
cellList[[24]] = cell24



c = names(cellList) = c('cell01','cell02','cell03','cell04','cell05','cell06','cell07','cell08','cell09','cell10','cell11','cell12','cell13','cell14','cell15','cell16','cell17','cell18','cell19','cell20','cell21','cell22','cell23','cell24')


save(cellList, file = 'f:/_R_WD/_data_frames/IL_juxta_GlyT2_fiber_stim/cellList.RData')

load('f:/_R_WD/_data_frames/IL_juxta_GlyT2_fiber_stim/cellList.RData')


    ### Corrected cellList: stims with no controls are removed
    
      cellList$cell01 = cellList$cell01[-1,]
      cellList$cell15 = cellList$cell15[-3:-4,]
      cellList$cell16 = cellList$cell16[-3:-4,]
      cellList$cell16 = cellList$cell16[-5:-6,]
      cellList$cell17 = cellList$cell17[-3:-4,]
      cellList$cell17 = cellList$cell17[-5:-6,]
      cellList$cell18 = cellList$cell18[-1,]
      cellList$cell19 = cellList$cell19[-1,]
    
      save(cellList,file = 'f:/_R_WD/_data_frames/IL_juxta_GlyT2_fiber_stim/Corrected_cellList.RData')

      
################################################

load('f:/_R_WD/_data_frames/IL_juxta_GlyT2_fiber_stim/Corrected_cellList.RData')

      
load(paste0(file_path,
            '/Corrected_cellList.RData'
            )
     )
   

v = formatC(seq(1,24), width = 2, flag = 0)
#sprintf()


ReArrange = function(v) {

#  cellName = paste0('cell', readline(prompt = 'Give me the input cell number (1 to 24): '))
  
  for (i in v) {
  cellName = paste0('cell', i)
  
  
  input = cellList[[cellName]]

  dim(input)

  rows = dim(input)[1]
  cols = dim(input)[2]   
  
  dim(input) = c(rows * cols,1)
  
  input = cbind(input, matrix(0, rows * cols , 1))
  input[1:rows,2] = 'b'
  input[(rows+1):(2*rows),2] = 'd'
  input[(2*rows+1):(3*rows),2] = 'a'
  
  input = cbind(input, matrix(as.character(cellName) ,rows * cols,1 ))
  colnames(input) = c('No.AP','condition','cell')

  assign(cellName, input, envir = globalenv())
  }
}  

ReArrange(v)

#allCells = mget(c)
c = c('cell01','cell02','cell03','cell04','cell05','cell06','cell07','cell08','cell09','cell10','cell11','cell12','cell13','cell14','cell15','cell16','cell17','cell18','cell19','cell20','cell21','cell22','cell23','cell24')


allCells = do.call(rbind, mget(c))

allCells = as.data.frame(allCells)
allCells$pinch = T

allCells$pinch[allCells$cell == 'cell18' |
                 allCells$cell == 'cell19' |
                 allCells$cell == 'cell20' |
                 allCells$cell == 'cell21' |
                 allCells$cell == 'cell23' |
                 allCells$cell == 'cell01' |
                 allCells$cell == 'cell06'] = FALSE

#subset(allCells$pinch, allCells$cell == 'cell18' | allCells$cell == 'cell19')


allCells$firing = 'regular'
allCells$firing[allCells$cell == 'cell02' |
                allCells$cell == 'cell03' |
                allCells$cell == 'cell04' |
                allCells$cell == 'cell05' |
                allCells$cell == 'cell10' |
                allCells$cell == 'cell14' |
                allCells$cell == 'cell15' |
                allCells$cell == 'cell16' |
                allCells$cell == 'cell17' |
                allCells$cell == 'cell24'] = 'irregular'

allCells$stimLength = 0
allCells$stimLength[allCells$cell == 'cell01' |
                  allCells$cell == 'cell02' |
                  allCells$cell == 'cell03'] = 33

allCells$stimLength[allCells$cell == 'cell04' |
                      allCells$cell == 'cell05' |
                      allCells$cell == 'cell06' |
                      allCells$cell == 'cell08' |
                      allCells$cell == 'cell09' ] = 16.47

allCells$stimLength[allCells$cell == 'cell16' |
                      allCells$cell == 'cell17' |
                      allCells$cell == 'cell18'] = 6.57

allCells$stimLength[allCells$cell == 'cell11' |
                  allCells$cell == 'cell12' |
                  allCells$cell == 'cell13' |
                  allCells$cell == 'cell14' |
                  allCells$cell == 'cell15' |
                  allCells$cell == 'cell19' |
                  allCells$cell == 'cell20' |
                  allCells$cell == 'cell21' |
                  allCells$cell == 'cell22' |
                  allCells$cell == 'cell23' |
                  allCells$cell == 'cell24' ] = 4.92



allCells$stimLength[allCells$cell == 'cell07' ] = 34.62
allCells$stimLength[allCells$cell == 'cell10' ] = 3.27

  
save(allCells, file = 'f:/_R_WD/_data_frames/IL_juxta_GlyT2_fiber_stim/allCells.RData')


##########################  Plotting  ###################################


load(file = 'f:/_R_WD/_data_frames/IL_juxta_GlyT2_fiber_stim/allCells.RData')
load('f:/_R_WD/useful_to_load/colorMatrix.RData')

load('C:/Users/Viktor/OneDrive - MTA KOKI/professional/_R_WD/useful_to_load/colorMatrix.RData')

load(paste0(file_path,
            '/allCells.RData'
            )
     )

str(allCells)

### convert to numeric
allCells$No.AP = as.numeric(levels(allCells$No.AP))[allCells$No.AP]
#allCells$No.AP = as.numeric(as.character(allCells$No.AP))

allCells$firing = as.factor(allCells$firing)


#all
ver.A = allCells

#pinched and regular firing
ver.B = subset(allCells, allCells$pinch == TRUE & allCells$firing == 'regular')

levels(ver.B$cell)
#ver.B$cell = factor(ver.B$cell)
length(levels(droplevels(ver.B$cell)))


#pinched and irregular firing
ver.C = subset(allCells, allCells$pinch == TRUE & allCells$firing == 'irregular')
levels(ver.C$cell)
length(levels(droplevels(ver.C$cell)))



#spontaneous and regular firing
ver.D = subset(allCells, allCells$pinch == FALSE & allCells$firing == 'regular')
levels(ver.D$cell)
levels(droplevels(ver.D$cell))



#spontaneous and irregular firing
ver.E = subset(allCells, allCells$pinch == FALSE & allCells$firing == 'irregular')

par(mfrow = c(2,3))
par(mfrow = c(1,1))

VER = ver.A

plot.No.AP = function(){
with(VER, 
     interaction.plot(factor(condition, levels = c('b','d','a')
                             ), 
                      cell, 
                      No.AP / stimLength,
                      fun = mean,
                      fixed = F,
                      type = 'b',
                      lty = 1,
                      pch = 1,
                      bty = 'L',
                      #col = as.factor(pinch),
                      col = colorMatrix['orange','dark'],
                      legend = F
                      )
     )


with(VER,
     boxplot(No.AP / stimLength ~ factor(condition,
                                     levels = c('b','d','a')
     ),
     boxwex = 0.1,
     add = T,
     outline = F,
     col= rgb(0,0,0, maxColorValue = 255 ,alpha=0)
     )
) 

}     
plot.No.AP()




with(subset(allCells, allCells$pinch == TRUE & allCells$firing == 'regular'), 
     stripchart(No.AP / stimLength ~ factor(condition, levels = c('b','d','a')
                                                 ),
                     vertical = T,
                     pch = 19,
                     method = 'jitter',
                     jitter = 0.1,
                     col = colorMatrix['blue','dark']
                     )
     )
     
with(subset(allCells, allCells$pinch == T & allCells$firing == 'irregular'), 
     stripchart(No.AP / stimLength ~ factor(condition, levels = c('b','d','a')
     ),
     vertical = T,
     pch = 19,
     method = 'jitter',
     jitter = 0.1,
     add =T,
     col = colorMatrix['blue','midle']
     )
)

with(subset(allCells, allCells$pinch == F & allCells$firing == 'regular'), 
     stripchart(No.AP / stimLength ~ factor(condition, levels = c('b','d','a')
     ),
     vertical = T,
     pch = 19,
     method = 'jitter',
     jitter = 0.1,
     add =T,
     col = colorMatrix['blue','light']
     )
)



### ggplot2 -------------------------------------------------------------

load(file = "f:/_R_WD/_data_frames/IL_juxta_GlyT2_fiber_stim/allCells.RData")
load(file = "f:/_R_WD/useful_to_load/colorMatrix.RData")
library(tidyverse)


#filter_condititon <- quote(pinch == F, firing == 'regular')
filter_condititon <- "pinch == T" #& firing == 'irregular'"


### transforming allCells data frame
allCells_for_plot <- allCells %>%
  dplyr::filter(eval(parse(text = filter_condititon))) %>% 
  mutate(No.AP = as.numeric(levels(No.AP))[No.AP]) %>%
  mutate(No.AP_p_stimLength = No.AP / stimLength) %>%
  dplyr::group_by(condition, cell) %>%
  summarise(mean_No.AP = mean(No.AP_p_stimLength))

### select cells to highlight on plot
cells_to_highlight <- allCells %>%
  filter(eval(parse(text = filter_condititon))) %>% 
  mutate(No.AP = as.numeric(levels(No.AP))[No.AP]) %>%
  mutate(No.AP_p_stimLength = No.AP / stimLength) %>%
  dplyr::group_by(condition, cell) %>%
  summarise(mean_No.AP = mean(No.AP_p_stimLength)) %>%
  filter(cell == "cell01" |
           cell == "cell06" |
           cell == "cell07" |
           cell == "cell08" |
           cell == "cell09" |
           cell == "cell11" |
           cell == "cell16")



### PLOT: MFR before/during/after (transformed allCells data frame)
ggplot(data = allCells_for_plot,
       mapping = aes(x = forcats::fct_relevel(condition, "b", "d", "a"),
                     y = mean_No.AP)) +
  #theme(panel.background = element_rect(fill = 0)) +
  theme_minimal() +
  theme(axis.text = element_text(size = 20),
        text = element_text(size = 20)) +
  geom_boxplot(width = 0.2, alpha = 0.5) +
  geom_point(shape = 21, 
             fill = "#EB8104", 
             #color = "white",
             size = 4) +
             #stroke = 2) +
  # geom_point(data = cells_to_highlight, 
  #            #shape = 21, 
  #            color = "#1D4871", 
  #            #color = "white",  
  #            size = 2) +
  #            #stroke = 2) +
  geom_line(aes(group = cell), color = "#EB8104") +
  # geom_line(data = cells_to_highlight, color = "#1D4871", aes(group = cell)) +
  # geom_text(data = cells_to_highlight, 
  #           aes(label = cell),
  #           color = "#1D4871",
  #           hjust = -0.2, vjust = -0.2) +
  scale_x_discrete(name = "Stimulus",labels = c("Before", "During", "After")) +
  labs(y = "Mean firing rate") 
# scale_y_continuous(sec.axis = sec_axis(~.*2, name = "proba axis"))

ggsave(file.path("output_data","inhibition_of_IL_cells_pinch.png"),
       width = 8,
       height = 12,
       dpi = 300)





#1D4871
#EB8104