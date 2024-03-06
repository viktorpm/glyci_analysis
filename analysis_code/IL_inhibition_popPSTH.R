### recording site: IL 
### stimulation: glycinergic fibers in the IL around the recorded cells
### to create population PSTHs across animals. APs are aligned to individual stimuli



########################## Reading and organizing data ##########################

### function: separated by CapitalLetters
### data from external source: separated.by '.'
### variable: separated_by '_'

library(tidyverse)

### source files from the utilities folder
walk(dir("supplementary_functions"), ~source(file.path("supplementary_functions", .x)))
source(file.path("supplementary_functions","ReadMatFile.R"))


ReadMatFile('cell21_13_right_IL_t3_3030_stim_control.mat')


### calculates the time difference between stimuli
d1 <- diff(stim)


### calculates and rounds the frequency of the stimuli
stim_freq <- round(1/d1, 0) 
length(stim_freq)
table(as.data.frame(stim_freq))


### removes gaps (inter stimulus train intervals, frequencies smaller than 1)????? same length as before removal???
stim_freq[which(stim_freq < 1)] = stim_freq[which(stim_freq < 1) -1]
table(as.data.frame(stim_freq))


### 
stim_freq <- append(stim_freq[,1], stim_freq[length(stim_freq),1])



levels(as.factor(stim_freq))
#count(as.factor(stim_freq))


### binds the rel_times matrix and the frequency matrix
mat_reltimes <- cbind(stim_freq, mat_reltimes) 
dim(mat_reltimes)

### renames columns
colnames(mat_reltimes) = c('HZ',
                           paste('AP', 
                                 1:length(mat_reltimes[1,2:length(mat_reltimes[1,])]),
                                 sep = '')
) 



### transposes reltimes matrix and converts it to data frame
mat_reltimes_df <- as.data.frame(t(mat_reltimes))

# run code on every cell, save results in R.Data file 


cell02 = named.list(mat_reltimes, d1, stim_freq)

save(cell01,
     cell02,
     cell03,
     cell04,
     cell05,
     cell06,
     cell07,
     cell08,
     cell09,
     cell10,
     cell11,
     cell12,
     cell13,
     cell14,
     cell15,
     cell16,
     cell17,
     cell18,
     cell19,
     cell20,
     cell21,
     cell22,
     cell23,
     cell24,
     file = 'IL_inhibition_popPSTH.RData')

############################################################################
################################ Plotting ##################################
############################################################################
path <- file.path('r_course_data')

load(file.path(path,
               'IL_inhibition_popPSTH.RData'
               )
     )



all.cell <- named.list(cell01$mat_reltimes[cell01$mat_reltimes > -1 & cell01$mat_reltimes < 1],
                       cell02$mat_reltimes[cell02$mat_reltimes > -1 & cell02$mat_reltimes < 1],
                       cell03$mat_reltimes[cell03$mat_reltimes > -1 & cell03$mat_reltimes < 1],
                       cell04$mat_reltimes[cell04$mat_reltimes > -1 & cell04$mat_reltimes < 1],
                       cell05$mat_reltimes[cell05$mat_reltimes > -1 & cell05$mat_reltimes < 1],
                       cell06$mat_reltimes[cell06$mat_reltimes > -1 & cell06$mat_reltimes < 1],
                       cell07$mat_reltimes[cell07$mat_reltimes > -1 & cell07$mat_reltimes < 1],
                       cell08$mat_reltimes[cell08$mat_reltimes > -1 & cell08$mat_reltimes < 1],
                       cell09$mat_reltimes[cell09$mat_reltimes > -1 & cell09$mat_reltimes < 1],
                       cell10$mat_reltimes[cell10$mat_reltimes > -1 & cell10$mat_reltimes < 1],
                       cell11$mat_reltimes[cell11$mat_reltimes > -1 & cell11$mat_reltimes < 1],
                       cell12$mat_reltimes[cell12$mat_reltimes > -1 & cell12$mat_reltimes < 1],
                       cell13$mat_reltimes[cell13$mat_reltimes > -1 & cell13$mat_reltimes < 1],
                       cell14$mat_reltimes[cell14$mat_reltimes > -1 & cell14$mat_reltimes < 1],
                       cell15$mat_reltimes[cell15$mat_reltimes > -1 & cell15$mat_reltimes < 1],
                       cell16$mat_reltimes[cell16$mat_reltimes > -1 & cell16$mat_reltimes < 1],
                       cell17$mat_reltimes[cell17$mat_reltimes > -1 & cell17$mat_reltimes < 1],
                       cell18$mat_reltimes[cell18$mat_reltimes > -1 & cell18$mat_reltimes < 1],
                       cell19$mat_reltimes[cell19$mat_reltimes > -1 & cell19$mat_reltimes < 1],
                       cell20$mat_reltimes[cell20$mat_reltimes > -1 & cell20$mat_reltimes < 1],
                       cell21$mat_reltimes[cell21$mat_reltimes > -1 & cell21$mat_reltimes < 1],
                       cell22$mat_reltimes[cell22$mat_reltimes > -1 & cell22$mat_reltimes < 1],
                       cell23$mat_reltimes[cell23$mat_reltimes > -1 & cell23$mat_reltimes < 1],
                       cell24$mat_reltimes[cell24$mat_reltimes > -1 & cell24$mat_reltimes < 1])



all.cell <- named.list(cell06$mat_reltimes[cell06$mat_reltimes > -1 & cell06$mat_reltimes < 1],
                       cell07$mat_reltimes[cell07$mat_reltimes > -1 & cell07$mat_reltimes < 1],
                       cell08$mat_reltimes[cell08$mat_reltimes > -1 & cell08$mat_reltimes < 1],
                       cell09$mat_reltimes[cell09$mat_reltimes > -1 & cell09$mat_reltimes < 1])


all.cell <- named.list(cell01$mat_reltimes[cell01$mat_reltimes > -1 & cell01$mat_reltimes < 1],
                       cell02$mat_reltimes[cell02$mat_reltimes > -1 & cell02$mat_reltimes < 1],
                       cell03$mat_reltimes[cell03$mat_reltimes > -1 & cell03$mat_reltimes < 1],
                       cell04$mat_reltimes[cell04$mat_reltimes > -1 & cell04$mat_reltimes < 1],
                       cell05$mat_reltimes[cell05$mat_reltimes > -1 & cell05$mat_reltimes < 1],
                       cell10$mat_reltimes[cell10$mat_reltimes > -1 & cell10$mat_reltimes < 1],
                       cell11$mat_reltimes[cell11$mat_reltimes > -1 & cell11$mat_reltimes < 1],
                       cell12$mat_reltimes[cell12$mat_reltimes > -1 & cell12$mat_reltimes < 1],
                       cell13$mat_reltimes[cell13$mat_reltimes > -1 & cell13$mat_reltimes < 1],
                       cell14$mat_reltimes[cell14$mat_reltimes > -1 & cell14$mat_reltimes < 1],
                       cell15$mat_reltimes[cell15$mat_reltimes > -1 & cell15$mat_reltimes < 1],
                       cell16$mat_reltimes[cell16$mat_reltimes > -1 & cell16$mat_reltimes < 1],
                       cell17$mat_reltimes[cell17$mat_reltimes > -1 & cell17$mat_reltimes < 1],
                       cell18$mat_reltimes[cell18$mat_reltimes > -1 & cell18$mat_reltimes < 1],
                       cell19$mat_reltimes[cell19$mat_reltimes > -1 & cell19$mat_reltimes < 1],
                       cell20$mat_reltimes[cell20$mat_reltimes > -1 & cell20$mat_reltimes < 1],
                       cell21$mat_reltimes[cell21$mat_reltimes > -1 & cell21$mat_reltimes < 1],
                       cell22$mat_reltimes[cell22$mat_reltimes > -1 & cell22$mat_reltimes < 1],
                       cell23$mat_reltimes[cell23$mat_reltimes > -1 & cell23$mat_reltimes < 1],
                       cell24$mat_reltimes[cell24$mat_reltimes > -1 & cell24$mat_reltimes < 1])


rm(list = ls(pattern = "cell0")) 
rm(list = ls(pattern = "cell1")) 
rm(list = ls(pattern = "cell2")) 

all.cell.vector <- unlist(all.cell)
all.cell.vector <- all.cell.vector + 0.0005


# numbers <- formatC(seq(1,24), width = 2, flag = 0)
# cell.names <- paste0('cell',numbers)
# 
# get(cell.names[1][[2]])
# 
# for (i in length(cell.names)){
#   all.cell[[i]] <- eval(parse(text = cell.names[i][[1]]))
# }
# 
# eval(parse(text = cell.names[1][[1]]))
# 
# for (i in length(cell.names)){
#   all.cell[[i]] <- mget(cell.names[i])
# }



hist_range_lower <- - 0.01
hist_range_upper <- 0.02
binsize <- 30
X_lim <- c(-0.01,0.02)
Y_lim <- c(0,140)

par(mfrow = c(1,2))

hist(all.cell.vector[all.cell.vector > hist_range_lower & all.cell.vector < hist_range_upper], 
     breaks = binsize,
     xlim = X_lim
     ,ylim = Y_lim
)
abline(v = 0)

#      #add = T, 
#      col = col_slow_stim,
#      lty = 0,
#      freq = T,
#      xlab = 'Latency (s)',
#      main = NULL
# )