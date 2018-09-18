### recording site: IL 
### stimulation: glycinergic fibers in the IL around the recorded cells
### to create population PSTHs across animals. APs are aligned to individual stimuli



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
                       'IL_juxta_GlyT2_fiber_stim',
                       'cell21_13_right_IL_t3_3030_stim_control.mat')


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
raw.rec <- readMat(file_path)



### exploring the list names
names(raw.rec)

### renaming lists (to make it easier to work with later). Different mat files can have different list names (each needs to be identified)
names(raw.rec) <- c('unit','level','stim','eeg') 
names(raw.rec) <- c('k','na','level','stim','unit','eeg')
names(raw.rec) <- c('level', 'stim', 'unit', 'unit2', 'eeg')



### AP: stores the time of action potentials in a matrix
AP = raw.rec$unit[[12]] 


### stim: stores the time of the stimuli in a matrix
stim = raw.rec$stim[[5]] 


### LFP: Stores the amplitude values of the cortical oscillation in a matrix. Amplitude changes over time and it is sampled at 20 kHz (can be handled as a time series object)
### scale function: standardization (mean = 0, SD = 1)
LFP = scale(raw.rec$eeg[[9]]) 


### calculates the sampling rate of the LFP from the time between two data points stored in raw.rec$eeg[[3]]. Stores the value in a vector
samp.rate.lfp <- c(1/raw.rec$eeg[[3]])


### calculates the sampling rate of the unit from the time between two data points stored in raw.rec$unit[[4]]. Stores the value in a vector
samp.rate.unit <- c(1/raw.rec$unit[[4]])


### calculates the length of the recording from the number of data points of the eeg and the sampling rate. Stores it in a vector
rec.length <- length(raw.rec$eeg[[9]]) / samp.rate.lfp


### creates a matrix from AP times (repeats them)
mat_AP <- matrix(AP, 
                 nrow = length(stim), 
                 ncol = length(AP), byrow = T
)
dim(mat_AP)



### creates a matrix from stim times (repeats them)
mat_stim <- matrix(stim, 
                   nrow = length(stim),
                   ncol = length(AP)
)
dim(mat_stim)


### subtracts the two matrices to calculate the times of every action potentials relative to each stimuli
mat_reltimes <- mat_AP - mat_stim


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
named.list <- function(...) { 
  l <- list(...)
  names(l) <- sapply(substitute(list(...)), deparse)[-1]
  l 
}

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

load_path <- file.path('f:',
                       '_R_WD',
                       '_data_frames',
                       'IL_juxta_GlyT2_fiber_stim',
                       'IL_inhibition_popPSTH.RData'
)

load(load_path)

named.list <- function(...) { 
  l <- list(...)
  names(l) <- sapply(substitute(list(...)), deparse)[-1]
  l 
}

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
binsize <- 15
X_lim <- c(-0.01,0.02)
Y_lim <- c(0,250)

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