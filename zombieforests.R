#Zombie forests

#This file runs the analyses and generates the plots for the main manuscript.


#Set up directory
library(here)
here::i_am("zombieforests.R")
main = here()

#Get numerics & functions
source("zombiefunctions.R")
source("zombieNumerics.R")

#To run analyses in parallel
library(parallel)
library(pbapply)

#If this is your first time running this file, set 'firstRun' to TRUE
#It will take hours to run the analyses the first time
firstRun = FALSE

################################################################################
# MAIN ANALYSES
################################################################################

#Runs one iteration of the code that finds the critical rate of climate change
#for a given parameterization

iterate <- function(subintervals, AInside, ABefore, AAfter, dispersalRate, 
                    precision, guess, domainStart, domainLength, patchLength, 
                    accuracy, main_dir){
  
  source("zombieNumerics.R")
  source("zombiefunctions.R")
  
  stages = dim(AInside)[1]
  
  K = matrix(data=0, nrow=stages, ncol=stages)
  K[1,stages] = dispersalRate
  
  critClimRate <- getCritClimRate(subintervals = subintervals,
                                  AInside = AInside,
                                  ABefore = ABefore,
                                  AAfter = AAfter,
                                  K = K,
                                  precision = precision,
                                  guess = guess,
                                  domainStart = domainStart,
                                  domainLength = domainLength,
                                  patchLength = patchLength,
                                  accuracy=accuracy) 
  return (critClimRate)
}

################################################################################
# FIGURE GENERATION
################################################################################

################################################################################
#Figure 1: Conceptual diagram <- <- <- <- <- <- <- <- <- <- <-<- <- <- <- <- <-
################################################################################

#Shell of structure

library(viridis)
colors = viridis(8)
time1 = seq(from=0, to = 10, length.out=100)
lowerPatchLimit = seq(from=0, to=10, length.out=100)
time2 = seq(from=0, to = 8, length.out=100)
upperPatchLimit = seq(from=2, to=10, length.out=100)

sub_dir = "Figures"
if (!dir.exists(sub_dir)){
  dir.create(sub_dir)
}
setwd(file.path(main, sub_dir))
svg(filename="fig1part1.svg", width=5, height=4, bg="transparent")
plot(time1, lowerPatchLimit, type="l", lty=2.0, col="darkgrey", xlab="time", 
     ylab="", yaxt='n', ann=FALSE, frame.plot=TRUE, lwd=2.0, cex.lab=0.8, 
     cex.axis=0.8)
polygon(c(0,10,10), c(0,0,10), col="lightgrey",border=NA)
polygon(c(0,8,0), c(2,10,10), col="lightgrey", border=NA)
lines(time1, lowerPatchLimit, lwd=2.0, lty=2.0, col="darkgrey")
lines(time2, upperPatchLimit, lty=2.0, col="darkgrey", lwd=2.0)
title(xlab="time", line=2, cex.lab=1.0)
title(ylab="x", line=1, cex.lab=1.0)

#Add in "patches"
for (i in 1:8){
  segments(i, i, i, i+2, col = colors[i], lwd = 3.0, xpd = FALSE)
}
segments(3,2,3,3,col="black", lty=3.0, lwd=2.0)
segments(2,2,3,2,col="black", lty=3.0, lwd=2.0)
dev.off()
setwd(main)

################################################################################

#Dispersal Kernel
sub_dir = "Figures"
if (!dir.exists(sub_dir)){
  dir.create(sub_dir)
}
setwd(file.path(main, sub_dir))
svg(filename="fig1part2.svg", width=2, height=4, bg="transparent")
xvalues=seq(from=-10,to=10, length.out=1000)
yvalues=-1*sapply(xvalues,laplace,y=0, b=1.0)
plot(yvalues,xvalues, type="l", axes=FALSE, xlab="", ylab="", lwd=1.5)
dev.off()
setwd(main)

################################################################################
#Figure 2: Example of critical speed of climate change <- <- <- <- <- <- <- <- <
################################################################################

#The example is for Sugar Maple, part a shows survival, part b shows extirpation

#Set up the dispersal kernel
K = matrix(data=0, nrow=3, ncol=3)
K[1,3] = 5.0

getCritClimRate(subintervals=1000, AInside=sugarMaple, ABefore=sugarMapleIV, 
                AAfter=sugarMapleI, K=K, precision=0.001, guess=1.0,
                domainStart=-100, domainLength=300, patchLength=10, 
                accuracy="accurate")

#C = 0.9 for survive,  0.9111328/0.916992 is critical speed
#C = 1.0 for die

survive = runModel(subintervals = 1000, 
                   timesteps = 25, 
                   timestepsCC = 100, 
                   c = 0.91,
                   AInside = sugarMaple,
                   ABefore = sugarMapleIV, 
                   AAfter = sugarMapleI,
                   K = K,
                   domainStart = -50, 
                   domainLength = 200, 
                   patchLength = 10,
                   accuracy="accurate")

die = runModel(subintervals = 1000,
               timesteps = 25,
               timestepsCC = 100,
               c = 1.1, 
               AInside = sugarMaple,
               ABefore = sugarMapleIV, 
               AAfter = sugarMapleI,
               K = K,
               domainStart = -50,
               domainLength = 200,
               patchLength = 10,
               accuracy="accurate")

#Chopping off excess space to make it look nicer
survive = survive[survive[,2]>=-25 & survive[,2]<=125, ]
die = die[die[,2]>=-25 & die[,2]<=125, ]

x11(width=5, height=4.5)
# Margins area
par(oma=c(1,1,1,1)) #outer margin area (bot, left, top, right) (space around)
par(mar=c(2,4,1,1)) #margin area (bot, left, top, right)
makePlot(survive, byGeneration=5, stage=4, xlab="")
plotA <- grab_grob()
dev.off()

x11(width=5, height=4.5)
par(oma=c(0,1,1,1)) #outer margin area (bot, left, top, right) (space around)
par(mar=c(6,4,1,1)) #margin area (bot, left, top, right)
makePlot(die, byGeneration=5, stage=4)
plotB <- grab_grob()
dev.off()


sub_dir = "Figures"
if (!dir.exists(sub_dir)){
  dir.create(sub_dir)
}
setwd(file.path(main, sub_dir))
pdf(file="fig2.pdf", width=5, height=8)
library(gridExtra)
grid.arrange(plotA,plotB, ncol=1) 
grid.text(c("(a)","(b)"), x=c(0.05, 0.05), y=c(0.95, 0.47), 
          vjust=0, hjust=0, gp=gpar(fontface=1))
dev.off()
setwd(main)

################################################################################
#Figure 3AB: Dispersal and reductions in survival/fecundity <- <- <- <- <- <- <- 
################################################################################

#In 3A (sugar maple) and 3B (white fir), the environment behind the patch 
#varies, but is constant and zero in front of the patch

#The other habitat configurations are examined in supplement S3

#Dispersal rates for x-axis, for sugar Maple
dispersalRates <- seq(0,20,0.5)

#One line for each environment type
lambdaSugarMapleI <- rep(0,41)
lambdaSugarMapleII <- rep(0,41)
lambdaSugarMapleIII <- rep(0,41)
lambdaSugarMapleIV <- rep(0,41)


#Environmental scenario 1
numCores = parallel::detectCores(logical=FALSE)
clus <- parallel::makeCluster(numCores-2)
lambdaSugarMapleI <- pbapply::pblapply(cl=clus, X=dispersalRates, FUN=iterate, 
                                       subintervals = 1000,
                                       AInside = sugarMaple, 
                                       ABefore = sugarMapleI,
                                       AAfter = sugarMapleI, 
                                       precision = 0.001, guess = 2.0,
                                       domainStart = -150, 
                                       domainLength = 500, 
                                       patchLength = 10,
                                       accuracy = "accurate",
                                       main_dir = main)
if (firstRun){
  sub_dir = "Data/Critical Speeds of Climate Change"
  if (!dir.exists(sub_dir)){
    dir.create(sub_dir)
  }
  setwd(file.path(main, sub_dir))
  saveRDS(lambdaSugarMapleI, file = "LSMI.rds")
}

#Environmental scenario 2
numCores = parallel::detectCores(logical=FALSE)
clus <- parallel::makeCluster(numCores-2)
lambdaSugarMapleII <- pbapply::pblapply(cl=clus, X=dispersalRates, FUN=iterate, 
                                        subintervals = 1000,
                                        AInside = sugarMaple, 
                                        ABefore = sugarMapleII,
                                        AAfter = sugarMapleI, 
                                        precision = 0.001, guess = 2.0,
                                        domainStart = -150, 
                                        domainLength = 500, 
                                        patchLength = 10,
                                        accuracy = "accurate",
                                        main_dir = main)
if (firstRun){
  sub_dir = "Data/Critical Speeds of Climate Change"
  if (!dir.exists(sub_dir)){
    dir.create(sub_dir)
  }
  setwd(file.path(main, sub_dir))
  saveRDS(lambdaSugarMapleII, file = "LSMII.rds")
}

#Environmental scenario 3
numCores = parallel::detectCores(logical=FALSE)
clus <- parallel::makeCluster(numCores-2)
lambdaSugarMapleIII <- pbapply::pblapply(cl=clus, X=dispersalRates, FUN=iterate, 
                                         subintervals = 1000,
                                         AInside = sugarMaple, 
                                         ABefore = sugarMapleIII,
                                         AAfter = sugarMapleI, 
                                         precision = 0.001, guess = 2.0,
                                         domainStart = -150, 
                                         domainLength = 500, 
                                         patchLength = 10,
                                         accuracy = "accurate",
                                         main_dir = main)
if (firstRun){
  sub_dir = "Data/Critical Speeds of Climate Change"
  if (!dir.exists(sub_dir)){
    dir.create(sub_dir)
  }
  setwd(file.path(main, sub_dir))
  saveRDS(lambdaSugarMapleIII, file = "LSMIII.rds")
}

#Environmental scenario 4
numCores = parallel::detectCores(logical=FALSE)
clus <- parallel::makeCluster(numCores-2)
lambdaSugarMapleIV <- pbapply::pblapply(cl=clus, X=dispersalRates, FUN=iterate, 
                                        subintervals = 1000,
                                        AInside = sugarMaple, 
                                        ABefore = sugarMapleIV,
                                        AAfter = sugarMapleI, 
                                        precision = 0.001, guess = 2.0,
                                        domainStart = -150, 
                                        domainLength = 500, 
                                        patchLength = 10,
                                        accuracy = "accurate",
                                        main_dir = main)
if (firstRun){
  sub_dir = "Data/Critical Speeds of Climate Change"
  if (!dir.exists(sub_dir)){
    dir.create(sub_dir)
  }
  setwd(file.path(main, sub_dir))
  saveRDS(lambdaSugarMapleIV, file = "LSMIV.rds")
}


#White fir analysis

#Set up dispersal rates, and one line for each environment type
dispersalRates <- seq(0,7.0,0.1)
lambdaWhiteFirI <- rep(0,71)
lambdaWhiteFirII <- rep(0,71)
lambdaWhiteFirIII <- rep(0,71)
lambdaWhiteFirIV <- rep(0,71)

#Environmental scenario 1
numCores = parallel::detectCores(logical=FALSE)
clus <- parallel::makeCluster(numCores-2)
lambdaWhiteFirI <- pbapply::pblapply(cl=clus, X=dispersalRates, FUN=iterate, 
                                     subintervals = 1000,
                                     AInside = whiteFir, 
                                     ABefore = whiteFirI,
                                     AAfter = whiteFirI, 
                                     precision = 0.0001, guess = 0.20,
                                     domainStart = -50, domainLength = 200, 
                                     patchLength = 1,
                                     accuracy = "accurate",
                                     main_dir = main)
if (firstRun){ 
  sub_dir = "Data/Critical Speeds of Climate Change"
  if (!dir.exists(sub_dir)){
    dir.create(sub_dir)
  }
  setwd(file.path(main, sub_dir))
  saveRDS(lambdaWhiteFirI, file = "LWFI.rds")
}

#Environmental scenario 2
numCores = parallel::detectCores(logical=FALSE)
clus <- parallel::makeCluster(numCores-2)
lambdaWhiteFirII <- pbapply::pblapply(cl=clus, X=dispersalRates, FUN=iterate, 
                                      subintervals = 1000,
                                      AInside = whiteFir, 
                                      ABefore = whiteFirII,
                                      AAfter = whiteFirI, 
                                      precision = 0.0001, guess = 0.20,
                                      domainStart = -50, domainLength = 200, 
                                      patchLength = 1,
                                      accuracy = "accurate",
                                      main_dir = main)
if (firstRun){
  sub_dir = "Data/Critical Speeds of Climate Change"
  if (!dir.exists(sub_dir)){
    dir.create(sub_dir)
  }
  setwd(file.path(main, sub_dir))
  saveRDS(lambdaWhiteFirII, file = "LWFII.rds")
}

#Environmental scenario 3
numCores = parallel::detectCores(logical=FALSE)
clus <- parallel::makeCluster(numCores-2)
lambdaWhiteFirIII <- pbapply::pblapply(cl=clus, X=dispersalRates, FUN=iterate, 
                                       subintervals = 1000,
                                       AInside = whiteFir, 
                                       ABefore = whiteFirIII,
                                       AAfter = whiteFirI, 
                                       precision = 0.0001, guess = 0.20,
                                       domainStart = -50, 
                                       domainLength = 200, 
                                       patchLength = 1,
                                       accuracy = "accurate",
                                       main_dir = main)
if (firstRun){
  sub_dir = "Data/Critical Speeds of Climate Change"
  if (!dir.exists(sub_dir)){
    dir.create(sub_dir)
  }
  setwd(file.path(main, sub_dir))
  saveRDS(lambdaWhiteFirIII, file = "LWFIII.rds")
}

#Environmental scenario 4
numCores = parallel::detectCores(logical=FALSE)
clus <- parallel::makeCluster(numCores-2)
lambdaWhiteFirIV <- pbapply::pblapply(cl=clus, X=dispersalRates, FUN=iterate, 
                                      subintervals = 1000,
                                      AInside = whiteFir, 
                                      ABefore = whiteFirIV,
                                      AAfter = whiteFirI, 
                                      precision = 0.0001, guess = 0.20,
                                      domainStart = -50, domainLength = 200, 
                                      patchLength = 1,
                                      accuracy = "accurate",
                                      main_dir = main)
if (firstRun){
  sub_dir = "Data/Critical Speeds of Climate Change"
  setwd(file.path(main, sub_dir))
  if (!dir.exists(sub_dir)){
    dir.create(sub_dir)
  }
  saveRDS(lambdaWhiteFirIV, file = "LWFIV.rds")
}


################################################################################

#Constructing plot
sub_dir = "Data/Critical Speeds of Climate Change"
setwd(file.path(main, sub_dir))
lambdaSugarMapleI = readRDS(file = "LSMI.rds")
lambdaSugarMapleII = readRDS(file = "LSMII.rds")
lambdaSugarMapleIII = readRDS(file = "LSMIII.rds")
lambdaSugarMapleIV = readRDS(file = "LSMIV.rds")
lambdaWhiteFirI = readRDS(file = "LWFI.rds")
lambdaWhiteFirII = readRDS(file = "LWFII.rds")
lambdaWhiteFirIII = readRDS(file = "LWFIII.rds")
lambdaWhiteFirIV = readRDS(file = "LWFIV.rds")

sub_dir = "Figures"
if (!dir.exists(sub_dir)){
  dir.create(sub_dir)
}
setwd(file.path(main, sub_dir))

svg(file="fig3.svg", width=7, height=10)

#(a) Sugar Maple

#x-axis
dispersalRates <- seq(0,20,0.5)

par(mfrow = c(2,1)) #2 rows, 1 column
par(oma=c(1,1,0,1)) #outer margin area (bot, left, top, right) (space around)
par(mar=c(4,4,0,0)) #margin area (bot, left, top, right)

col1 = "#a6cee3" #light blue, environmental scenario 1
col2 = "#33a02c" #dark green, environmental scenario 2
col3 = "#b2df8a" #light green, environmental scenario 3
col4 = "#1f78b4" #dark blue, environmental scenario 4

#Set up plot and lines
plot(dispersalRates, lambdaSugarMapleII, type='l', col=col2, xlab='',
     ylab='critical speed of climate change (m/yr)', lwd=3.0, cex.axis=0.8) 
lines(dispersalRates, lambdaSugarMapleIV, col=col4, lwd=3.0) 
lines(dispersalRates, lambdaSugarMapleIII, col=col3, lwd=3.0) 
lines(dispersalRates, lambdaSugarMapleI, col=col1, lwd=3.0) 

#Add symbols (because some lines coincide)
points(dispersalRates, lambdaSugarMapleI, pch=1)
points(dispersalRates, lambdaSugarMapleII, pch=2)
points(dispersalRates, lambdaSugarMapleIII, pch=3)
points(dispersalRates, lambdaSugarMapleIV, pch=4)

#legend
legend(x="topright", bty='n', xpd=NA, cex=1.0,
       c("scenario 1","scenario 2", "scenario 3", "scenario 4"), 
       col = c(col1, col2, col3, col4), lwd=3.0) 
        # gives the legend appropriate colors

legend(x="topright", bty='n', xpd=NA, cex=1.0,
       c("scenario 1","scenario 2", "scenario 3", "scenario 4"), 
       lwd=1.0, lty=0,
       pch=c(1,2,3,4)) 
        # gives the legend appropriate symbols 

#(b) White Fir

#x-axis
dispersalRates <- seq(0,7,0.1)

#Set up plot and plot lines
plot(dispersalRates, lambdaWhiteFirII, type='l', col=col2, 
     xlab='average dispersal distance (m)',
     ylab='critical speed of climate change (m/yr)', 
     lwd=3.0, cex.axis=0.8, xlim=c(0,7), ylim=c(0, 0.08)) 
lines(dispersalRates, lambdaWhiteFirIV, col=col4, lwd=3.0) 
lines(dispersalRates, lambdaWhiteFirIII, col=col3, lwd=3.0) 
lines(dispersalRates, lambdaWhiteFirI, col=col1, lwd=3.0) 

#Add in points
points(dispersalRates, lambdaWhiteFirI, pch=1)
points(dispersalRates, lambdaWhiteFirII, pch=2)
points(dispersalRates, lambdaWhiteFirIII, pch=3)
points(dispersalRates, lambdaWhiteFirIV, pch=4)

#legend
legend(x="topright", bty='n', xpd=NA, cex=1.0,
       c("scenario 1","scenario 2", "scenario 3", "scenario 4"), 
       col = c(col1, col2, col3, col4), lwd=3.0) 
        # gives the legend appropriate symbols 

legend(x="topright", bty='n', xpd=NA, cex=1.0,
       c("scenario 1","scenario 2", "scenario 3", "scenario 4"), 
       lwd=1.0, lty=0,
       pch=c(1,2,3,4)) 
        # gives the legend appropriate symbols 

dev.off()
setwd(main)


################################################################################
#Figures 4 & 5: Contributions from zombie forest <- <- <- <- <- <- <- <- <- <- <
################################################################################

# Is there any difference (in contributions from the back of the patch) when the 
# front of the patch has no survival/growth/reproduction?

#Sugar Maple First

#Set up K:
K = matrix(data=0, nrow=3, ncol=3)
K[1,3] = 5.0
#c = 0.5 
#K = 5.0 for highest critical rate of climate change
# (max c: I:0.7, II:1.1, III:0.7, IV:0.9)
#K = 10.0 for halfway through tested values from fig. 3 
# (max c: I:0.5, II:1.1, III:0.5, IV:0.7)

#Environmental scenario 1
data1 = runModel(subintervals = 2000,
                 timesteps = 50,
                 timestepsCC = 300,
                 c = 0.5,
                 AInside = sugarMaple,
                 ABefore = sugarMapleI,
                 AAfter = sugarMapleI,
                 K = K,
                 domainStart = -150,
                 domainLength = 500,
                 patchLength = 10,
                 accuracy = "accurate",
                 contributions = TRUE)

#Environmental scenario 2
data2 = runModel(subintervals = 2000,
                 timesteps = 50,
                 timestepsCC = 300,
                 c = 0.5,
                 AInside = sugarMaple,
                 ABefore = sugarMapleII,
                 AAfter = sugarMapleI,
                 K = K,
                 domainStart = -150,
                 domainLength = 500,
                 patchLength = 10,
                 accuracy = "accurate",
                 contributions = TRUE)

#Environmental scenario 3
data3 = runModel(subintervals = 2000,
                 timesteps = 50,
                 timestepsCC = 300,
                 c = 0.5,
                 AInside = sugarMaple,
                 ABefore = sugarMapleIII,
                 AAfter = sugarMapleI,
                 K = K,
                 domainStart = -150,
                 domainLength = 500,
                 patchLength = 10,
                 accuracy = "accurate",
                 contributions = TRUE)

#Environmental scenario 4
data4 = runModel(subintervals = 2000,
                 timesteps = 50,
                 timestepsCC = 300,
                 c = 0.5,
                 AInside = sugarMaple,
                 ABefore = sugarMapleIV,
                 AAfter = sugarMapleI,
                 K = K,
                 domainStart = -150,
                 domainLength = 500,
                 patchLength = 10,
                 accuracy = "accurate",
                 contributions = TRUE)

#Environmental scenario 1
summary1 = summarize(data1)
percent.table1 = density.by.descent(summary1)

#Environmental scenario 2
summary2 = summarize(data2)
percent.table2 = density.by.descent(summary2)

#Environmental scenario 3
summary3 = summarize(data3)
percent.table3 = density.by.descent(summary3)

#Environmental scenario 4
summary4 = summarize(data4)
percent.table4 = density.by.descent(summary4)

if (firstRun){
  sub_dir = "Data/Zombie Forest Contributions"
  if (!dir.exists(sub_dir)){
    dir.create(sub_dir)
  }
  setwd(file.path(main, sub_dir))
  
  saveRDS(percent.table1, file = "SMcontrib1.rds")
  saveRDS(percent.table2, file = "SMcontrib2.rds")
  saveRDS(percent.table3, file = "SMcontrib3.rds")
  saveRDS(percent.table4, file = "SMcontrib4.rds")
  
  saveRDS(data4, file = "SMexample.rds")
}

################################################################################
#White fir
#No survival/growth/recruitment at the front of the patch

K = matrix(data=0, nrow=5, ncol=5)
K[1,5] = 3.0 
#c = 0.01
#K = 3.0 for peak (max c: I:0.0, II:0.08, III:0.01, IV:0.02)
#K = 3.0 for halfway through tested values
#K (to avoid zeros 1.0) (max c: I:0.01. II:0.06, III:0.02, IV:0.025)

#Environmental scenario 1
wf.data1 = runModel(subintervals = 2000,
                    timesteps = 50,
                    timestepsCC = 300,
                    c = 0.02,
                    AInside = whiteFir,
                    ABefore = whiteFirI,
                    AAfter = whiteFirI,
                    K = K,
                    domainStart = -50,
                    domainLength = 200,
                    patchLength = 1,
                    accuracy = "accurate",
                    contributions = TRUE)

#Environmental scenario 2
wf.data2 = runModel(subintervals = 2000,
                    timesteps = 50,
                    timestepsCC = 300,
                    c = 0.02,
                    AInside = whiteFir,
                    ABefore = whiteFirII,
                    AAfter = whiteFirI,
                    K = K,
                    domainStart = -50,
                    domainLength = 200,
                    patchLength = 1,
                    accuracy = "accurate",
                    contributions = TRUE)

#Environmental scenario 3
wf.data3 = runModel(subintervals = 2000,
                    timesteps = 50,
                    timestepsCC = 300,
                    c = 0.02,
                    AInside = whiteFir,
                    ABefore = whiteFirIII,
                    AAfter = whiteFirI,
                    K = K,
                    domainStart = -50,
                    domainLength = 200,
                    patchLength = 1,
                    accuracy = "accurate",
                    contributions = TRUE)

#Environmental scenario 4
wf.data4 = runModel(subintervals = 2000,
                    timesteps = 50,
                    timestepsCC = 300,
                    c = 0.02,
                    AInside = whiteFir,
                    ABefore = whiteFirIV,
                    AAfter = whiteFirI,
                    K = K,
                    domainStart = -50,
                    domainLength = 200,
                    patchLength = 1,
                    accuracy = "accurate",
                    contributions = TRUE)

#Environmental scenario 1
summary1 = summarize(wf.data1)
wf.percent.table1 = density.by.descent(summary1)

#Environmental scenario 2
summary2 = summarize(wf.data2)
wf.percent.table2 = density.by.descent(summary2)

#Environmental scenario 3
summary3 = summarize(wf.data3)
wf.percent.table3 = density.by.descent(summary3)

#Environmental scenario 4
summary4 = summarize(wf.data4)
wf.percent.table4 = density.by.descent(summary4)

if (firstRun){
  sub_dir = "Data/Zombie Forest Contributions"
  if (!dir.exists(sub_dir)){
    dir.create(sub_dir)
  }
  setwd(file.path(main, sub_dir))
  
  saveRDS(wf.percent.table1, file = "WFcontrib1.rds")
  saveRDS(wf.percent.table2, file = "WFcontrib2.rds")
  saveRDS(wf.percent.table3, file = "WFcontrib3.rds")
  saveRDS(wf.percent.table4, file = "WFcontrib4.rds")
}

################################################################################
# Constructing plot

sub_dir = "Data/Zombie Forest Contributions"
if (!dir.exists(sub_dir)){
  dir.create(sub_dir)
}
setwd(file.path(main, sub_dir))
exampleSM = readRDS(file = "SMexample.rds")


sub_dir = "Figures"
if (!dir.exists(sub_dir)){
  dir.create(sub_dir)
}
setwd(file.path(main, sub_dir))

#Figure 4
svg(filename="fig4part1.svg", width=15, height=6, bg="transparent")

sz = 1.3

layout(mat = matrix(c(1, 2, 3, 3, 4, 4), 
                    nrow = 2, 
                    ncol = 3),
       heights = c(1, 1, 1),    # Heights of the two rows
       widths = c(1, 1, 1))     # Widths of the three columns

par(oma=c(0,0,0,0), #outer margin area (bot, left, top, right) (space around)
    mar=c(5,5,2,1))#margin area (bot, left, top, right))

# (a) Example of contributions before climate change

#Plot timestep 25

inMat = exampleSM[[1]]
outMat = exampleSM[[2]]
descendentMat = exampleSM[[3]]

inBefore = inMat[inMat[,1]==50,]
outBefore = outMat[outMat[,1]==50,]
descendentBefore = descendentMat[descendentMat[,1]==50,]

locations = inBefore[,2]
colors = viridis(3)

plot(locations, inBefore[,6], type="l", col=colors[3], lwd=2.0,
     xlab="location", ylab="population density", xlim=c(-30,75), cex=sz, 
     cex.axis=0.8*sz, cex.lab=1.1*sz)
lines(locations, outBefore[,6], col=colors[1], lwd=2.0)
lines(locations, descendentBefore[,6], col=colors[2], lwd=2.0)
abline(v=-5)
abline(v=5)

legend(x="top", bty='n', xpd=NA, cex=sz,
       c("zombies","descendants","others"), # puts text in the legend
       lwd=c(2.0,2.0,2.0),
       col=colors)

mtext('(a)', side=2, line=3, at=9800, cex=0.8*sz, las=1)

# (b) Example of contributions after climate change

#Plot timestep 150

inAfter = inMat[inMat[,1]==150,]
outAfter = outMat[outMat[,1]==150,]
descendentAfter = descendentMat[descendentMat[,1]==150,]

locations = inAfter[,2]
colors = viridis(3)

plot(locations, descendentAfter[,6], type="l", col=colors[2], lwd=2.0,
     xlab="location", ylab="population density", xlim=c(-30,75), cex=sz, 
     cex.axis=0.8*sz, cex.lab=1.1*sz)
lines(locations, inAfter[,6], col=colors[3], lwd=2.0)
lines(locations, outAfter[,6], col=colors[1], lwd=2.0)
abline(v=-5+100*0.5)
abline(v=5+100*0.5)

legend(x="top", bty='n', xpd=NA, cex=sz,
       c("zombies","descendants","others"), # puts text in the legend
       lwd=c(2.0,2.0,2.0),
       col=colors)

mtext('(b)', side=2, line=3, at=1.02e8, cex=0.8*sz, las=1)

dev.off()
setwd(main)

################################################################################
# Figure 5: Zombie forest contributions to the core population
################################################################################

sub_dir = "Data/Zombie Forest Contributions"
if (!dir.exists(sub_dir)){
  dir.create(sub_dir)
}
setwd(file.path(main, sub_dir))
percent.table1 = readRDS(file = "SMcontrib1.rds")
percent.table2 = readRDS(file = "SMcontrib2.rds")
percent.table3 = readRDS(file = "SMcontrib3.rds")
percent.table4 = readRDS(file = "SMcontrib4.rds")

sub_dir = "Figures"
if (!dir.exists(sub_dir)){
  dir.create(sub_dir)
}
setwd(file.path(main, sub_dir))
svg(file="fig5.svg", width=10, height=7, bg="transparent")

par(mfrow = c(2, 2)) #1 row, 2 columns
par(oma=c(1,1,0,1)) #outer margin area (bot, left, top, right) (space around)
par(mar=c(4,4,0,0)) #margin area (bot, left, top, right)

col1 = "#a6cee3" #light blue
col2 = "#33a02c" #dark green 
col3 = "#b2df8a" #light green
col4 = "#1f78b4" #dark blue

# (a) Percent originating from zombie forest Sugar Maple

sz = 1.3

times=percent.table1$timestep

#Calculate proportion of the non-zombie forest that is descended from zombies
percent.table1 = cbind(percent.table1, percent.table1$percent.descent/
                         (percent.table1$percent.descent + 
                            percent.table1$percent.in))
percent.table2 = cbind(percent.table2, percent.table2$percent.descent/
                         (percent.table2$percent.descent + 
                            percent.table2$percent.in))
percent.table3 = cbind(percent.table3, percent.table3$percent.descent/
                         (percent.table3$percent.descent + 
                            percent.table3$percent.in))
percent.table4 = cbind(percent.table4, percent.table4$percent.descent/
                         (percent.table4$percent.descent + 
                            percent.table4$percent.in))
colnames(percent.table1)[5] = "proportion.descent"
colnames(percent.table2)[5] = "proportion.descent"
colnames(percent.table3)[5] = "proportion.descent"
colnames(percent.table4)[5] = "proportion.descent"

plot(times, percent.table1$proportion.descent, type="l", lty=1.0, ylim=c(0,1.0), lwd=3.0,
     xlim=c(50,350), col=col1, xlab="", ylab="proportion", cex=sz, cex.lab=1.1*sz)
lines(times, percent.table2$proportion.descent, lty=1.0, lwd=3.0, col=col2)
lines(times, percent.table3$proportion.descent, lty=1.0, lwd=3.0, col=col3)
lines(times, percent.table4$proportion.descent, lty=1.0, lwd=3.0, col=col4)

points(times[seq(1, 350, 10)], percent.table1$proportion.descent[seq(1, 350, 10)], pch=1)
points(times[seq(1, 350, 10)], percent.table3$proportion.descent[seq(1, 350, 10)], pch=3)
points(times[seq(1, 350, 10)], percent.table4$proportion.descent[seq(1, 350, 10)], pch=4)
points(times[seq(1, 350, 10)], percent.table2$proportion.descent[seq(1, 350, 10)], pch=2)

#legend
legend(x="right", bty='n', xpd=NA, cex=1.0,
       c("scenario 1","scenario 2", "scenario 3", "scenario 4"), # puts text in the legend
       col = c(col1, col2, col3, col4), lwd=3.0) # gives the legend appropriate symbols 

legend(x="right", bty='n', xpd=NA, cex=1.0,
       c("scenario 1","scenario 2", "scenario 3", "scenario 4"), # puts text in the legend
       lwd=1.0, lty=0,
       pch=c(1,2,3,4)) # gives the legend appropriate symbols 


#(b) Percent zombie forest, Sugar maple

plot(times, percent.table1$percent.out, type="l", lty=1.0, lwd=3.0, col=col1, xlim=c(50,350),
     xlab="", ylab="", cex=sz, cex.lab=1.1*sz, ylim=c(0,1.0))
lines(times, percent.table2$percent.out, lty=1.0, lwd=3.0, col=col2)
lines(times, percent.table3$percent.out, lty=1.0, lwd=3.0, col=col3)
lines(times, percent.table4$percent.out, lty=1.0, lwd=3.0, col=col4)

points(times[seq(1, 350, 10)], percent.table1$percent.out[seq(1, 350, 10)], pch=1)
points(times[seq(1, 350, 10)], percent.table2$percent.out[seq(1, 350, 10)], pch=2)
points(times[seq(1, 350, 10)], percent.table3$percent.out[seq(1, 350, 10)], pch=3)
points(times[seq(1, 350, 10)], percent.table4$percent.out[seq(1, 350, 10)], pch=4)

#legend
legend(x="topright", bty='n', xpd=NA, cex=1.0,
       c("scenario 1","scenario 2", "scenario 3", "scenario 4"), # puts text in the legend
       col = c(col1, col2, col3, col4), lwd=3.0) # gives the legend appropriate symbols 

legend(x="topright", bty='n', xpd=NA, cex=1.0,
       c("scenario 1","scenario 2", "scenario 3", "scenario 4"), # puts text in the legend
       lwd=1.0, lty=0,
       pch=c(1,2,3,4)) # gives the legend appropriate symbols 

# (c) Percent originating from zombie forest White Fir

# Construct white fir plot

sub_dir = "Data/Zombie Forest Contributions"
if (!dir.exists(sub_dir)){
  dir.create(sub_dir)
}
setwd(file.path(main, sub_dir))

percent.table1 = readRDS(file = "WFcontrib1.rds")
percent.table2 = readRDS(file = "WFcontrib2.rds")
percent.table3 = readRDS(file = "WFcontrib3.rds")
percent.table4 = readRDS(file = "WFcontrib4.rds")

times = percent.table1$timestep

#Calculate proportion of the non-zombie forest that is descended from zombies
percent.table1 = cbind(percent.table1, percent.table1$percent.descent/
                         (percent.table1$percent.descent + 
                            percent.table1$percent.in))
percent.table2 = cbind(percent.table2, percent.table2$percent.descent/
                         (percent.table2$percent.descent + 
                            percent.table2$percent.in))
percent.table3 = cbind(percent.table3, percent.table3$percent.descent/
                         (percent.table3$percent.descent + 
                            percent.table3$percent.in))
percent.table4 = cbind(percent.table4, percent.table4$percent.descent/
                         (percent.table4$percent.descent + 
                            percent.table4$percent.in))
colnames(percent.table1)[5] = "proportion.descent"
colnames(percent.table2)[5] = "proportion.descent"
colnames(percent.table3)[5] = "proportion.descent"
colnames(percent.table4)[5] = "proportion.descent"

plot(times, percent.table1$proportion.descent, type="l", lty=1.0, ylim=c(0,1.0), lwd=3.0,
     xlim=c(50,350), col=col1,xlab="year", ylab="proportion", cex=sz, cex.lab=1.1*sz)
lines(times, percent.table2$proportion.descent, lty=1.0, lwd=3.0, col=col2)
lines(times, percent.table3$proportion.descent, lty=1.0, lwd=3.0, col=col3)
lines(times, percent.table4$proportion.descent, lty=1.0, lwd=3.0, col=col4)

points(times[seq(1, 350, 10)], percent.table2$proportion.descent[seq(1, 350, 10)], pch=2)
points(times[seq(1, 350, 10)], percent.table1$proportion.descent[seq(1, 350, 10)], pch=1)
points(times[seq(1, 350, 10)], percent.table3$proportion.descent[seq(1, 350, 10)], pch=3)
points(times[seq(1, 350, 10)], percent.table4$proportion.descent[seq(1, 350, 10)], pch=4)

legend(x="right", bty='n', xpd=NA, cex=1.0,
       c("scenario 1","scenario 2", "scenario 3", "scenario 4"), # puts text in the legend
       col = c(col1, col2, col3, col4), lwd=3.0) # gives the legend appropriate symbols 

legend(x="right", bty='n', xpd=NA, cex=1.0,
       c("scenario 1","scenario 2", "scenario 3", "scenario 4"), # puts text in the legend
       lwd=1.0, lty=0,
       pch=c(1,2,3,4)) # gives the legend appropriate symbols 


# (d) Percent zombie forest white fir
plot(times, percent.table1$percent.out, lwd=3.0, col=col1, type="l",
     lty=1.0, ylim=c(0,1.0), xlim=c(50,350), ylab="", xlab="year", cex=sz, cex.lab=1.1*sz)
lines(times, percent.table2$percent.out, lty=1.0, lwd=3.0, col=col2)
lines(times, percent.table3$percent.out, lty=1.0, lwd=3.0, col=col3)
lines(times, percent.table4$percent.out, lty=1.0, lwd=3.0, col=col4)

points(times[seq(1, 350, 10)], percent.table1$percent.out[seq(1, 350, 10)], pch=1)
points(times[seq(1, 350, 10)], percent.table2$percent.out[seq(1, 350, 10)], pch=2)
points(times[seq(1, 350, 10)], percent.table3$percent.out[seq(1, 350, 10)], pch=3)
points(times[seq(1, 350, 10)], percent.table4$percent.out[seq(1, 350, 10)], pch=4)

#legend
legend(x="right", bty='n', xpd=NA, cex=1.0,
       c("scenario 1","scenario 2", "scenario 3", "scenario 4"), # puts text in the legend
       col = c(col3, col2, col4, col1), lwd=3.0) # gives the legend appropriate symbols 

legend(x="right", bty='n', xpd=NA, cex=1.0,
       c("scenario 1","scenario 2", "scenario 3", "scenario 4"), # puts text in the legend
       lwd=1.0, lty=0,
       pch=c(1,2,3,4)) # gives the legend appropriate symbols 

sub_dir = "Figures"
if (!dir.exists(sub_dir)){
  dir.create(sub_dir)
}
setwd(file.path(main, sub_dir))
dev.off()





