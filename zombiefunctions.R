# zombiefunctions.R
# File to accompany zombieforests.R

#COMPADRE data can be downloaded from the Compadre website,
# https://compadre-db.org/
#This code requires COMPADRE_v.4.0.1.RData

library(gridGraphics)
grab_grob <- function(){
  grid.echo()
  grid.grab()
}

library(colorRamps) #For colorful plots
library(viridis) #For yellow/blue color ramps

################################################################################
## FUNCTION DEFINITIONS
################################################################################

#A function that plots the population density distribution.
#The population density distribution is plotted every 'n'th generation, where 
#the 'byGeneration' parameter sets 'n'.

#The size or age class plotted is specified by the parameter 'stage'. To plot 
#the total population density distribution summed across all stages, specify 
#stage as 'n+1', where 'n' is the total number of stages/size classes/ages 
#classes in the population.

makePlot <- function(data, byGeneration = 10, stage=4, title="", 
                     xlab="location", ylab="population density"){
  
  #Make initial plot, set up axes
  colors = viridis(max(data[,"time"])+1) 
  column = stage+2
  
  plotdata <- data[ which(data[,"time"] == 0), ]
  plot(plotdata[,"location"], plotdata[,column], type= "l", xlab = xlab, 
       ylab = ylab, main = title,
       xlim=c(min(data[,"location"]), max(data[, "location"])),
       ylim=c(min(data[,column]), max(data[,column])),
       col = colors[1])
  
  #Adds all the other distributions, one generation at a time
  for (i in 1:(as.integer(max(data[,"time"])/byGeneration))){
    plotdata <- data[ which(data[,"time"] == i*byGeneration), ]
    lines(plotdata[,"location"], plotdata[,column],
          col = colors[i*byGeneration +1])
  }
}

################################################################################

#This function plots the population density distribution over space for a single
#timestep of the model. The timestep is specified by the parameter 'genNum' and
#the stage/ size class/ age class is specified by the parameter 'stage', with 
#behavior identical to the 'makePlot' function.

plotStep <- function(data, genNum = 0, stage = 1, title=""){
  
  column = stage+2
  
  plotdata <- data[ which(data[,1] == genNum), ]
  plot(plotdata[,2], plotdata[,column], type= "l", xlab = "location", 
       ylab = "population density", main = title,
       xlim=c(min(data[,2]), max(data[, 2])),
       ylim=c(min(plotdata[,column]), max(plotdata[,column]))
  )
  
}

################################################################################

#This is just like plotstep, but it adds a line to the plot so that multiple
# distributions can be compared.

linesStep <- function(data, genNum = 0, stage = 1, title="", color="black"){
  
  column = stage+2
  
  plotdata <- data[ which(data[,"time"] == genNum), ]
  lines(plotdata[,"location"], plotdata[,column], type= "l", xlab = "location", 
        ylab = "population density", main = title,
        col = color,
        xlim=c(min(data[,"location"]), max(data[, "location"])),
        ylim=c(min(plotdata[,column]), max(plotdata[,column])),
        lwd=2.0, xlab="", ylab=""
  )
  
}

################################################################################

# A function that determines if given data structure represents a population 
# that is below the critical speed of climate change or above the critical speed 
# of climate change

#Find the maximum pop density across all generations (since we assume population 
#will grow for the first few generations, since there is no climate change).

#Find generation number at maximum pop density

#Boolean function:  YES/TRUE = population survives
# No/FALSE = population dies 

persist <- function(data){
  maxtimesteps = max(data[,1])
  return (maxtimesteps == data[which.max(data[,ncol(data)]), 1])
}
#Note that this function give misleading answers if the population grows for a 
#little while, then declines AND the model run doesn't last long enough to 
#capture the decline. This is why it is critical to investigate the dynamics and 
#run the model long enough to capture the dynamics, as I have done in the 
#analyses below

################################################################################

#A function that finds the critical rate of climate change, given subintervals, 
#AInside, ABefore, AAfter, and K, and precision interval, guess, domainStart, 
#domainLength, patchLength
#
# subintervals: this controls the spatial resolution of the population density 
#               data. Appropriate values are in the 1000-2000 range.
#
# AInside: the transition matrix inside of the patch
# ABefore: the transition matrix before (left of) the patch
# AAfter: the transition matrix after (right of) the patch
# K: matrix of dispersal kernels
# precision: desired precision of the critical rate of climate change
# guess: what do you think the critical rate of climate change could be?
# domainStart: position value of leftmost part of entire domain under 
# consideration
# domainLength: length of entire domain
# patchLength: length of the moving habitat patch, centered at zero


getCritClimRate <- function(subintervals, AInside, ABefore, AAfter, K, 
                            precision=0.001, guess=1.0,
                            domainStart, domainLength, patchLength, accuracy){
  
  #Test if guess is too big for domain
  if (-(domainStart)+(patchLength/2)+100*guess > domainLength){
    stop("The domainLength needs to increase to accomodate the large guess value.")
  }
  
  #Based on the standard bisection method
  b = guess
  a = 0.0
  dataA = runModel(subintervals = subintervals,
                   timesteps = 50,
                   timestepsCC = 100,
                   c = a, 
                   AInside = AInside,
                   ABefore = ABefore,
                   AAfter = AAfter,
                   K = K,
                   domainStart = domainStart,
                   domainLength = domainLength,
                   patchLength = patchLength,
                   accuracy = accuracy)
  
  dataB = runModel(subintervals = subintervals,
                   timesteps = 50,
                   timestepsCC = 100,
                   c = b, 
                   AInside = AInside,
                   ABefore = ABefore,
                   AAfter = AAfter,
                   K = K,
                   domainStart = domainStart,
                   domainLength = domainLength,
                   patchLength = patchLength,
                   accuracy = accuracy)
  
  if (!persist(dataA)){
    #stop("Population cannot persist, even with no climate change.")
    return (0.0)
  }
  if (persist(dataB)){
    stop("Guess needs to be higher.")
  }
  
  #Use bisection method to narrow in on the value
  while(b-a>precision){
    print(c) #Could remove print line later
    c = (a+b)/2
    dataC = runModel(subintervals = subintervals,
                     timesteps = 50,
                     timestepsCC = 100,
                     c = c, 
                     AInside = AInside,
                     ABefore = ABefore,
                     AAfter = AAfter,
                     K = K,
                     domainStart = domainStart,
                     domainLength = domainLength,
                     patchLength = patchLength,
                     accuracy = accuracy)
    if (persist(dataC) == persist(dataA)){
      a = c #check this
    } else {
      b = c
    }
  }
  return(c)
}

################################################################################
#Function to find multiplier to reduce eigenvalue

# Say you have a matrix and you want to know what scalar to multiply the 
# fecundity entries in the matrix by to change the dominant eigenvalue of that 
# matrix to something  new.
# This function finds that scalar multiplier, given the matrix and the desired
# dominant eigenvalue 'desiredLambda'. 
# Will find the scalar multiplier within the specified error tolerance 
# 'errorTol'.
# A similar function that operates on the survival/growth entries is defined 
# below.

#Uses bisection method

findMultiplierFecundity <- function(mat, desiredLambda, errorTol = 0.01){
  
  #Matrix size
  size = dim(mat)[1]
  a = 0.0;
  b = 1.0;
  c = 1.0;
  currentLambda = Re(eigen(mat)$values[1])
  
  #Test first
  vec = matrix(data=1, nrow=size, ncol=size)
  vec[1,2:size]=0
  newMat = vec*mat
  if(Re(eigen(newMat)$values[1]) > desiredLambda){
    stop("It's not possible to reduce lambda to the desired level")
  }
  
  while(abs(currentLambda-desiredLambda) > errorTol & c > 0){
    c = (a+b)/2
    #print(c)
    vec = matrix(data=1, nrow=size, ncol=size)
    vec[1,2:size]=c
    newMat = vec*mat
    currentLambda = Re(eigen(newMat)$values[1])
    if (currentLambda > desiredLambda){
      b = c 
    } else {
      a = c
    }
  }
  return(c)
}

################################################################################

# Finding survival multiplier

# This function is similarly to 'findMultiplierFecundity', but it finds a 
# multiplier by which to reduce all of the survival entries in the matrix, 
# instead of the fecundity entries.

findMultiplierSurvival <- function(mat, desiredLambda, errorTol = 0.01){
  
  #Matrix size
  size = dim(mat)[1]
  a = 0.0;
  b = 1.0;
  c = 1.0;
  currentLambda = Re(eigen(mat)$values[1])
  
  #Test first
  vec = matrix(data=0.0, nrow=size, ncol=size)
  vec[1,2:size] = 1.0
  newMat = vec*mat
  if(Re(eigen(newMat)$values[1]) > desiredLambda){
    stop("It's not possible to reduce lambda to the desired level")
  }
  
  #while (b-a > errorTol){
  while(abs(currentLambda-desiredLambda) > errorTol & c > 0){
    c = (a+b)/2
    #print(c)
    vec = matrix(data=c, nrow=size, ncol=size)
    vec[1,2:size]=1.0
    newMat = vec*mat
    currentLambda = Re(eigen(newMat)$values[1])
    if (currentLambda > desiredLambda){
      b = c 
    } else {
      a = c
    }
  }
  return(c)
}

################################################################################
# CONTRIBUTIONS FROM BACK OF THE PATCH (FUNCTION DEFINITIONS)
################################################################################

#This function summarizes the data by summing across size/stage classes to
#determine the total population density across all sizes/stages.

#Uses the standard data table produced by contrib (function that calculates
#contributions from the back of the patch) to create a summary table that 
#totals the population densities of individuals behind in patch, individuals
#behind the patch, and individuals descendent from individuals behind the patch
#at each timestep and spatial location.
summarize <- function(data){
  
  inMat = data[[1]]
  outMat = data[[2]]
  descendentMat = data[[3]]
  
  timestep = outMat[,1]
  location = outMat[,2]
  
  startClass = 3
  maxClass = ncol(outMat)
  
  pop.den.in = rowSums(inMat[,startClass:maxClass])
  pop.den.out = rowSums(outMat[,startClass:maxClass])
  pop.den.descendent = rowSums(descendentMat[,startClass:maxClass])
  pop.den.total = pop.den.in + pop.den.out + pop.den.descendent
  
  summary = data.frame(timestep,
                       location,
                       pop.den.in,
                       pop.den.out,
                       pop.den.descendent,
                       pop.den.total)
  
  return(summary)
}

################################################################################

#This function uses the summary data created from the 'summarize' function to 
#produce a plot that partitions total population density over space into 
#individuals in the patch, individuals behind the patch, and individuals 
#descendent from those behind the patch at the specified timesteps 'time'.

pop.den.contrib.plot <- function(summary, time, title=""){
  
  #Make initial plot (total population density), set up axes
  plotdata <- summary[ which(summary$"timestep"== time), ]
  plot(plotdata$"location", plotdata$"pop.den.total", type= "l", 
       xlab = "location", ylab = "population density", main = title)
  
  #Add in plot of descendents
  lines(plotdata$"location", plotdata$"pop.den.descendent", col = "blue")
  
  #Add in plot of individuals behind the patch
  lines(plotdata$"location", plotdata$"pop.den.out", col = "red")
  
  #Add in plot of individuals in the patch
  lines(plotdata$"location", plotdata$"pop.den.in", col="green")
}

################################################################################

#This plots the percent of individuals in the patch, behind the 
#patch, and descended from those behind the patch over time, using the 
#'percent.table' generated by the 'density.by.descent' table.

contributions.plot <- function(percent.table){
  
  total.contrib = percent.table$"percent.descent" + percent.table$"percent.out"
  
  #Make initial plot, set up axes
  plot(percent.table$"timestep", total.contrib, 
       xlab = "time", ylab="proportion of population",
       col="black", type='l')
  
  #add in other lines
  lines(percent.table$"timestep", percent.table$"percent.out", col="red")
  lines(percent.table$"timestep", percent.table$"percent.descent", col="blue")
}

################################################################################

#This function uses the 'summary' data generated by the 'summarize' function to
#make a table of percent of individuals in the patch, percent of individuals
#behind the patch, and percent of individuals descendent from those behind the 
#patch at each model timestep.

density.by.descent <- function(summary){
  steps = max(summary$"timestep")
  timestep = 0:steps
  percent.in = vector(length=steps)
  percent.out = vector(length=steps)
  percent.descent = vector(length=steps)
  for (i in timestep){
    subdata <- summary[which(summary$"time" == i), ]
    sum.total = sum(subdata$"pop.den.total")
    if (sum.total > 0){
      sum.in = sum(subdata$"pop.den.in")
      sum.out = sum(subdata$"pop.den.out")
      sum.descent = sum(subdata$"pop.den.descendent")
      percent.in[i+1] = sum.in/sum.total
      percent.out[i+1] = sum.out/sum.total
      percent.descent[i+1] = sum.descent/sum.total
    } else {
      percent.in[i+1] = 0
      percent.out[i+1] = 0
      percent.descent[i+1] = 0
    }
  }
  
  percent.table = data.frame(timestep,
                             percent.in,
                             percent.out,
                             percent.descent)
  
  return (percent.table)
  
}

################################################################################
# LOAD IN ALL DATA MATRICES

#Load in both matrices from Compadre
load(file = "COMPADRE_v.4.0.1.RData")

#Sugar Maple
sugarMaple = compadre$mat[[6413]]$matA

#White Fir
whiteFir = compadre$mat[[6360]]$matA

# Here is a list of the various scenarios outside the patch:
# I: Matrix is 0 outside of patch
# II: Matrix has no recruitment (but otherwise is not changed)
# III: reduction in fecundity to 0
# IV: Reduction in survival so that Lambda is 0.973 (white fir) 
#       or 0.81 (sugar maple)

#Find suitable reductions
################################################################################

#I: Matrix is 0 outside of patch

whiteFirI = matrix(data=0, nrow=5, ncol=5)
sugarMapleI = matrix(data=0, nrow=3, ncol=3)

################################################################################

#II: No recruitment 

whiteFirII = whiteFir
whiteFirII[2,1] = 0
Re(eigen(whiteFirII)$values[1]) #0.973

sugarMapleII = sugarMaple
sugarMapleII[2,1] = 0
sugarMapleII[3,1] = 0
Re(eigen(sugarMapleII)$values[1]) #0.81

################################################################################

#III: reduction in fecundity to 0

vec = matrix(data=1, nrow=3, ncol=3)
vec[1,2:3]=0

sugarMapleIII = vec*sugarMaple
Re(eigen(sugarMapleIII)$values[1]) #0.81

vec = matrix(data=1, nrow=5, ncol=5)
vec[1,2:5]=0

whiteFirIII = vec*whiteFir
Re(eigen(whiteFirIII)$values[1]) #0.973

################################################################################

#IV: Reduction in survival so that Lambda is 0.973 (white fir) 
#     or 0.81 (sugar maple)

findMultiplierFecundity(whiteFir, 0.973, errorTol=0.001) #0
findMultiplierSurvival(whiteFir, 0.973, errorTol=0.000000001) #0.8856162
findMultiplierFecundity(sugarMaple, 0.81, errorTol=0.001) #0
findMultiplierSurvival(sugarMaple, 0.81, errorTol=0.0000000001) #0.5696808

vec = matrix(data=0.8856162, nrow=5, ncol=5)
vec[1,2:5] = 1.0

whiteFirIV = whiteFir*vec
Re(eigen(whiteFirIV)$values[1]) #0.973

vec = matrix(data=0.5696808, nrow=3, ncol=3)
vec[1,2:3] = 1.0

sugarMapleIV = sugarMaple*vec
Re(eigen(sugarMapleIV)$values[1]) #0.81

