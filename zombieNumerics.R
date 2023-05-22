#zombieNumerics.R

#This program numerically solves an integrodifference equation model of a 
#habitat shifting in response to climate change, where the boundaries are still 
#sharp, but survival outside of the ideal patch is possible. The model includes 
#stage structure.

#Options for each model run:
#(1) Accuracy -- either "accurate" or "fast" (see notes below)
#(2) Demography over space -- the transition matrices inside the patch, to the
# left of the patch, and to the right of the patch can each be varied 
# independently.

#Global parameters
MU = 0.0; #Location of initial point release of the population
SIGMA = 2.0; #Square root of the variance of the initial population distribution


################################################################################
#Main Function
################################################################################

#Call runModel from the main code to get the data used to construct the plots.
### Parameters ###
#subintervals: increase for higher resolution; decrease for higher speed
#timesteps: these are the number of timesteps (or generations) prior to climate 
#           change, to allow the population to equilibriate
#timestepsCC: these are the number of timesteps with climate change
#c: the rate of climate change
#AInside: the transition matrix describing population growth inside the patch
#ABefore: the transition matrix describing population growth "behind" the patch
#AAfter: the transition matrix describing population growth "in front of" the 
#        patch
#K: matrix of dispersal kernel parameter values
#domainStart: spatial coordinate of the beginning of the domain
#domainLength: length the entire domain in which the patch will travel
#patchLength: length of the moving habitat patch
#accuracy: options are "accurate" and "fast" -- "accurate" is slower, but uses 
#          R's built-in integration routine to get an accurate result for the 
#          integration.
#          "fast" uses rectangles to approximate the integral, but is much 
#          faster than the "accurate" method. In practice, they both produce the 
#          same results given "enough" (~1000) subintervals.
#contributions: this parameter controls the type of data that is outputted. To 
#               tally the contributions of the zombie forest to individuals 
#               inside the patch, this boolean value should be set to TRUE. When 
#               contributions is set to TRUE, 3 data tables are outputted: one 
#               describing the population behind the patch, one describing the
#               population inside/in front of the patch and not descended from 
#               the back of the patch, and one describing the population inside/ 
#               in front of the patch and descended from the population behind 
#               the patch. When contributions is FALSE, these populations are
#               lumped together and reported in a single table. Setting 
#               contributions to TRUE slows the speed of computation.
runModel <- function(subintervals, timesteps, timestepsCC, c, AInside,
                     ABefore, AAfter, K, domainStart, domainLength, patchLength,
                     accuracy="accurate", contributions=FALSE){
  
  #Test that dimensions make sense
  check(AInside, ABefore, AAfter, K);
  
  #Test that the domain size makes sense
  if (domainStart > -(patchLength/2)){
    stop("Need smaller negative domainStart value.")
  }
  
  if (-(domainStart)+(patchLength/2)+timestepsCC*c > domainLength){
    stop("Need larger domainLength (or smaller c).")
  }
  
  #Patch start -- left spatial coordinate of patch
  patchStart = -patchLength/2;
  
  #Initialize main matrix
  #number of rows in data table
  numRows = subintervals*(timesteps + timestepsCC + 1); 
  numCols = ncol(AInside) + 2; #number of columns in data table
  
  #If we are counting contributions from the back of the patch
  if (contributions){ 
    #Set up the three blank matrices
    #Inside the patch, not descended from zombies
    inMatrix <- matrix(0, numRows, numCols); 
    outMatrix <- matrix(0, numRows, numCols); #Outside of the patch
    #Inside, descended from zombies
    descendentMatrix <- matrix(0, numRows, numCols); 
    
    #Set up the data tables with time/location/inital population values
    allMatrices = setInitialMatrices(subintervals, (timesteps+timestepsCC), 
                                     AInside.ncol(), inMatrix, outMatrix, 
                                     descendentMatrix, domainStart, 
                                     domainLength, patchStart);
    inMatrix = allMatrices[[1]]
    outMatrix = allMatrices[[2]]
    descendentMatrix = allMatrices[[3]]
    
  } else { #If we aren't counting contributions from the back of the patch
    mainMatrix <- matrix(0,numRows, numCols); #Fills matrix with zeros
    #set up initial matrix
    mainMatrix = setInitialMatrix(subintervals, timesteps + timestepsCC, 
                                  ncol(AInside), mainMatrix, domainStart, 
                                  domainLength);
  }
  
  #Go timestep by timestep, computing the dynamics
  for(t in 1:(timesteps+timestepsCC)){ #total number of timesteps
    
    #First, find the location of the traveling patch
    
    if (t <= timesteps){ #before climate change
      patchStartCurrent = patchStart; #left endpoint of patch
    } else { #after climate change
      patchStartCurrent = patchStart + (t-timesteps)*c; #left endpoint of patch
    }
    
    patchEnd = patchStartCurrent + patchLength; #right endpoint of patch
    
    #Then, calculate population density in each stage class
    if (accuracy == "accurate" || accuracy == "slow"){ #using integration
      #Uses R's built in numerical integration routines to solve the integral
      if (contributions){ #If measuring contributions from the back of the patch
        inMatrix = tryCatch({
          inMatrix = doStep_accurate(inMatrix, t, subintervals, AInside, 
                                     ABefore, AAfter, K, patchStartCurrent, 
                                     patchEnd, DomainLength);
        }, error = function(err) {
          inMatrix = doStep_fast(inMatrix, t, subintervals, AInside, ABefore, 
                                 AAfter, K, patchStartCurrent, patchEnd, 
                                 domainLength);
          return(inMatrix)
        })
        outMatrix = tryCatch({
          outMatrix = doStep_accurate(outMatrix, t, subintervals, AInside, 
                                      ABefore, AAfter, K, patchStartCurrent, 
                                      patchEnd, DomainLength);
        }, error = function(err){
          outMatrix = doStep_fast(outMatrix, t, subintervals, AInside, ABefore, 
                                  AAfter, K, patchStartCurrent, patchEnd, 
                                  domainLength);
          return(outMatrix)
        })
        descendentMatrix = tryCatch({
          descendentMatrix = doStep_accurate(descendentMatrix, t, subintervals, 
                                             AInside, ABefore, AAfter, K, 
                                             patchStartCurrent, patchEnd, 
                                             DomainLength);
        }, error = function(err){
          descendentMatrix = doStep_fast(descendentMatrix, t, subintervals, 
                                         AInside, ABefore, AAfter, K, 
                                         patchStartCurrent, patchEnd, 
                                         domainLength);
          return(descendentMatrix)
        })
      } else { #No contributions from the back
        mainMatrix = tryCatch({
          mainMatrix = doStep_accurate(mainMatrix, t, subintervals, AInside, 
                                       ABefore, AAfter, K, patchStartCurrent, 
                                       patchEnd, domainLength);
        }, error = function(err){
          mainMatrix = doStep_fast(mainMatrix, t, subintervals, AInside, 
                                   ABefore, AAfter, K, patchStartCurrent, 
                                   patchEnd, domainLength);
          return(mainMatrix)
        })
      }
    } else { #Fast, uses rectangles to integrate
      if (contributions){ #If measuring contributions from the back of the patch
        inMatrix = doStep_fast(inMatrix, t, subintervals, AInside, ABefore, 
                               AAfter, K, patchStartCurrent, patchEnd, 
                               domainLength);
        outMatrix = doStep_fast(outMatrix, t, subintervals, AInside, ABefore, 
                                AAfter, K, patchStartCurrent, patchEnd, 
                                domainLength);
        descendentMatrix = doStep_fast(descendentMatrix, t, subintervals, 
                                       AInside, ABefore, AAfter, K, 
                                       patchStartCurrent, patchEnd, 
                                       domainLength);
      } else { #No contributions from the back
        mainMatrix = doStep_fast(mainMatrix, t, subintervals, AInside, ABefore, 
                                 AAfter, K, patchStartCurrent, patchEnd, 
                                 domainLength);
        #Uses rectangles to integrate
      }
    }
    
    #Before ending this timestep, if we are measuring contributions from the 
    #back of the patch, we need to do transfers
    
    if (contributions){
      # (1) Transfer everything that was in the patch that landed out of the 
      #     patch to "outMatrix"
      output = transfer(outMatrix, inMatrix, patchStartCurrent, subintervals, 
                        t, 0);
      #0 == transer from in to out towards the left
      outMatrix = output[[1]]
      inMatrix = output[[2]]
      
      # (2) Transfer the seedlings from the back of the patch that landed in the 
      #     patch to the descendant class
      output = transfer(descendentMatrix, outMatrix, patchStartCurrent, 
                        subintervals, t, 1);
      #1 == transfer from outMatrix to descendentMatrix to the right
      descendentMatrix = output[[1]]
      outMatrix = output[[2]]
      
      # (3) Transfer everything that was in the "descendent" portion of the 
      #     patch that landed out of the patch to "outMatrix"
      output = transfer(outMatrix, descendentMatrix, patchStartCurrent, 
                        subintervals, t, 0);
      #0 == transfer from descendent to out towards the left
      outMatrix = output[[1]]
      descendentMatrix = output[[2]]
    }
    
  } #timestep by timestep loop
  
  #Process the matrix
  
  if(contributions){ #Measuring contributions from the back of the patch
    inMatrix = processMatrix(inMatrix)
    outMatrix = processMatrix(outMatrix)
    descendentMatrix = processMatrix(descendentMatrix)
    
  } else { #Not measuring contributions from the back of the patch
    mainMatrix = processMatrix(mainMatrix)
  }
  
  #Return the matrix
  #If we are measuring contributions from the back of the patch
  if (contributions){ 
    returnList = list(inMat = inMatrix, outMat = outMatrix, 
                      descendentMat = descendentMatrix)
    return(returnList);
  } else { #If we aren't measuring contributions from the back of the patch
    return (mainMatrix);
  }
  
}

#Add a column of total population density across all stage classes
#Add labels to the data
processMatrix <- function(matrix){
  numCols = ncol(matrix);
  
  #Total N column
  matrix = cbind(matrix, rowSums(matrix[,3:numCols]))
  
  matNames = list()
  matNames[1:2] = c("time", "location")
  for (pos in 3:numCols){
    matNames[pos] = paste("N", pos-2, sep="")
  }
  matNames[(numCols+1)] = "N"
  
  paste("N", 1, sep="")
  
  #Name columns
  colnames(matrix) <- matNames
  
  return(matrix)
}

#Sets up the initial data, assuming all individuals pooled (contributions FALSE)
setInitialMatrix <- function(subintervals, totalTimesteps, numClasses, 
                             mainMatrix, domainStart, domainLength){
  #Set up T (the first column)
  #Put the time values back in the matrix
  mainMatrix[,1] = rep(0:totalTimesteps, each=subintervals); 
  
  #Set up X (location) (the second column)
  for (i in (0:totalTimesteps)){
    for (j in (1:subintervals)) {
      mainMatrix[i*subintervals + j, 2] = domainStart + 
        (j-1+0.5)*(domainLength/subintervals); 
    }
  }
  
  #Set up N1 (the third column)
  mainMatrix[(1:subintervals),3] = initialN(mainMatrix[(1:subintervals),2]);
  #Only set population for first time step
  
  #All other stages start with zero, so they don't need to be set up
  
  #return the matrix
  return(mainMatrix)
}

#Sets up initial data tables for when contributions = TRUE
setInitialMatrices <- function(subintervals, totalTimesteps, numClasses, 
                               inMatrix, outMatrix, descendentMatrix, 
                               domainStart, domainLength, patchStart){
  
  #Set up the time column in each of the three matrices
  inMatrix[,1] = rep(0:totalTimesteps, each=subintervals)
  outMatrix[,1] = inMatrix[,1]
  descendentMatrix[,1] = inMatrix[,1]
  
  #Set up x (location) (the second column)
  for (i in (0:totalTimesteps)){
    for (j in (1:subintervals)){
      inMatrix[i*subintervals+j, 2] = domainStart + 
        (j-1+0.5)*(domainLength/subintervals);
    }
  }
  outMatrix[,2] = inMatrix[,2]
  descendentMatrix[,2] = inMatrix[,2]
  
  #Set up N1 (the third column) "seedlings"
  #Only set the population for timestep 0 (the first timestep)
  inMatrix[(1:subintervals),3] = initialN(inMatrix[(1:subintervals),2]);
  
  #Move anything outside the patch to the outside,
  #Transferring from donor to recipient
  lowerPatchBoundary = patchStart
  output = transfer(receivingMat=outMatrix, 
                    donorMat=inMatrix, lowerPatchBoundary, subintervals, 
                    timestep=0, 
                    direction=0); #to the left
  outMatrix = output[[1]]
  inMatrix = output[[2]]
  
  #Since this is the beginning of time, there aren't any descendants yet, 
  #so third column remains zeros in descendentMatrix
  
  #All other stages starte with zero, so they don't need to be set up
  
  matrices = list(inMatrix, outMatrix, descendentMatrix)
  return(matrices)
}

#Normally distributed initial population
initialN <- function(x){
  return (100/(sqrt(2*pi)*SIGMA)*exp(-(x-(MU))^2/(2*SIGMA^2)))
}

#Computes one timestep from the previous
#Uses rectangles to approximate integrals for speed
doStep_fast <- function(mainMatrix, timestep, subintervals, AInside, ABefore, 
                        AAfter, K, patchStart, patchEnd, domainLength){
  
  #Find the starting subinterval in the matrix and the ending subinterval in the 
  #matrix, for previous step
  #starting row position in mainMatrix for previous frame:
  startPrior = (timestep-1)*subintervals +1 ;  
  #ending row position in mainMatrix for previous frame:
  endPrior = timestep*subintervals; 
  #starting row position for current step:
  startCurrent = startPrior + subintervals; 
  #ending row position for current step:
  endCurrent = endPrior + subintervals; 
  
  #spatial values for entire domain:
  yvalues = mainMatrix[startPrior:endPrior,2] 
  #spatial values behind the patch:
  x1values = yvalues[which(yvalues < patchStart)] 
  #spatial values in the patch:
  x2values = yvalues[which(yvalues >= patchStart & yvalues <= patchEnd)]
  #spatial values in front of the patch:
  x3values = yvalues[which(yvalues > patchEnd)] 
  
  #Location of patch in the previous timestep:
  patchStartPrior = startPrior+length(x1values)-1
  patchEndPrior = endPrior-length(x3values)+1
  
  #Calculate distribution for each stage class one at a time
  numStages = ncol(mainMatrix) - 2; #number of stages or size classes
  
  #source stages in the plant transition matrix
  for (sourceStage in (1: numStages)){ 
    
    #destination stages in the plant transition matrix
    for (destStage in (1: numStages)){ 
      
      b = K[destStage, sourceStage]; #dispersal parameter
      
      #Previous population density values in area behind, in, and in front of 
      #the patch
      N1values =mainMatrix[startPrior:patchStartPrior, sourceStage+2]
      N2values =mainMatrix[(patchStartPrior+1):(patchEndPrior-1), sourceStage+2]
      N3values =mainMatrix[patchEndPrior:endPrior, sourceStage+2]
      
      if (b == 0){ #Case 1: no dispersal from sourceStage to destStage
        
        #Find current patch boundaries
        patchStartCurrent = patchStartPrior + subintervals
        patchEndCurrent = patchEndPrior + subintervals
        
        #behind the patch
        mainMatrix[startCurrent:patchStartCurrent, destStage+2] = 
          mainMatrix[startCurrent:patchStartCurrent, destStage+2] +
          ABefore[destStage, sourceStage]*N1values;
        
        #in the patch
        mainMatrix[(patchStartCurrent+1):(patchEndCurrent-1), destStage+2] =
          mainMatrix[(patchStartCurrent+1):(patchEndCurrent-1), destStage+2] +
          AInside[destStage, sourceStage]*N2values;
        
        #in front of the patch
        mainMatrix[patchEndCurrent: endCurrent, destStage+2] =
          mainMatrix[patchEndCurrent: endCurrent, destStage+2] +
          AAfter[destStage, sourceStage]*N3values;
        
      #if dispersal kernel is laplace distribution, have to perform integration
      } else { 
        
        #Here, we use simple rectangles to integrate
        #Fast, but not accurate
        
        KBefore = outer(yvalues, x1values, "laplace", b)
        KInside = outer(yvalues, x2values, "laplace", b)
        KAfter = outer(yvalues, x3values, "laplace", b)
        
        dx = domainLength/subintervals
        
        #Add on new area
        mainMatrix[startCurrent:endCurrent, destStage + 2] =
          mainMatrix[startCurrent:endCurrent, destStage + 2] + dx*
          ( KBefore %*% (ABefore[destStage, sourceStage]*N1values) +
              KInside %*% (AInside[destStage, sourceStage]*N2values) +
              KAfter %*% (AAfter[destStage, sourceStage]*N3values))
        
      } #End if/else
    } #End dest stage loop
  } #End source stage loop
  
  return(mainMatrix);
} #end doStep function

# Accurate version is slow
# Uses R-'s built-in integration routines to solve the integrodifference 
# equations. Parameters and output are equivalent to the fast version. The 
# calculation method is the only difference in this doStep function.
doStep_accurate <- function(mainMatrix, timestep, subintervals, AInside, 
                            ABefore, AAfter, K, patchStart, patchEnd, 
                            domainLength){
  
  #Find the starting subinterval in the matrix and the ending subinterval in the 
  #matrix, for previous step
  #starting row position in mainMatrix for previous frame
  startPrior = (timestep-1)*subintervals +1 ;  
  #ending row position in mainMatrix for previous frame
  endPrior = timestep*subintervals; 
  startCurrent = startPrior + subintervals;
  endCurrent = endPrior + subintervals;
  
  #All spatial values in domain:
  yvalues = mainMatrix[startPrior:endPrior,2] 
  #Behind the patch:
  x1values = yvalues[which(yvalues < patchStart)] 
  #In the patch:
  x2values = yvalues[which(yvalues >= patchStart & yvalues <= patchEnd)] 
  #In front of the patch:
  x3values = yvalues[which(yvalues > patchEnd)] 
  #x source
  
  #Patch boundaries in previous timestep
  patchStartPrior = startPrior+length(x1values)-1 #left patch endpoint
  patchEndPrior = endPrior-length(x3values)+1 #right patch endpoint
  
  #Calculate distribution for each stage class one at a time
  numStages = ncol(mainMatrix) - 2; #number of stages or size classes
  
  #Patch boundaries in current timestep
  patchStartCurrent = patchStartPrior + subintervals #Left patch endpoint
  patchEndCurrent = patchEndPrior + subintervals #right patch endpoint
  
  #source stages in the plant transition matrix
  for (sourceStage in (1: numStages)){ 
    
    #destination stages in the plant transition matrix
    for (destStage in (1: numStages)){ 
      
      #Population density in 3 sections of the domain
      #behind the patch:
      N1values = mainMatrix[startPrior:patchStartPrior, sourceStage+2] 
      #in the patch:
      N2values = mainMatrix[(patchStartPrior+1):(patchEndPrior-1), 
                            sourceStage+2]
      #in front of the patch:
      N3values = mainMatrix[patchEndPrior:endPrior, sourceStage+2] 
      
      b = K[destStage, sourceStage]; #dispersal parameter
      
      
      if (b == 0){ #Case 1: no dispersal from sourceStage to destStage
        
        #behind the patch
        mainMatrix[startCurrent:patchStartCurrent, destStage+2] = 
          mainMatrix[startCurrent:patchStartCurrent, destStage+2] +
          ABefore[destStage, sourceStage]*N1values;
        
        #in the patch
        mainMatrix[(patchStartCurrent+1):(patchEndCurrent-1), destStage+2] =
          mainMatrix[(patchStartCurrent+1):(patchEndCurrent-1), destStage+2] +
          AInside[destStage, sourceStage]*N2values;
        
        #in front of the patch
        mainMatrix[patchEndCurrent: endCurrent, destStage+2] =
          mainMatrix[patchEndCurrent: endCurrent, destStage+2] +
          AAfter[destStage, sourceStage]*N3values;
        
        #if dispersal kernel is laplace distribution, have to perform 
        #integration
      } else { 
        
        for (i in (startPrior:endPrior)){ #i tracks positions in the mainMatrix
          
          y = mainMatrix[i, 2]; #dest
          
          zeroTol = 0.001 #For integrating zero-functions
          
          
          integrand = function(x){ #dispersal
            d = function(x){
              return (laplace(y, x, b));
            }
            
            #growth
            f = splinefun(mainMatrix[startPrior:endPrior,2],
                          mainMatrix[startPrior:endPrior,(sourceStage+2)]);
            
            return(d(x)*f(x)) #dispersal*growth
          } #End integrand function
          
          #New population density values (this is the integrand)
          newN1values = integrand(x1values) #behind the patch
          newN2values = integrand(x2values) #in the patch
          newN3values = integrand(x3values) #in front of the patch
          
          #Built-in integration function has trouble with zero functions,
          #So we ensure it doesn't have to try to integrate any zero functions 
          #with some tests.
          
          #Case 1: Back of patch
          #as long as integrand is not zero everywhere
          if (!all(newN1values < zeroTol)){ 
            lowerIntegrationLimit = 
              x1values[min( which ( newN1values > zeroTol ))]
            back = ABefore[destStage, sourceStage]*
              integrate(integrand, lowerIntegrationLimit, patchStart)[[1]] 
              #takes just the value, not the error reports
          } else { #If integral is zero
            back = 0
          }
          
          #Case 2: Middle of patch
          #as long as integrand is not zero everywhere
          if (!all(newN2values < zeroTol)){  
            patch = AInside[destStage, sourceStage]*
              integrate(integrand, patchStart, patchEnd)[[1]]
          } else { #If integral is zero
            patch = 0
          }
          
          #Case 3: Front of patch
          #as long as integrand is not zero everywhere
          if (!all(newN3values < zeroTol)){ 
            upperIntegrationLimit = 
              x3values[max( which( newN3values > zeroTol ))]
            front = AAfter[destStage, sourceStage]*
              integrate(integrand, patchEnd, upperIntegrationLimit)[[1]]
          } else { #If integral is zero
            front = 0
          }
          
          
          #add the results of this integration to the total
          mainMatrix[i+subintervals, destStage+2] = 
            mainMatrix[i+subintervals, destStage+2] + back + patch + front
        }#end looping y-value by y-value
      } #end if/else statement for dispersal 0 or not 0
    } #destination stage loop
  } #source stage loop
  
  return(mainMatrix);
}#end doStep function

#Standard laplace distribution, used for dispersal
laplace = function(y,x,b){
  return ( (1/(2*b))*exp(-abs(x-y)/b) );
}

#Checks to ensure that dimensions of input match appropriately
check= function(AInside, ABefore, AAfter, K){
  
  if (nrow(ABefore) != ncol(ABefore)){
    stop("ABefore must be a square matrix.")
  }
  if (nrow(AAfter) != ncol(AAfter)){
    stop("AAfter must be a square matrix.")
  }
  if (nrow(AInside) != ncol(AInside)){
    stop("AInside must be a square matrix.")
  }
  if (nrow(AInside) != nrow(ABefore)){
    stop("AInside and ABefore must be the same size.")
  }
  if (nrow(AInside) != nrow(AAfter)){
    stop("AInside and AAfter must be the same size.")
  }
  if (nrow(K) != ncol(K)){
    stop("K should be a square matrix of dispersal parameters.")
  }
  if (nrow(K) != nrow(AInside)){
    stop("K must be the same size as AInside.")
  }
  
  
}


#Performs transfers between 3 matrices when computing contributions from the
#back of the patch
transfer = function(receivingMat, donorMat, lowerPatchBoundary, subintervals, 
                    timestep, direction){
  # bool direction from donor to recipient: 0 == left; 1 == right
  
  #Find number of stages
  numStages = ncol(donorMat)-2
  
  for (i in 1:subintervals){
    
    if ((!direction && receivingMat[i,2] < lowerPatchBoundary) 
        || (direction && receivingMat[i,2] >= lowerPatchBoundary)){ 
      #checks the direction
      
      for (stage in 1:numStages){ #transfers all size classes
        
        receivingMat[(timestep*subintervals+i), (stage+2)] <- 
          receivingMat[(timestep*subintervals+i), (stage+2)] + 
          donorMat[(timestep*subintervals+i), (stage+2)];
        donorMat[(timestep*subintervals+i), (stage+2)] = 0;
      } #End stage by stage
    } #if statement
  } #End subinterval by subinterval
  
  output = list(receivingMat, donorMat)
  
  return(output)
}

