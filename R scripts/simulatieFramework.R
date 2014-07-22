#usage; clear de global environment en run schatter en functies.R voordat je deze file laat lopen

packages <- c('rje','pim','snowfall')
lapply(packages, require, character.only=T)

warningWrapper <- function(expr) {
  myWarnings <- 0
  wHandler <- function(w) {
    myWarnings <<-1
  }
  val <- withCallingHandlers(expr, warning = wHandler)
  c(val,myWarnings)
} 

wrapperFunction <- function(data,param){
    rbind(warningWrapper(IPW("y","a",c("l"),data,estimatePS=TRUE)),
          warningWrapper(imputation("y","a",c("l"),data,"probit")),
          warningWrapper(imputation("y","a",c("l"),data,"probit",robust=TRUE)),
          warningWrapper(doubleRobust("y","a",c("l"),data,"probit")),
          warningWrapper(doubleRobust("y","a",c("l"),data,"probit",robust=TRUE))
          #,warningWrapper(doubleRobust("y","a",c("l"),data,"probit",parametersPS = c(0,param["gamma"]),estimatePS=FALSE))
        )
}

sfInit(parallel=TRUE, cpus=2, type="SOCK")
sfExportAll()
sfLibrary(rje)
sfLibrary(pim)

simulatedData <-  function(param, nObs){
  l <- rnorm(nObs, mean = 0, sd = param["sigmL"])
  a <- rbinom(nObs, size = 1, prob = expit(param["gamma"]*l + param["gammaSq"]*l^2) )
  y <- rnorm(nObs, mean = param["alphaA"]*a + param["alphaL"]*l +param["alphaLSq"]*l^2, sd=param["sigmE"] )
  data <- data.frame(cbind(y,a,l))
  return(data)
}

simulatedDataList <- function(nruns, param, nObs){
   replicate(nruns,simulatedData(param,nObs),simplify = FALSE)
}

simulation <- function(param, nObs, nruns){
    tempDataList <- simulatedDataList(nruns, param, nObs)
    sfSapply(tempDataList,wrapperFunction, param=param)
}

summaryFunction <- function(matrix1,matrix2,theoreticalProb){
    I <- matrix(nrow = nrow(matrix1),ncol = 6)
    I[,1] <- rowMeans(matrix1 - theoreticalProb,na.rm=TRUE)
    I[,2] <- apply(matrix1,1,sd,na.rm=TRUE)
    I[,3] <- rowSums((matrix1 - theoreticalProb)^2,na.rm=TRUE)/rowSums(!is.na(matrix1))
    I[,4] <- apply(abs(matrix1-theoreticalProb),1,median,na.rm=TRUE)
    I[,5] <- rowSums(is.na(matrix1))
    I[,6] <- rowSums(matrix2)
    colnames(I) <- c("bias","MCSD","MSE","MAE","number of failures","non-convergence warnings")
    return(I)
}
  
fillFrame <-function(simFrame,nObs,nruns){
  ptm <- proc.time()
    for(i in 1:nrow(simFrame)){
      paramVect <- as.vector(t(simFrame[,-(ncol(simFrame)-2):-ncol(simFrame)][i,]))
      names(paramVect) <- names(simFrame)[-(ncol(simFrame)-2):-ncol(simFrame)]
      tempMatrix <- simulation(paramVect,nObs,nruns)
      resultMatrix <- tempMatrix[1:(nrow(tempMatrix)/2),]
      warningMatrix <- tempMatrix[(1+nrow(tempMatrix)/2):nrow(tempMatrix),]      
      simFrame$matList[[i]] <- tempMatrix
      simFrame$summaryList[[i]] <- summaryFunction(matrix1 = resultMatrix, matrix2 = warningMatrix,
                                                   theoreticalProb = simFrame$theoreticalProb[i])
      print(i)
      print(proc.time()-ptm)
    }
    simFrame
}

