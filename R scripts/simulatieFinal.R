#beschouwde functies:
#invers wegen (geschatte param), imputatie(zowel lexico als reduced poset)
#doublerobust (geschatte param en zowel lexico als reduced) en (gegeven param en lexico)
nRuns <- 5000
nObs <- 100
paramNames <- c("alphaA","alphaL","sigmE","sigmL","gamma","alphaLSq","gammaSq","theoreticalProb")
simulatieFrame <- data.frame(t(rep(NA,length(paramNames))))
names(simulatieFrame) <- paramNames
simulatieFrame[1,] <- c(0,1,1,1,1,0,0,0.5)
simulatieFrame[2,] <- c(0,3,3,1,2,0,0,0.5)
simulatieFrame[3,] <- c(3,1,1,1,1,0,0,theoreticalProb(3,1,1,1))
simulatieFrame[4,] <- c(3,3,1,1,1,0,0,theoreticalProb(3,3,1,1))
simulatieFrame[5,] <- c(3,3,3,1,1,0,0,theoreticalProb(3,3,3,1))
simulatieFrame[6,] <- c(3,1,1,1,2,0,0,theoreticalProb(3,1,1,1))
simulatieFrame[7,] <- c(3,3,1,1,2,0,0,theoreticalProb(3,3,1,1))
simulatieFrame[8,] <- c(3,3,3,1,2,0,0,theoreticalProb(3,3,3,1))
simulatieFrame$matList <- vector("list",8)
simulatieFrame$summaryList <- vector("list",8)
set.seed(222)
simFrame <- fillFrame(simulatieFrame,nObs,nRuns)
save(simFrame,file="n100.Rda")
load("n100.Rda")

nRuns <- 1000
nObs <- 500
set.seed(111)
simFrameN500 <- fillFrame(simulatieFrame,nObs,nRuns)
save(simFrame,file="n500.Rda")



#beschouwde functies:
#invers wegen (geschatte param), imputatie(zowel lexico als reduced poset)
#doublerobust (geschatte param en zowel lexico als reduced)
nRuns <- 5000
nObs <- 100
paramNames <- c("alphaA","alphaL","sigmE","sigmL","gamma","alphaLSq","gammaSq","theoreticalProb")
simulatieFrame <- data.frame(t(rep(NA,length(paramNames))))
names(simulatieFrame) <- paramNames
simulatieFrame[1,] <- c(0,1,1,1,1,0,3,0.5)
simulatieFrame[2,] <- c(0,1,1,1,3,0,3,0.5)
simulatieFrame[3,] <- c(1,1,1,1,3,0,1,theoreticalProb(1,1,1,1))
simulatieFrame[4,] <- c(1,1,1,1,1,0,3,theoreticalProb(1,1,1,1))
simulatieFrame[5,] <- c(1,1,1,1,3,0,3,theoreticalProb(1,1,1,1))
simulatieFrame[6,] <- c(3,1,3,1,3,0,1,theoreticalProb(3,1,3,1))
simulatieFrame[7,] <- c(3,1,3,1,1,0,3,theoreticalProb(3,1,3,1))
simulatieFrame[8,] <- c(3,1,3,1,3,0,3,theoreticalProb(3,1,3,1))
simulatieFrame$matList <- vector("list",8)
simulatieFrame$summaryList <- vector("list",8)
set.seed(333)
simFrameN100AlphaLSQ <- fillFrame(simulatieFrame,nObs,nRuns)
save(simFrameN100AlphaLSQ,file="n100AlphaLSq.Rda")

nRuns <- 1000
nObs <- 500
set.seed(350)
simFrameN500GammaSQ <- fillFrame(simulatieFrame,nObs,nRuns)
save(simFrameN500GammaSQ,file="n500AlphaLSq.Rda")


nRuns <- 100
nObs <- 10000
set.seed(123456)
paramNames <- c("alphaA","alphaL","sigmE","sigmL","gamma","alphaLSq","gammaSq","theoreticalProb")
simulatieFrame <- data.frame(t(rep(NA,length(paramNames))))
names(simulatieFrame) <- paramNames
simulatieFrame[1,] <- c(0,1,1,1,1,3,0,0.5)
simulatieFrame[2,] <- c(0,3,1,1,1,3,0,0.5)
simulatieFrame[3,] <- c(1,3,1,1,1,1,0,0)
simulatieFrame[4,] <- c(1,3,1,1,1,3,0,0)
simulatieFrame[5,] <- c(1,1,1,1,1,3,0,0)
simulatieFrame[6,] <- c(3,3,3,1,1,1,0,0)
simulatieFrame[7,] <- c(3,3,3,1,1,3,0,0)
simulatieFrame[8,] <- c(3,1,3,1,1,3,0,0)
simulatieFrame$matList <- vector("list",8)
simulatieFrame$summaryList <- vector("list",8)

pnormWrapper <- function(l,alphaL,alphaLSq,alphaA,sigmE){
  pnorm((alphaL*(l[,2]-l[,1])+alphaLSq*(l[,2]^2-l[,1]^2)+alphaA)/(sqrt(2*sigmE^2)))
}

estimatedProbability<- list(vector("list",8),vector("list",8))

for(i in 1:nrow(simulatieFrame)){
  paramVect <- as.vector(t(simulatieFrame[,-(ncol(simulatieFrame)-2):-ncol(simulatieFrame)][i,]))
  names(paramVect) <- names(simulatieFrame)[-(ncol(simulatieFrame)-2):-ncol(simulatieFrame)]
  x <- replicate(nRuns,matrix(rnorm(2*nObs, mean = 0, sd = paramVect["sigmL"]),ncol=2),simplify = FALSE)
  results <- colMeans(sapply(x,pnormWrapper,alphaL=paramVect["alphaL"],alphaLSq=paramVect["alphaLSq"]
                             ,alphaA=paramVect["alphaA"],sigmE=paramVect["sigmE"]))
  simulatieFrame$theoreticalProb[i] <- mean(results)
  summaryVector <- c(mean(results),sd(results))
  names(summaryVector) <- c("mean","sd")
  estimatedProbability[[1]][[i]] <- results
  estimatedProbability[[2]][[i]] <- summaryVector
}
simulatieFrame$theoreticalProb[1:2]  <- 0.5
nRuns <- 5000
nObs <- 100
set.seed(444)
simFrameN100AlphaLSQ <- fillFrame(simulatieFrame,nObs,nRuns)
save(simFrameN100AlphaLSQ,file="n100AlphaLSq2.Rda")
save(estimatedProbability,file="estProbnAlphaLSq.Rda")

nRuns <- 1000
nObs <- 500
set.seed(555)
simFrameN500AlphaLSQ <- fillFrame(simulatieFrame,nObs,nRuns)
save(simFrameN500AlphaLSQ,file="n500AlphaLSqtemp.Rda")
save(simFrameN500AlphaLSQ,file="n500AlphaLSq2.Rda")

nRuns <- 100
nObs <- 10000
set.seed(123)
paramNames <- c("alphaA","alphaL","sigmE","sigmL","gamma","alphaLSq","gammaSq","theoreticalProb")
simulatieFrame <- data.frame(t(rep(NA,length(paramNames))))
names(simulatieFrame) <- paramNames
simulatieFrame[1,] <- c(0,3,3,1,3,3,3,0)
simulatieFrame[2,] <- c(0,3,3,1,3,3,1,0)
simulatieFrame[3,] <- c(0,3,3,1,3,1,3,0)
simulatieFrame[4,] <- c(0,3,3,1,3,1,1,0)
simulatieFrame[5,] <- c(3,3,3,1,3,3,3,0)
simulatieFrame[6,] <- c(3,3,3,1,3,3,1,0)
simulatieFrame[7,] <- c(3,3,3,1,3,1,3,0)
simulatieFrame[8,] <- c(3,3,3,1,3,1,1,0)


simulatieFrame$matList <- vector("list",8)
simulatieFrame$summaryList <- vector("list",8)

estimatedProbability2<- list(vector("list",8),vector("list",8))

for(i in 1:nrow(simulatieFrame)){
  paramVect <- as.vector(t(simulatieFrame[,-(ncol(simulatieFrame)-2):-ncol(simulatieFrame)][i,]))
  names(paramVect) <- names(simulatieFrame)[-(ncol(simulatieFrame)-2):-ncol(simulatieFrame)]
  x <- replicate(nRuns,matrix(rnorm(2*nObs, mean = 0, sd = paramVect["sigmL"]),ncol=2),simplify = FALSE)
  results <- colMeans(sapply(x,pnormWrapper,alphaL=paramVect["alphaL"],alphaLSq=paramVect["alphaLSq"]
                             ,alphaA=paramVect["alphaA"],sigmE=paramVect["sigmE"]))
  simulatieFrame$theoreticalProb[i] <- mean(results)
  summaryVector <- c(mean(results),sd(results))
  names(summaryVector) <- c("mean","sd")
  estimatedProbability2[[1]][[i]] <- results
  estimatedProbability2[[2]][[i]] <- summaryVector
}
simulatieFrame$theoreticalProb <- as.vector(cbind(c(0.5,0.5,0.5,0.5),
    as.vector(cbind(rep(estimatedProbability2[[2]][[5]][1],2),rep(estimatedProbability2[[2]][[7]][1],2)))))

nRuns <- 5000
nObs <- 100
set.seed(450)
simFrameN100AlphaLSQgammaSQ <- fillFrame(simulatieFrame,nObs,nRuns)
save(simFrameN100AlphaLSQgammaSQ,file="n100AlphaLSqGammaSq.Rda")
save(estimatedProbability2,file="estProbnAlphaLSqGammaSq.Rda")

nRuns <- 1000
nObs <- 500
set.seed(666)
simFrameN500AlphaLSQgammaSQ <- fillFrame(simulatieFrame,nObs,nRuns)
save(simFrameN500AlphaLSQgammaSQ,file="n500AlphaLSqGammaSq.Rda")


