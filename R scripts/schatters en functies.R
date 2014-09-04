library(pim)
library(rje)

#mogelijke todo: nettere error afhandeling (gaat moeilijk wegens try ipv tryCatch in pim package)
#     in imputatie en dubbelrobust,
#reduced poset kan normaal sneller berekend worden


getInverse <- function(x, type = c("logit","probit")){
  #past de inverse van type toe op x
  #Arg:
  #   x: vector waarop de inverse functie toegepast moet worden
  #   type: naam van de functie  
  switch(type,
         logit = expit(x),
         probit = pnorm(x))
}

pseudoObs <- function(outcome) {
  #berekend de pseudoObservatie matrix van een gegeven vector
  #Arg: 
  #   outcome: de vector met observaties
  I <- t(matrix(rep(outcome,length(outcome)),ncol = length(outcome)))
  I <- ifelse(I < outcome,1, 0.5*(I==outcome) + 0)
  return(t(I))
}

differenceMatrix <-function(mat){
  #bereken het verschil (van rijen) van alle combinaties van 2 rijen,
  #de output geeft een matrix met nrow(mat)^2 rijen terug met als k-de rij
  #het verschill tussen de  (1 + (k-1 (mod nrow(mat)))ste rij
  #en  de  (ceiling((k)/ncol(mat)))ste rij
  #Arg: 
  #   mat: matrix of dataframe met de te vergelijken rijen  
  
  #indien de input een data.frame is wordt het naar een matrix gecast.
  if(class(mat) == "data.frame"){
    mat <- data.matrix(mat)
  }
  #maak een matrix aan bestande uit de originale matrix nrow(mat) keer onder elkaar 'geplakt'
  tempMat <- do.call("rbind",rep(list(mat),nrow(mat)))
  #maak een matrix waar iedere originele rij nrow keer herhaald wordt in dezelfde volgorde als in mat
  mat <- mat[rep(1:nrow(mat),rep(nrow(mat),nrow(mat))),] 
  return(tempMat-mat)
}

theoreticalProb <- function(alphaA,alphaL,sigmaE,sigmaL){
  #bereken de theoretische waarde van de marginale probabilistische index,
  #onder het model Y = alpha0 + alphaA * A + alphaL*L + E
  #met voor de verdelingen: (L,E) = (N(muL,sigmaL),N(0,sigmaE))
  #Arg: 
  #   zie hierboven
  return(pnorm(alphaA/(sqrt(2*(sigmaE^2 + alphaL^2*sigmaL^2)))))
}

reducedPoset <- function(data,treatment){
  #geeft een matrix van de poset vorm zoals die in het PIM packet gebruikt wordt
  #deze poset bevat enkel de elt'n waarvoor A_i=0 en A_j=1
  #Arg:
  #   data: dataframe waarvoor de poset gemaakt moet worden
  #   treatment: string die met de treatment correspondeert
  
  seq <- seq(nrow(data))
  posComb <- cbind(rep(seq,nrow(data)),rep(seq,each = nrow(data)))
  col1 <-posComb[,1][data[[treatment]][posComb[,1]]==0&data[[treatment]][posComb[,2]]==1]
  col2 <-posComb[,2][data[[treatment]][posComb[,1]]==0&data[[treatment]][posComb[,2]]==1]  
  return(cbind(col1,col2))
}

IPW <- function(outcome, treatment, confounders, dataframe,parametersPS,estimatePS = TRUE){
  #bereken de invers wegen schatter voor een marginale probabilistische index
  #waarbij de parameters voor de logistische regressie geschat kunnen worden
  #Arg: 
  #   dataframe: dataframe met alle benodigde waarden in
  #   outcome: String met de naam van de outcome in het dataframe
  #   treatment: String met de naam van de treatment in het dataframe
  #   confounders: (vector van) String(s) met de naam(en) van de confounders in het dataframe
  #   estimatePS: optionele optie voor indien false moeten er parameters meegegeven worden
  #   parametersPS: vector met als eerste element de 'intercept', en de andere parameters voor de PS
  #           in dezelfde volgorde als de confounders
  
  #bereken de fitted values van de logistische regressie
  if(estimatePS){
    glmFormula <- as.formula(paste(treatment,paste(confounders,collapse="+"),sep="~"))
    fit <- glm(glmFormula, family = "binomial",data=dataframe)
    fitted <- predict(fit, type = "response")
  }
  else{
    fitted <- as.vector((expit(parameters[1] + data.matrix(dataframe[confounders])%*%parameters[-1])))
  }
  #bereken twee vectoren met als waarden (1 -) treatment gedeeld door (1 -) fitted   
  factor1 <- (1-dataframe[[treatment]])/(1-fitted)
  factor2 <- dataframe[[treatment]]/fitted
  #maak een matrix van pseudoObservaties
  pseudoMat <- pseudoObs(dataframe[[outcome]])
  #zet de diagonaal elementen op nul zodat er uiteindelijk enkel over i != j gesommeerd wordt
  diag(pseudoMat) <- 0
  #bereken de vooropgestelde som en return het uiteindelijke resultaat
  sumValue <- as.numeric(t(factor1)%*%pseudoMat%*%factor2)
  return((sumValue/(nrow(dataframe)*(nrow(dataframe)-1))))
}

imputation <-function(outcome ,treatment ,confounders, dataframe, linkfunction, robust=FALSE){
  #bereken de imputatie schatter voor een marginale probabilistische index
  #Arg: 
  #   dataframe: dataframe waar de (numeric) waarden horende bij de formule terug te vinden zijn
  #   outcome: String met de naam van de outcome in het dataframe
  #   treatment: String met de naam van de treatment in het dataframe
  #   confounders: (vector van) String(s) met de naam(en) van de confounders in het dataframe
  #   robust: boolean, indien vals gebruikt hij de alle pseudo-observaties
  #           indien true enkel degene zoals gespecifieerd in reducedposet().
  
  #fit het PIM model en haal de coefficienten op
  formula <- as.formula(paste(outcome, paste(treatment, paste(confounders,collapse = " + "), sep= " + "), sep="~"))
  if(robust){
    fit <- pim(formula, link=linkfunction, data = dataframe, poset = reducedPoset(dataframe,treatment))
  }
  else{
    fit <- pim(formula, data = dataframe, poset = lexiposet, link=linkfunction)
  }
  if(is.null(fit$coeff)){
    return(NA)
  }
  coeff <- fit$coeff
  #maak een matrix met het verschil tussen de verschillende rijen van de matrix met confounders
  #en verwijder de ongewenste rijen met het verschil tussen dezelfde rij (dus = 0)
  diffMat <- differenceMatrix(dataframe[confounders])
  indicesToRemove <- (rep(nrow(dataframe),nrow(dataframe))+1)*c(0:(nrow(dataframe)-1))+1
  diffMat[indicesToRemove,] <- NA
  #fit P(Y <= Y'|A=0, A'=1,L,L') en bereken uiteindelijk de gewenste schatter
  fitted <- getInverse((coeff[1] - diffMat%*%coeff[-1]),linkfunction)
  return(sum(fitted, na.rm=TRUE)/(nrow(dataframe)*(nrow(dataframe)-1)))
}


doubleRobust <- function(outcome ,treatment ,confounders, dataframe,linkfunction,robust=FALSE
                             ,parametersPS,estimatePS = TRUE){
  #bereken de dubbel robuste schatter voor de marginale probabilistische index
  #Arg: 
  #   dataframe: dataframe waar de (numeric) waarden horende bij de formule terug te vinden zijn
  #   outcome: String met de naam van de outcome in het dataframe
  #   treatment: String met de naam van de treatment in het dataframe
  #   confounders: (vector van) String(s) met de naam(en) van de confounders in het dataframe
  #   robust: boolean, indien vals gebruikt hij de alle pseudo-observaties
  #           indien true enkel degene zoals gespecifieerd in reducedposet()
  #   parametersPS: vector met als eerste elemnt de 'intercept', en de andere parameters in de 
  #           zelfde volgorde als de confounders
  
  
  #bereken de fitted values van de logistische regressie
  if(estimatePS){
    glmFormula <- as.formula(paste(treatment,paste(confounders,collapse="+"),sep="~"))
    glmFit <- glm(glmFormula, family = "binomial",data = dataframe)
    glmFitted <- predict(glmFit, type = "response")
  }
  else{
    glmFitted <- as.vector((expit(parametersPS[1] + data.matrix(dataframe[confounders])%*%parametersPS[-1])))    
  }
  #bereken twee vectoren met als waarden (1 -) treatment gedeeld door (1 -) fitted (van logistisch model)  
  factor1 <- (1-dataframe[[treatment]])/(1-glmFitted)
  factor2 <-dataframe[[treatment]]/glmFitted
  #bereken het PIM model en bijhorende coefficienten
  pimFormula <- as.formula(paste(outcome, paste(treatment,paste(confounders,collapse="+"),sep="+"),sep="~"))
  if(robust){
    pimFit <- pim(pimFormula, link = linkfunction, data = dataframe, poset = reducedPoset(dataframe,treatment))
  }
  else{
    pimFit <- pim(pimFormula, data = dataframe, poset = lexiposet, link = linkfunction)
  }
  if(is.null(pimFit$coeff)){
    return(NA)
  }
  pimCoeff <- pimFit$coeff
  #bereken de difference matrix en zet de rijen horend bij L = L' op NA
  diffMat <- differenceMatrix(dataframe[confounders])
  indicesToRemove <- (rep(nrow(dataframe),nrow(dataframe))+1)*c(0:(nrow(dataframe)-1))+1
  diffMat[indicesToRemove,] <- NA
  #maak een matrix van fitted values aan die dezelfde vorm als de pseudomatrix heeft met NA op de diagonaal
  fittedMat <- matrix(getInverse(pimCoeff[1] - diffMat%*%pimCoeff[-1],linkfunction),nrow = nrow(dataframe))
  #bereken uiteindelijk de dubbel robuste schatter
  tempMat <- (pseudoObs(dataframe[[outcome]])-fittedMat)*factor1
  diag(tempMat) <- 0
  return((sum(tempMat%*%factor2,na.rm=T) + sum(fittedMat,na.rm=T))/((nrow(dataframe)*(nrow(dataframe)-1))))
}
