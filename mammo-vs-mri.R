#library(plyr)

#clears workspace and setup
set.seed(4)

ls() 

#population size
n <- 10000
ages <- sample(40:70,n,replace=TRUE) #uniform distribution, will change later <Task>
class <- sample(0:1,n,replace=TRUE) #0 is mammo arm & 1 is MRI arm

mriBIRADS <- function(class){
  
  if (class== 0)
    return(NA);
  
  pBR1 <- .863 * (1/2)
  pBR2 <- .863 * (1/2)
  pBR3 <- .057
  
  
  #Turn Off for Testing to increase cancer rate
  pBR4 <- 0.081 * (1/2)
  pBR5 <- 0.081 * (1/2)
  
  #Turn On for Testing to increase cancer rate
  # pBR4 <- 0.9 * (1/2) 
  # pBR5 <- 0.9 * (1/2) 
  
  eventProb <- sample(c('BR1','BR2','BR3', 'BR4', 'BR5'),size=1, replace = TRUE, c(pBR1, pBR2, pBR3, pBR4, pBR5))

  return (eventProb);
}

mriBenignVsCA <- function(grade){
  pCA <- -1 
  
  if(is.na(grade)) {
    return (NA)
  }
  
  else if (grade == 'BR1' | grade == 'BR2'){
    pCA <- 0
  }
  else if (grade == 'BR3'){ ## DOUBLE CHECK THIS DATA POINT
    pCA <- 0
  }
  else if (grade == 'BR4' | grade == 'BR5'){
    
    #Turn Off for testing to increase cancer rate
    pCA <- 48/2120
    
    #Turn On for testing to increase cancer rate
    # pCA <- 1
    
  }
  else {
    return (NA)
  }

  eventProb <- sample(c('Benign','CA'),size=1, replace = TRUE, c(1-pCA, pCA))
  
  return (eventProb);
}


mriStageCA <- function(cancerStatus){
  
  if(is.na(cancerStatus)) {
    return (NA)
  }
  
  if (cancerStatus == "Benign"){
    return(NA);
  }
  
  pTis <- .38
  pT1a <- .08
  pT1b <- .25
  pT1c <- .23
  pT2 <-  .06
  
  eventProb <- sample(c('Tis','T1a','T1b','T1c','T2'),size=1, replace = TRUE, c(pTis, pT1a, pT1b, pT1c, pT2))
  
  return (eventProb);
}

mriIncidentBIRADS <- function(cancerStatus){
  
  if(is.na(cancerStatus)) {
    return (NA)
  }
  
  if (cancerStatus == 'CA') #Only Benign Patients Screened for Incidence
    return(NA);
  
  pBR1 <- .96 * (1/2)
  pBR2 <- .96 * (1/2)
  pBR3 <- .032
  
  #Turn Off for testing to increase cancer rate
  pBR4 <- 0.007 * (1/2)
  pBR5 <- 0.007 * (1/2)
  
  #Turn On for testing to increase cancer rate
  # pBR4 <- 0.7 * (1/2)
  # pBR5 <- 0.7 * (1/2)
  
  eventProb <- sample(c('BR1','BR2','BR3', 'BR4', 'BR5'),size=1, replace = TRUE, c(pBR1, pBR2, pBR3, pBR4, pBR5))
  
  return (eventProb);
}

mriIncidentBiopsy <- function(grade){
  
  if(is.na(grade) | grade == 'BR1' | grade == 'BR2' ) {
    return (NA)
  }
  
  #Turn Off for testing to increase biopsy rate
  pBiopsy = .333
  
  #Turn Onn for testing to increase biopsy rate
  # pBiopsy = 1
  
  eventProb <- sample(c('YesBiopsy','NoBiopsy'),size=1, replace = TRUE, c(pBiopsy, 1-pBiopsy))
  
  return (eventProb);
}

data <- data.frame(ages,class)

data["MRI_Prev_BIRADS"] <- mapply(mriBIRADS, data$class)
data["MRI_BenignVsCA"] <- mapply(mriBenignVsCA, data$MRI_Prev_BIRADS)
data["MRI_StageCA"] <- mapply(mriStageCA, data$MRI_BenignVsCA)

data["MRI_Inc_BIRADS"] <-mapply(mriIncidentBIRADS, data$MRI_BenignVsCA)
data["MRI_Inc_Biopsy"] <-mapply(mriIncidentBiopsy, data$MRI_Inc_BIRADS)

data[sapply(data, is.character)] <- lapply(data[sapply(data, is.character)], as.factor)

summary(data)

write.table(data, na = "", file = "mammo-MRI-screening.csv", sep = ",", col.names = NA, qmethod = "double")

