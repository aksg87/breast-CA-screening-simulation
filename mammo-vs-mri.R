#library(plyr)

#clears workspace and setup
set.seed(4)

ls() 

#population size
n <- 10
ages <- sample(40:70,n,replace=TRUE) #uniform distribution
class <- sample(0:1,n,replace=TRUE) #0 is mammo 1 is MRI

mriBIRADS <- function(class){
  
  if (class== 0)
    return(NA);
  
  pBR1 <- .863 * (1/2)
  pBR2 <- .863 * (1/2)
  pBR3 <- .057
  pBR4 <- 0.081 * (1/2)
  pBR5 <- 0.081 * (1/2)

  eventProb <- sample(c('BR1','BR2','BR3', 'BR4', 'B53'),size=1, replace = TRUE, c(pBR1, pBR2, pBR3, pBR4, pBR5))

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
  else if (grade == 'BR3'){ ## CHECK THIS DATA POINT
    pCA <- 0
  }
  else if (grade == 'BR4'){
    print('BR4')
    pCA <- .02
  }
  else if (grade == 'BR5'){ ## CHECK THIS DATA POINT
    pCA <- .02
  }
  else {
    return (NA)
  }

  eventProb <- sample(c('Benign','CA'),size=1, replace = TRUE, c(1-pCA, pCA))
  
  return (eventProb);
}


mriStageCA <- function(x){
 
  pTis <- .38
  pT1a <- .08
  pT1b <- .25
  pT1c <- .23
  pT2 <-  .06
  
  eventProb <- sample(c('Tis','T1a','T1b','T1c','T2'),size=1, replace = TRUE, c(pTis, pT1a, pT1b, pT1c, pT2))
  
  return (eventProb);
}

data <- data.frame(ages,class)

data["mriBIRADS"] <- mapply(mriBIRADS, data$class)
data["mribenignVsCA"] <- mapply(mriBenignVsCA, data$mriBIRADS)
data["StageCA"] <- mapply(mriStageCA, data$mribenignVsCA)

data
