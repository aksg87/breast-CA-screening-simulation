set.seed(4)
library(plyr)

ls()

#population size
n <- 10000
ages <-
  sample(40:70, n, replace = TRUE) #uniform distribution
class <- sample(0:1, n, replace = TRUE) #0 is mammo arm & 1 is MRI arm

mriBIRADS <- function(class) {
  if (class == 0)
    return(NA)

  eventProb <- sample(
    c('BR3', 'BR45', 'BR12'),
    size = 1,
    replace = TRUE,
    c(.057, 0.063, (1 -.057 - 0.063) )
  )
  return (eventProb)
}


mriBenignVsCA <- function(grade) {
  pCA <- -1
  
  if (is.na(grade)) {
    return (NA)
  }
  
  else if (grade == 'BR3') {
    pCA <- 1 / 100
  }
  
  else if (grade == 'BR45') {
    pCA <- 35.8 / 100
  }
  
  else {
    pCA <- 0
  }
  
  eventProb <-
    sample(c('Benign', 'CA'),
           size = 1,
           replace = TRUE,
           c(1 - pCA, pCA))
  
  return (eventProb)
}


mriIncidentBIRADS <- function(grade, cancerStatus) {
  if (is.na(grade)) {
    return (NA)
  }
  
  else if ((!is.na(cancerStatus) &  cancerStatus == 'Benign')) {
    eventProb <-
      sample(
        c('BR3', 'BR45', 'BR12'), #DOUBLE LOGIC with BR3
        size = 1,
        replace = TRUE,
        c(.032, 0.023, (1 - 0.32 - 0.023))
      )
    
    return (eventProb)
  
  }
  
  return (NA)
}


mriIncidentBenignVsCA <- function(grade) {
  pCA <- -1
  
  if (is.na(grade)) {
    return (NA)
  }
  
  else if (grade == 'BR3') {
    pCA <- 1 / 100
    
  }
  else if (grade == 'BR45') {
    pCA <- 32.5 / 100
    
  }
  
  else {
    pCA <- 0
  }
  
  eventProb <-
    sample(c('Benign', 'CA'),
           size = 1,
           replace = TRUE,
           c(1 - pCA, pCA))
  
  return (eventProb)
}

data <- data.frame(ages, class)

data["MRI_Prev_BIRADS"] <- mapply(mriBIRADS, data$class)

data["MRI_BenignVsCA"] <-
  mapply(mriBenignVsCA, data$MRI_Prev_BIRADS)

data["MRI_Inc_BIRADS_0"] <-
  mapply(mriIncidentBIRADS, data$MRI_Prev_BIRADS, data$MRI_BenignVsCA)
data["MRI_Inc_BenignVsCA_0"] <-
  mapply(mriIncidentBenignVsCA, data$MRI_Inc_BIRADS_0)

for (year in 0:10) {
  BR <- paste(c("MRI_Inc_BIRADS", year), collapse = "_")
  B_vs_Ca <- paste(c("MRI_Inc_BenignVsCA", year), collapse = "_")
  
  BR_new <- paste(c("MRI_Inc_BIRADS", year + 1), collapse = "_")
  B_vs_Ca_new <-
    paste(c("MRI_Inc_BenignVsCA", year + 1), collapse = "_")
  
  data[, BR_new] <-
    mapply(mriIncidentBIRADS, data[, BR], data[, B_vs_Ca])
  data[,B_vs_Ca_new] <-mapply(mriIncidentBenignVsCA, data[,BR_new])
}

# data[paste(c("MRI_Inc_BIRADS", 0+1), collapse = "_")] <-mapply(mriIncidentBIRADS, data$paste(c("MRI_Inc_BenignVsCA", 0), collapse = "_"))
data[sapply(data, is.character)] <-
  lapply(data[sapply(data, is.character)], as.factor)

toCost = function(x){
  if (is.na(x)) {
    return (0)
  }
  
  else if (x == 'BR12') {
    return (500)
  }
  else if (x == 'BR3'){
    return (500+500)
  }
  else if (x == 'BR45') {
    return (500+1100)
  }
  else {
    return (0)
  }
}

selectCA = function(x){
  if (is.na(x)) {
    return (0)
  }
  
  else if (x == 'CA') {
    return (1)
  }
  else {
    return (0)
  }
}

dataAsDollars <- as.data.frame(lapply(data, FUN = function(x) {sapply(x, FUN = toCost)}))

selectCancers <-as.data.frame(lapply(data, FUN = function(x) { sum(sapply(x, FUN = selectCA) )  }))

dataSum <- as.data.frame(lapply(dataAsDollars, FUN = sum))

dataCount <- as.data.frame(sapply(data[], FUN = count))


write.table(
  data,
  na = "",
  file = "MRI_Intervals.csv",
  sep = ",",
  col.names = NA,
  qmethod = "double"
)


write.table(
  dataAsDollars,
  na = "",
  file = "MRI_Dollar_Conversion.csv",
  sep = ",",
  col.names = NA,
  qmethod = "double"
)

write.table(
  dataSum,
  na = "",
  file = "MRI_Dollar_Sum.csv",
  sep = ",",
  col.names = NA,
  qmethod = "double"
)

write.table(
  selectCancers,
  na = "",
  file = "MRI_Cancer_Count.csv",
  sep = ",",
  col.names = NA,
  qmethod = "double"
)

summary(data)

# Reference
# df <- as.data.frame(matrix(runif(5 * 5), 5, 5))
# convert <- function(x) {if ( x < 0.5) return (0) else return (1)}
# df2 <- as.data.frame(lapply(df[ ], FUN = function(x) {sapply(x, FUN = convert)}))
# ****************************************************************************************************************



mammoBIRADS <- function(class) {
  if (class == 1)
    return(NA)
  
  eventProb <- sample(
    c('BR3', 'BR45', 'REC', 'BR12'),
    size = 1,
    replace = TRUE,
    c(.023, 0.015, 0.185, (1 - 0.23 - 0.015 - 0.185))
    
    # RECALL – 10% ($188->467)         N and $
    #   BR3 – 2.3%    ($300->650)                     N and $
    #   BR4/5 – 1.5%   ($1000)                    N and $
    
  )
  return (eventProb)
}


mammoBenignVsCA <- function(grade) {
  pCA <- -1
  
  if (is.na(grade)) {
    return (NA)
  }
  
  else {
    return ('Benign')
  }
  
  # Assigning Benign to Allow MRI code to work
  # eventProb <-
  #   sample(c('Benign', 'CA'),
  #          size = 1,
  #          replace = TRUE,
  #          c(1 - pCA, pCA))
  # 
  # return (eventProb)
}


mammoIncidentBIRADS <- function(grade, cancerStatus) {
  if (is.na(grade)) {
    return (NA)
  }

  eventProb <- sample(
    c('BR3', 'BR45', 'REC', 'BR12'),
    size = 1,
    replace = TRUE,
    c(.023, 0.015, 0.085, (1 - 0.23 - 0.015 - 0.085))
    
    # RECALL – 7%  ($188)         N and $
    #   BR3 – 2.3%    ($300)                     N and $
    #   BR4/5 – 1.5%   ($1000)                    N and $
    
  )
  return (eventProb)

}


mammoIncidentBenignVsCA <- function(grade) {
  pCA <- -1
  
  if (is.na(grade)) {
    return (NA)
  }
  else {
    return ('Benign')
  }
  # Assigning Benign to Allow MRI code to work
}







data <- data.frame(ages, class)

data["MAMMO_Prev_BIRADS"] <- mapply(mammoBIRADS, data$class)

data["MAMMO_BenignVsCA"] <-
  mapply(mammoBenignVsCA, data$MAMMO_Prev_BIRADS)

data["MAMMO_Inc_BIRADS_0"] <-
  mapply(mammoIncidentBIRADS, data$MAMMO_Prev_BIRADS, data$MAMMO_BenignVsCA)
data["MAMMO_Inc_BenignVsCA_0"] <-
  mapply(mammoIncidentBenignVsCA, data$MAMMO_Inc_BIRADS_0)

year <- 0
for (year in 0:30) {
  BR <- paste(c("MAMMO_Inc_BIRADS", year), collapse = "_")
  B_vs_Ca <- paste(c("MAMMO_Inc_BenignVsCA", year), collapse = "_")
  
  BR_new <- paste(c("MAMMO_Inc_BIRADS", year + 1), collapse = "_")
  B_vs_Ca_new <-
    paste(c("MAMMO_Inc_BenignVsCA", year + 1), collapse = "_")
  
  data[, BR_new] <-
    mapply(mammoIncidentBIRADS, data[, BR], data[, B_vs_Ca])
  data[,B_vs_Ca_new] <-mapply(mammoIncidentBenignVsCA, data[,BR_new])
}

# data[paste(c("MAMMO_Inc_BIRADS", 0+1), collapse = "_")] <-mapply(mammoIncidentBIRADS, data$paste(c("MAMMO_Inc_BenignVsCA", 0), collapse = "_"))
data[sapply(data, is.character)] <-
  lapply(data[sapply(data, is.character)], as.factor)

toCost = function(x){
  if (is.na(x)) {
    return (0)
  }
  else if (x == 'BR12') {
    return (155)
  }
  
  else if (x == 'REC') {
    return (155+467)
  }
  else if (x == 'BR3'){
    return (155+650)
  }
  else if (x == 'BR45') {
    return (155 + 1000)
  }
  else {
    return (0)
  }
  
  # MAMMO ARM ($80)
  # FIRST
  # RECALL – 10% ($188)         N and $
  #   BR3 – 2.3%    ($300)                     N and $
  #   BR4/5 – 1.5%   ($1000)                    N and $
  #   
  #   SUBSEQUENT (12 month interval) x 11
  # RECALL – 7%  ($188)         N and $
  #   BR3 – 2.3%    ($300)                     N and $
  #   BR4/5 – 1.5%   ($1000)                    N and $
  
}

selectCA = function(x){
  if (is.na(x)) {
    return (0)
  }
  
  else if (x == 'CA') {
    return (1)
  }
  else {
    return (0)
  }
}

dataAsDollars <- as.data.frame(lapply(data, FUN = function(x) {sapply(x, FUN = toCost)}))

selectCancers <-as.data.frame(lapply(data, FUN = function(x) { sum(sapply(x, FUN = selectCA) )  }))

dataSum <- as.data.frame(lapply(dataAsDollars, FUN = sum))

dataCount <- as.data.frame(sapply(data[], FUN = count))


write.table(
  data,
  na = "",
  file = "Mammo_Intervals.csv",
  sep = ",",
  col.names = NA,
  qmethod = "double"
)


write.table(
  dataAsDollars,
  na = "",
  file = "Mammo_Dollar_Conversion.csv",
  sep = ",",
  col.names = NA,
  qmethod = "double"
)

write.table(
  dataSum,
  na = "",
  file = "Mammo_Dollar_Sum.csv",
  sep = ",",
  col.names = NA,
  qmethod = "double"
)

write.table(
  selectCancers,
  na = "",
  file = "Mammo_Cancer_Count.csv",
  sep = ",",
  col.names = NA,
  qmethod = "double"
)

summary(data)