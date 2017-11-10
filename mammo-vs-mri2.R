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
    c(.057, 0.23, 0.92)
  )
  return (eventProb)
}


mriBenignVsCA <- function(grade) {
  pCA <- -1
  
  if (is.na(grade)) {
    return (NA)
  }
  
  else if (grade == 'BR45') {
    pCA <- 2.26 / 100
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
        c(.032, 0.023, 0.945)
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
  
  else if (grade == 'BR45') {
    pCA <- 0.69 / 100
    
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

for (year in 0:2) {
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
    return (626)
  }
  else if (x == 'BR3'){
    return (626+626)
  }
  else if (x == 'BR45') {
    return (626+1100)
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

selectCancers

table(data)

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
  file = "Dollar_Conversion.csv",
  sep = ",",
  col.names = NA,
  qmethod = "double"
)

write.table(
  dataSum,
  na = "",
  file = "Dollar_Sum.csv",
  sep = ",",
  col.names = NA,
  qmethod = "double"
)

write.table(
  selectCancers,
  na = "",
  file = "Cancer_Count.csv",
  sep = ",",
  col.names = NA,
  qmethod = "double"
)

summary(data)

# Reference
# df <- as.data.frame(matrix(runif(5 * 5), 5, 5))
# convert <- function(x) {if ( x < 0.5) return (0) else return (1)}
# df2 <- as.data.frame(lapply(df[ ], FUN = function(x) {sapply(x, FUN = convert)}))

