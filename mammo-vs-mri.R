library(plyr)

#clears workspace and setup
set.seed(4)

ls() 

#population size
n <- 1000
ages <- sample(40:70,n,replace=TRUE) #uniform distribution
class <- sample(0:1,n,replace=TRUE) #0 is mammo 1 is MRI


