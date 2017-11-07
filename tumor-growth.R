# {α1, α2, β1, β2} = {1.07, 1.31, 1.47, 6.51}  / assuming cell size = 0.001 mm
set.seed(4)

vol = function(t, k){ 128 / (1 + ((128/.001)^.25 - 1) * exp(-.25*k*t))^4}

generate_K <- function(α1, α2) {
  #https://msalganik.wordpress.com/2017/01/21/making-sense-of-the-rlnorm-function-in-r/
  
  m <- α1
  s <- α2
  location <- log(m^2 / sqrt(s^2 + m^2))
  shape <- sqrt(log(1 + (s^2 / m^2)))
  print(paste("location:", location))
  print(paste("shape:", shape))
  k_Values <- rlnorm(n=1, location, shape)
  
  return(k_Values)
}

(4/3)*pi*(4)^3

generate_tumor <- function(α1, α2) {
  
  k <- generate_K(α1, α2)
    
  for (tumor_year in 0:20) {
    size <- vol(tumor_year, k)
    if (size >= 8)
      break;
  }
  
  if (tumor_year >= 20){
    generate_tumor(α1, α2) #search for another tumor
  }
  else
    return( c(tumor_year=tumor_year, size= size , vol(tumor_year+1, k), vol(tumor_year+2, k),vol(tumor_year+3, k), k=k))
    #return( c(α1=α1, α2=α2,tumor_year=tumor_year, k=k))

}

ks <- generate_K(1.07,1.31)

generate_tumor(1.07,1.31)





# *************

n <- 1000
ages <-sample(40:70, n, replace = TRUE)
pCA <- 0.02
gen_Ca <- function(age) { sample(c('Benign', 'CA'),size = 1, replace = TRUE, c(1 - pCA, pCA))}

apply_genTumor <- function(cancerStatus) {
    
  if (cancerStatus != 'CA')
    return (NA)
  else
    return (generate_tumor(1.07,1.31))
}

data = data.frame(ages)

data$BenignVsCA <- mapply(gen_Ca, ages)

data$Tumors<- mapply(apply_genTumor, data$BenignVsCA)

data$Tumor[[6]][1]

summary(data$Tumors)

summary(data$BenignVsCA)


results <- data[data$BenignVsCA == "CA",]

results



#lapply(results, function(x) write.table( data.frame(x), 'test.csv'  , append= T, sep=',' ))



write.table(
  data$BenignVsCA,
  na = "",
  file = "tumorSizes.csv",
  sep = ",",
  col.names = NA,
  qmethod = "double"
)





vol = function(t, k){ 128 / (1 + ((128/.001)^.25 - 1) * exp(-.25*k*t))^4}

volumes <-mapply(vol, t=3, k_Values)
mean(volumes)

plot(volumes, ylim=c(0,128))

mriDetect <- function(k, t_0, t_1){vol(3,k) < t_0 & vol(4,k) > t_1 }

filter <- mapply(mriDetect, k_Values, t_0 = 3, t_1 = 4)

plot(volumes[filter], ylim=c(0,128))






#check change in tumor after 1.7 years
ratio = function(k) {vol(2.7,k)/vol(1, k)}
ratio_Values <- mapply(ratio,k_Values)
mean(ratio_Values)








Vmax <- 128
cell <- .01

t <- 0
denominator <- (1 + ((Vmax/cell)^.25 - 1) * exp(-.25*k*t))^4

Vi_0 <- Vmax / denominator
Vi_0 

t <- 1.7
denominator <- (1 + ((Vmax/cell)^.25 - 1) * exp(-.25*k*t))^4

Vi_1 <- Vmax / denominator
Vi_1

#Overall, the mean time taken to grow from 10 mm to 20 mm was estimated as 1.7 years



#curve(eq, from=0, to=1.5, xlab="time", ylab="growth")

