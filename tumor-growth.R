#{α1, α2, β1, β2} = {1.07, 1.31, 1.47, 6.51}
set.seed(4)

m <- 1.07
s <- 1.31
location <- log(m^2 / sqrt(s^2 + m^2))
shape <- sqrt(log(1 + (s^2 / m^2)))
print(paste("location:", location))
print(paste("shape:", shape))
k_Values <- rlnorm(n=1000, location, shape)

k_Values

#https://msalganik.wordpress.com/2017/01/21/making-sense-of-the-rlnorm-function-in-r/
#one million draws from a log-normal distribution with a mean of 7 and a standard deviation of 75

vol = function(t, k){ 128 / (1 + ((128/.01)^.25 - 1) * exp(-.25*k*t))^4}

vol_afterTime = function(k) {vol(6,k)}
volumes <-mapply(vol_afterTime,k_Values)
mean(volumes)
plot(volumes)

#check change in tumor after 1.7 years
ratio = function(k) {vol(7.7,k)/vol(6, k)}
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

