#{α1, α2, β1, β2} = {1.07, 1.31, 1.47, 6.51}
set.seed(4)

m <- 1.07
s <- 1.31
location <- log(m^2 / sqrt(s^2 + m^2))
shape <- sqrt(log(1 + (s^2 / m^2)))
print(paste("location:", location))
print(paste("shape:", shape))
k_Values <- rlnorm(n=10, location, shape)

k_Values

#https://msalganik.wordpress.com/2017/01/21/making-sense-of-the-rlnorm-function-in-r/
#one million draws from a log-normal distribution with a mean of 7 and a standard deviation of 75

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

