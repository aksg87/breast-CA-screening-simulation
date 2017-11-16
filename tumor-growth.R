# {α1, α2, β1, β2} = {1.07, 1.31, 1.47, 6.51}  / assuming cell size = 0.001 mm
set.seed(4)

intervals <- seq(1,20)

vol = function(t, k){ 128 / (1 + ((128/.001)^.25 - .001) * exp(-.25*k*t))^4}

generate_K <- function(α1, α2) {
  #https://msalganik.wordpress.com/2017/01/21/making-sense-of-the-rlnorm-function-in-r/
  
  m <- α1
  s <- α2 
  location <- log(m^2 / sqrt(s^2 + m^2))
  shape <- sqrt(log(1 + (s^2 / m^2)))
  k_Values <- rlnorm(n=1, location, shape)
  
  return(k_Values)
}


generate_tumor <- function(α1, α2, interval) {
  
  k <- generate_K(α1, α2)
  tumor_year <- sample(0:interval, 1, replace = FALSE)
  tumor_time_upper <-30
  
  while (TRUE) {
    size <- vol(tumor_year, k)
    if (size >= 5 | tumor_year >= tumor_time_upper)
      break;
    tumor_year <- tumor_year + 1
  }
  
  if (tumor_year >= tumor_time_upper){
    generate_tumor(α1, α2, interval) #search for another tumor
  }
  else
    return(size)
    #return( c(tumor_year=tumor_year, size= size , vol(tumor_year+1, k), vol(tumor_year+2, k),vol(tumor_year+3, k), k=k))
    #return( c(α1=α1, α2=α2,tumor_year=tumor_year, k=k))

}

pCA <- 0.02
gen_Ca <- function(age) { sample(c('Benign', 'CA'), size = 1, replace = TRUE, c(1 - pCA, pCA))}

apply_genTumor <- function(cancerStatus, interval, α1, α2) {
    
  if (cancerStatus != 'CA')
    return (NA)
  else
    return (generate_tumor(α1, α2,interval))
}


intervalApply <- function(interval, α1, α2, seed){
  
  print(interval)
  
  n <- 100000
  
  set.seed(seed)
  ages <-sample(40:70, n, replace = TRUE)
  data = data.frame(ages)
  data$BenignVsCA <- mapply(gen_Ca, ages)
  data$Tumors<- mapply(apply_genTumor, interval=interval, α1 = α1 , α2 = α2, data$BenignVsCA)
  m <- mean(data[data$BenignVsCA == "CA",]$Tumors)
  sd <- sd(data[data$BenignVsCA == "CA",]$Tumors)
  
  q99 <- quantile(data[data$BenignVsCA == "CA",]$Tumors, .99)
  q95 <- quantile(data[data$BenignVsCA == "CA",]$Tumors, .95)
  q90 <- quantile(data[data$BenignVsCA == "CA",]$Tumors, .90)
  q85 <- quantile(data[data$BenignVsCA == "CA",]$Tumors, .85)
  q80 <- quantile(data[data$BenignVsCA == "CA",]$Tumors, .80)
  q75 <- quantile(data[data$BenignVsCA == "CA",]$Tumors, .75)
  q50 <- quantile(data[data$BenignVsCA == "CA",]$Tumors, .50)
  q25 <- quantile(data[data$BenignVsCA == "CA",]$Tumors, .25)
  q10 <- quantile(data[data$BenignVsCA == "CA",]$Tumors, .10)
    
  
  return(list(mean=m, sd=sd, q99=q99, q95=q95, q90=q90,q85=q85,q80=q80,q75=q75, q50=q50, q25=q25, q10=q10))
}



generateData <- function(α1, α2, title){

  results <- mapply(intervalApply, intervals, α1 = α1 , α2 = α2, 5)
  
  results <- sapply(results, function(x)unlist(x))
  
  mean_volume <- results[c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)]
  sd_volume <- results[c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)]
  p99 <- results[c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)]
  p95 <- results[c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)]
  p90 <- results[c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)]
  p85 <- results[c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE)]
  p80 <- results[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE)]
  p75 <- results[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE)]
  p50 <- results[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE)]    
  p25 <- results[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE)]   
  p10 <- results[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE)] 
      
  png(filename=paste0("/Users/akshaygoel/Desktop/α1-", α1,"-α2-",α2,"-mean.png"), width = 4, height = 5, units = 'in', res = 150)
  plot(p99, cex = 0.30, ylim=c(0,50), col = "white", main="tumor volume percentiles", xlab="MRI interval length (years)", ylab="mean tumor volume (mm diameter)")
  grid()
  lines(p99, lty = 1)
  lines(p95, lty = 2)
  lines(p90, lty = 3)
  lines(p85, lty = 4)
  lines(p80, lty = 5)
  lines(p75, lty = 6)
  lines(p50, lty = 7)
  lines(p25, lty = 8)
  # lines(p10, lty = 9)
       
  mtext(paste0(title))

  text(intervals[c(TRUE, TRUE)], p99[c(TRUE, TRUE)]+.5, paste(round(p99[c(TRUE, TRUE)], digits = 1)), cex=0.40, font=4)
  text(intervals[c(TRUE, TRUE)], p95[c(TRUE, TRUE)]+.5, paste(round(p95[c(TRUE, TRUE)], digits = 1)), cex=0.40, font=4)
  text(intervals[c(TRUE, TRUE)], p90[c(TRUE, TRUE)]+.5, paste(round(p90[c(TRUE, TRUE)], digits = 1)), cex=0.40, font=4)
  text(intervals[c(TRUE, TRUE)], p85[c(TRUE, TRUE)]+.5, paste(round(p85[c(TRUE, TRUE)], digits = 1)), cex=0.40, font=4)
  text(intervals[c(TRUE, TRUE)], p80[c(TRUE, TRUE)]+.5, paste(round(p80[c(TRUE, TRUE)], digits = 1)), cex=0.40, font=4)
  text(intervals[c(TRUE, TRUE)], p75[c(TRUE, TRUE)]+.5, paste(round(p75[c(TRUE, TRUE)], digits = 1)), cex=0.40, font=4)
  text(intervals[c(TRUE, TRUE)], p50[c(TRUE, TRUE)]+.5, paste(round(p50[c(TRUE, TRUE)], digits = 1)), cex=0.40, font=4)
  text(intervals[c(TRUE, TRUE)], p25[c(TRUE, TRUE)]+.5, paste(round(p25[c(TRUE, TRUE)], digits = 1)), cex=0.40, font=4)

  legend(16, 25, legend=c("99th %", "95th %", "90th %", "85th %", "80th %", "75th %", "50th %", "25th %"), lty=c(1:8), cex=0.45)
  
  dev.off()
  
  # png(filename=paste0("/Users/akshaygoel/Desktop/α1-", α1,"-α2-",α2,"-sd.png"),, width = 4, height = 5, units = 'in', res = 150)
  # plot(sd_volume, ylim=c(0,50), xlim=c(0,20), cex = 0.30, main="σ of tumor volumes on detection",  xlab="MRI interval length (years)", ylab="σ tumor volume (mm diameter)")
  # mtext(paste0(title))
  # text(intervals, sd_volume + 1.2, paste(round(sd_volume, digits = 1)), cex=0.30)
  # dev.off()
  
  return(list(mean_volume,sd_volume))

}

x <- seq(1,10)

plot(x)
grid()

dev.off()

generateTumorGraph <- function(α1, α2, title){
  apply_gen = function(x){generate_tumor(α1, α2,1)}
  
  sample <-seq(1, 1000) 
  
  r <- mapply(apply_gen, sample)
  
  png(filename=paste0("/Users/akshaygoel/Desktop/tumors-visual-",α1,"-α2-",α2,".png"),, width = 8, height = 10, units = 'in', res = 300)
  plot(r, main="individual tumor volumes",  xlab="patient index", ylab="volume (mm diameter)")
  mtext(paste0(title, " (300 patients,  α1 ", α1, "  α2 ", α2,")"))
  dev.off()

}



all <- generateData(1.07, 1.31, "all patients")
young <- generateData(1.38, 1.36,  "age group 50 to 59")
old <- generateData(0.70, 1.18, "age group 60 to 69")







young_vs_old_mean <- unlist(young[1])-unlist(old[1])
young_vs_old_sd <- unlist(young[2])-unlist(old[2])



png(filename=paste0("/Users/akshaygoel/Desktop/young-mean-vs-old-mean.png"), width = 4, height = 5, units = 'in', res = 150)
plot(young_vs_old_mean, main="difference in mean tumor volume", ylim=c(0,50), xlim=c(0,20),  xlab="MRI interval of screening (years)", ylab="50-to-59 yrs - 60-to-69 yrs")
text(intervals, young_vs_old_mean + 1.2, paste(round(young_vs_old_mean, digits = 1)), cex=0.35)
dev.off()

png(filename=paste0("/Users/akshaygoel/Desktop/young-mean-vs-old-sd.png"), width = 4, height = 5, units = 'in', res = 150)
plot(young_vs_old_sd, main="difference in σ tumor volume", ylim=c(0,50), xlim=c(0,20),  xlab="MRI interval of screening (years)", ylab="50 to 59 yrs - 60-to-69 yrs")
text(intervals, young_vs_old_sd + 1.2, paste(round(young_vs_old_sd, digits = 1)), cex=0.35)
dev.off()


generateTumorGraph(1.07, 1.31, "all patients")
generateTumorGraph(1.38, 1.36, "age group 50 to 59")
generateTumorGraph(0.70, 1.18, "age group 60 to 69")


results <- mapply(intervalApply, intervals, α1 = α1 , α2 = α2, seed = 3)
results <- sapply(results, function(x) unlist(x))
mean_volume <- results[c(TRUE, FALSE, FALSE, FALSE)]
sd_volume <- results[c(FALSE, TRUE, FALSE, FALSE)]






plot(sd_volume, ylim=c(0,50), main="σ of tumor volumes on MRI detection",  xlab="MRI interval of screening (years)", ylab="sd tumor volume mm")


intervalApply3 <- function(interval, α1, α2){
  n <- 1000
  set.seed(4)
  ages <-sample(40:70, n, replace = TRUE)
  data = data.frame(ages)
  data$BenignVsCA <- mapply(gen_Ca, ages)
  data$Tumors<- mapply(apply_genTumor, interval=interval, data$BenignVsCA, α1, α2)
  return(data[data$BenignVsCA == "CA",]$Tumors)
}






# For all age groups combined, model parameters were estimated as 
# {α1, α2, β1, β2} = {1.07, 1.31, 1.47, 6.51}, while the two age 
# groups 50 to 59 years and 60 to 69 years gave estimates of
# {1.38, 1.36, 1.50, 6.33} and {0.70, 1.18, 1.46, 6.65}, respectively. 




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

