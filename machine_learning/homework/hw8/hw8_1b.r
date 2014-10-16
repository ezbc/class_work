
n_list <- c(32, 64, 128, 256, 512)
sigma_list <- c(0.1, 0.01)

results <- list()
R1_list <- c()
R2_list <- c()

for (i in 1:length(n_list)){
  n <- n_list[i]
  x <- seq(1, 100, length=n)
  y <- dnorm(x, m=10, sd=10) + 10*dnorm(x, m=90, sd=40) + 5*dnorm(x, m=33, sd=20)
  for (j in 1:2){
    I_results <- array(c(20))
    R_gcvs <- c()
    for (k in 1:20){
      y_noisey <- y + rnorm(n, m=0, sd=sigma_list[j])
      y_noisey <- y_noisey / max(y_noisey)
      
      library(pspline)
      
      fit <- smooth.Pspline(x, y_noisey, method=3)
      #plot(x, y_noisey, xlab="x", ylab="y")
      #lines(fit$x, fit$ysmth)
      
      R <- 1/length(x) * sum((y_noisey - fit$ysmth)^2)
      R_gcvs = c(R_gcvs, R)
      
      # Now find optimal lambda
      lambdas <- seq(0.1*fit$spar, 10*fit$spar, by=0.01*fit$spar)
      Rs <- array(c(1000))
      for (lambda in 1:length(lambdas)){
        fit <- smooth.Pspline(x, y_noisey, method=1, spar=lambdas[lambda])
        Rs[lambda] <- 1/length(x) * sum((y_noisey - fit$ysmth)^2)
      }
      R_optimal = min(Rs)
      
      # Calculate inefficiency
      I_results[k] <- R/R_optimal
    }
    results <- c(results, list(I_results, n_list[i], sigma_list[j]))
    if (j == 1){
      R1_list <- c(R1_list, mean(R_gcvs))
    } 
    if (j == 2){
      R2_list <- c(R2_list, mean(R_gcvs))
    }
  }
}

# Box plots
png(paste(c("hw8_1b_Rdist.png"), collapse="", sep=""), height=600, width=300)
par(mfrow=c(4,2))

for (i in 1:8) {
  n <- results[(i-1)*3 + 2][[1]]
  sigma <- results[(i-1)*3 + 3][[1]]
  name.pic <- paste("n = ", as.character(round(n, digits=0)), ", sigma = ", as.character(round(sigma*100, digits=2)), "%")
  boxplot(results[(i-1)*3 + 1], main=name.pic, ylab='I', asp=1)
}
graphics.off()

# R vs n
png(paste(c("hw8_1b_logR_vs_logn.png"), collapse="", sep=""), height=800, width=500)
plot(0.0001,0.0001,xlim = c(20,600),ylim = c(0.001,0.1), xlab="log(n)", ylab="log(R)", type = "n", log="xy")
cl <- c("red", "blue")
R_lists <- list(R1_list, R2_list)
for (i in 1:2){
  lines(n_list, R_lists[[i]], col=cl[i], type = 'b')
}
legend("topleft", legend = c("sigma = 1%", "sigma = 10%"), col=1:2, pch=1) # optional legend
graphics.off()

