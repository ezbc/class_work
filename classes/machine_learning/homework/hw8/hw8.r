
x <- seq(1, 100)
y <- dnorm(x, m=10, sd=5) + 10*dnorm(x, m=90, sd=40) + 5*dnorm(x, m=33, sd=20)

y_noisey <- y + rnorm(100, m=0, sd=0.001)
y_noisey <- y_noisey / max(y_noisey)

library(pspline)

fit <- smooth.Pspline(x, y_noisey, method=3)
plot(x, y_noisey, xlab="x", ylab="y")
lines(fit$x, fit$ysmth)

R <- 1/length(x) * sum((y_noisey - fit$ysmth)^2)

# Now find optimal lambda value
lambdas = seq(0.1*fit$spar, 10*fit$spar, by=0.01*fit$spar)
Rs = array(c(100))
for (lambda in 1:length(lambdas)){
  fit <- smooth.Pspline(x, y_noisey, method=1, spar=lambdas[lambda])
  Rs[lambda] <- 1/length(x) * sum((y_noisey - fit$ysmth)^2)
}
R_optimal = min(Rs)

fit <- smooth.Pspline(x, y_noisey, method=1, spar=lambdas[which.min(Rs)])
plot(x, y_noisey, xlab="x", ylab="y")
lines(fit$x, fit$ysmth)

V = fit$spar
V
R/R_optimal



