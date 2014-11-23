
#Bayesian confidence interval experiment generator.
library(gss)
library(tuneR)

sample_size <-32
x<-((1:sample_size)-.5)/sample_size
function_type = 'smooth'
y_true <- dnorm(x, m=0.5,sd=0.1)
#function_type = 'pointy'
#y_true <- dnorm(x, m=0.5,sd=0.01)/10
rep = 10 # 10 times

nlarep = seq (1: rep)
scorerep = seq(1:rep)
countrep = seq(1:rep)
sigma = seq(1:rep)*0.1
totalcount = 0

for (looprep in 1:rep) {
set.seed(1+looprep-1) #change the random seed if you want
ss <- .5  #true sigma for simulated noise 
#generates 100 data points on 0,1
y <- y_true + sigma[looprep] * rnorm(x)
#rnorm is a standard normal rv at each value of x
#can be multiplied by standard deviation other than 1
cubic.fit <- ssanova(y~x,type="cubic",method="v")
#method="v" is gcv est of lambda
est <- predict(cubic.fit,data.frame(x=x),se.fit=TRUE)

psfilename = paste(function_type, "_sig", sigma[looprep], "n", sample_size, ".png", sep = "")

# "postcript with filename sends everything to
#that file until graphics off
png(file= psfilename, height=400,width=400)
# height is inches centimeters?.
plot(x,y, xlim=c(0,1), ylim=c(-3.5,5.5))
# col=6 is color code 6 = cyan, 5 =purple 4=blue 
#3=green 2=red, 1=black (default)
lines(x,y_true,col=6)
#plots true fn, connected by lines
lines(x,est$fit,col=4)
#plots est fit in blue
lines(x,est$fit+1.96*est$se.fit,col=5)
lines(x,est$fit-1.96*est$se.fit,col=5)

#plots cconfidence interval in purple
############################
graphics.off()

trufun<-y_true
countrep[looprep] = sum(abs(trufun-est$fit)<1.96*est$se)
sigma[looprep]=summary(cubic.fit)$sigma
#nextone is mine
sigrat <- sigma/ss  #sigrat = ratio of estimated to "true" noise sd.
#counts no of trufn inside ci
totalcount = totalcount + countrep[looprep]

summary(cubic.fit)
#provides an estimate of the error standard deviation
nlarep[looprep] = cubic.fit$nla

#cubic.fig$nla gives log_10 lambda at the c.fit$score
scorerep[looprep] = cubic.fit$score
#cubic.fig$score gives value of GCV score at the lambda being used


}#end of loop

avecount =  totalcount/rep # average the count values


result = cbind (nlarep, scorerep, sigrat,countrep)
# a matrix of 10 rows and 3 columns for nla, score and count
# colnames(result) =c("nla", "gcvscore", "count") # give the column name

averow = cbind ('average', 'count','=', avecount)
result = rbind (result, averow)

namerow = cbind ('nla', 'gcvscore', 'sigrat','count')
result = rbind (namerow, result)
write (t(result), file = 'result.txt', ncolumns = 4)
#q()
