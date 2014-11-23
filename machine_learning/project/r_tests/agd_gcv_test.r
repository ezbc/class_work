### first you should check your working directory buy using the following command:###
getwd()
### then copy the hw6.dat into your working directory.###

####load your data into R.  Remember to put data in the working directory!!!##
data.hw6 <- read.table("hw6.dat", header=TRUE)
attach(data.hw6)

data = read.csv('spectrum0.csv',sep=" ", head=FALSE)
x = data$V1
y = data$V2

####now you can use x, y####

###LET'S START!###

### install the packages "pspline" ###
#install.packages("pspline")
### choose a minor, say 80 ###
### you should be able to see messages telling you R is installing the package ###

###load the package### (Make sure you install the package before loading) ###
library(pspline)

###Method 2:   GCV Method###

###We can use GCV method by setting the argument: method=3###

GCVfit <- smooth.Pspline(x, y, method=3, norder=6)
#GCVfit <- smooth.spline(x, y=ydata)

pred.prime <- predict(GCVfit, deriv=4)
yprime <- diff(GCVfit$ysmth[2:1958])/diff(GCVfit$x)


### to check out how the plot looks like, type in: ###
name.pic <- paste("lambda.GCV=",as.character(round(GCVfit$spar,digits=8)))
#plot(x,y,xlab="x",ylab="y",main=name.pic)
#lines(GCVfit$x, GCVfit$y)
#plot((0,100),(0,2),xlab="x",ylab="y",main=name.pic)
#lines(pred.prime$x, pred.prime$ysmth)

### then create a ps file for LaTeX: ###
postscript("GCVfit.ps", height=6, width=6, horizo=F)
name.pic <- paste("lambda.GCV=",as.character(round(GCVfit$spar,digits=8)))
#plot(x,y,xlab="x",ylab="y",main=name.pic)
#lines(GCVfit$x, GCVfit$y)
plot(GCVfit$x[2:1960], yprime, xlab="x",ylab="y",main=name.pic, type='l')
#lines(pred.prime$x, pred.prime$ysmth)

graphics.off()

###lambda by GCV###

GCVfit$spar

###FINISH!###

#This demo is made by Xiwen Ma. Edited by Zhigeng Geng. #
#Questions, comment?  Please send email to zgeng@stat.wisc.edu#
#HAVE FUN!#

