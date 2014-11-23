### first you should check your working directory buy using the following command:###
getwd()
### then copy the hw6.dat into your working directory.###

####load your data into R.  Remember to put data in the working directory!!!##
data.hw6 <- read.table("hw6.dat", header=TRUE)
attach(data.hw6)

####now you can use x, y####

###LET'S START!###

### install the packages "pspline" ###
#install.packages("pspline")
### choose a minor, say 80 ###
### you should be able to see messages telling you R is installing the package ###

###load the package### (Make sure you install the package before loading) ###
#library(pspline)

###Method 1:  EYE BALL METHOD###

###Let's start from big lambda.###
###In the function of smooth.Pspline, spar corresponds to lambda###

fit1 <- smooth.Pspline(x, y, spar=0.001, method=1)

### to check out how the plot looks like, type in: ###
name.pic <- paste("lambda =",as.character(fit1$spar))
plot(x,y,xlab="x",ylab="y",main=name.pic)
lines(fit1$x, fit1$y)

### then create a ps file for LaTeX: ###
postscript("manual_fit.ps", height=6, width=6, horizo=F)
name.pic <- paste("lambda =",as.character(fit1$spar))
plot(x,y,xlab="x",ylab="y",main=name.pic)
lines(fit1$x, fit1$y)
graphics.off()


###Try smaller lambda###
fit2 <- smooth.Pspline(x, y, method=3)

### to check out how the plot looks like, type in: ###
name.pic <- paste("lambda=",as.character(fit2$spar))
plot(x,y,xlab="x",ylab="y",main=name.pic)
lines(fit2$x, fit2$y)

### then create a ps file for LaTeX: ###
postscript("gcv_fit.ps", height=6, width=6, horizo=F)
name.pic <- paste("lambda=",as.character(fit2$spar))
plot(x,y,xlab="x",ylab="y",main=name.pic)
lines(fit2$x, fit2$y)
graphics.off()

###repeat trying different lambda (0.01,0.001,0.0001,0.00001) see the difference###
###Don't forget to change name 'fitx' at each run ###
###my best pic is of lambda=0.002###

###Method 2:   GCV Method###

###We can use GCV method by setting the argument: method=3###

GCVfit <- smooth.Pspline(x, y, method=3)

### to check out how the plot looks like, type in: ###
name.pic <- paste("lambda.GCV=",as.character(round(GCVfit$spar,digits=8)))
plot(x,y,xlab="x",ylab="y",main=name.pic)
lines(GCVfit$x, GCVfit$y)

### then create a ps file for LaTeX: ###
postscript("GCVfit.ps", height=6, width=6, horizo=F)
name.pic <- paste("lambda.GCV=",as.character(round(GCVfit$spar,digits=8)))
plot(x,y,xlab="x",ylab="y",main=name.pic)
lines(GCVfit$x, GCVfit$y)
graphics.off()

###lambda by GCV###

GCVfit$spar

###FINISH!###


#This demo is made by Xiwen Ma. Edited by Zhigeng Geng. #
#Questions, comment?  Please send email to zgeng@stat.wisc.edu#
#HAVE FUN!#

