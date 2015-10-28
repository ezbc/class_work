.BG
.FN predict.smooth.Pspline
.TL
Smoothing Spline of Arbitrary Order at New Data
.DN
Uses an object of class `smooth.Pspline' to evaluate a polynomial smoothing
spline of arbitrary order or one of its derivatives at new argument values.
.CS
predict.smooth.Pspline(splobj, xarg, 
.OA
nderiv=0)
.RA
.AG splobj
a fitted 'smooth.Pspline' object.
.AG xarg
the argument values at which the spline or its derivative is to be evaluated.
.AG deriv
the order of the derivative required -- default is 0, the function itself.
.RT
a list with components 'xarg' and 'dy'; the 'xarg' component is identical
to the input 'xarg' sequence, the 'dy' component is the evaluated derivative
of order 'deriv'. 
.DT
The method produces results similar to function `predict.smooth.spline', but
the smoothing function is a natural smoothing spline rather than a B-spline
smooth, and the order of the spline can be chosen freely, where order
in this case means the order of the derivative that is penalized.  
'smooth.spline' penalizes the second derivative, and consequently only
derivatives or order 0 or 1 are useful, but because 'smooth.Pspline'
penalizes a derivative of order m, derivatives up to order m-1 are useful.
The general recommendation is to penalize the derivative two beyond the
highest order derivative to be evaluated.

.SH REFERENCES
Heckman, N. and Ramsay, J. O. (1996) Spline smoothing with model based
penalties.  McGill University, unpublished manuscript.
.SA
`predict.smooth.Pspline, plot.smooth.Pspline, smooth.spline'
.EX
# order 4 smooth by minimizing the gcv criterion:

fit <- smooth.Pspline(x, y, norder=4, method=3)

# evaluate the second derivative or acceleration

D2fit <- predict.smooth.Pspline(fit, x, 2)

.KW smooth
.WR
