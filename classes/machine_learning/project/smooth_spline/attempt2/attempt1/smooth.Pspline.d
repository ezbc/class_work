.BG
.FN smooth.Pspline
.TL
Fit a Polynomial Smoothing Spline of Arbitrary Order
.DN
Returns an object of class `smooth.Pspline' which is a natural 
polynomial smooth of the input data of order fixed by the user.
.CS
smooth.Pspline(x, y, 
.OA
w=rep(1, length(x)), 
.OA
norder=2, 
.OA
df=norder + 2, 
.OA
spar=0, 
.OA
method=1)
.RA
.AG x
values of the predictor variable.  These must be strictly increasing, and
there must be at least `2*norder + 1' of them.
.AG y
one or more sets of response variable values.  If there is one response
variable, `y' is an array of the same length as `x'; if more than one, then
`y' is a matrix with `length(x)' rows and number of columns equal to the number
of variables.
.AG w
vector of positive weights for smoothing of the same length as `x'.  
If measurements at different values of x have different variances, 
`w' should be 
inversely proportional to the variances.  The default is that all weights
are one.
.AG norder
the order of the spline.  `norder = 2' gives the cubic smoothing spline, 
and more generally the smoothing function is a piecewise polynomial of 
degree `2*norder - 1'.  If derivatives are to be computed from the smoothing
using `predict.smooth.Pspline', the order should be one or two more than the
highest order of derivative.
.AG df
a number which specifies the degrees of freedom = trace(S).  
Here S is the implicit smoothing matrix.  `df' controls the amount of smoothing
if `method = 2'.
.AG spar
the usual smoothing parameter for smoothing splines,
which is the coefficient of the integrated squared derivative of order `norder'.
`spar' controls the amount of smoothing if `method = 1'.  
.AG method
the method for controlling the amount of smoothing.  `method = 1' uses the
value supplied for spar.  `method = 2' adjusts `spar' so that the degrees of
freedom is equal to `df'.  `method = 3' adjusts `spar' so that the 
generalized cross-validation criterion is minimized.  `method = 4' adjusts
`spar' so that the ordinary cross-validation criterion is minimized.
If 'method=3' or 'method=4', 'spar' defines the initial value for the 
minimization algorithm if positive; otherwise an internally generated value 
is used.
.RT
an object of class `smooth.Pspline' is returned, consisting of the fitted
smoothing spline evaluated at the supplied data, some fitting criteria
and constants.  This object contains the information necessary to evaluate
the smoothing spline or one of its derivatives at arbitrary argument
values using `predict.smooth.Pspline'.  The components of the returned
list are
.AG norder
the order of the spline  
.AG x
values of the predictor variable
.AG ysmth
a matrix with `length(x)' rows, each column of which contains
the smoothed response variable values for the corresponding column of `y'.  
.AG lev
leverage values, which are the diagonal elements of the smoother matrix S.
.AG gcv
generalized cross-validation criterion value
.AG cv
ordinary cross-validation criterion value
.AG df
a number which supplies the degrees of freedom = trace(S) rather than
a smoothing parameter.
.AG spar
the final smoothing parameter for smoothing splines.  This
is unchanged if `method = 1', but adjusted otherwise.
.AG call
the call that produced the fit
.DT
The method produces results similar to function `smooth.spline', but
the smoothing function is a natural smoothing spline rather than a B-spline
smooth, and as a consequence will differ slightly for `norder = 2' over the
initial and final intervals.  

The main extension is the possibility of 
setting the order of derivative to be penalized, so that derivatives
of any order can be computed using the companion function
`predict.smooth.Pspline'.  The algorithm is of order N, meaning that the 
number of floating point operations is proportional to the number of values
being smoothed.  Note that the argument values must be strictly increasing,
a condition that is not required by `smooth.spline'.  

Note that the appropriate or minimized value of the smoothing parameter
`spar' will depend heavily on the order; the larger the order, the smaller
this parameter will tend to be.
.SH REFERENCES
Heckman, N. and Ramsay, J. O. (1996) Spline smoothing with model based
penalties.  McGill University, unpublished manuscript.
.SA
`predict.smooth.Pspline, plot.smooth.Pspline, smooth.spline'
.EX
# order 2 smooth with fixed smoothing parameter:

fit <- smooth.Pspline(x, y, spar=1e-4)

# order 4 smooth by minimizing the gcv criterion:

fit <- smooth.Pspline(x, y, norder=4, method=3)

# order 3 smooth fixing degrees of freedom to 10

fit <- smooth.Pspline(x, y, norder=3, df= 10, method=2)

.KW smooth
.WR
