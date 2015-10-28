#!/usr/bin/python

#import numpy.f2py as f2py

#with open('Pspline.f', 'r') as f:
#    code = f.read()

#f2py.compile(code, modulename='Pspline')

from Pspline import pspline as fit_pspline

def pspline(x, y, w=None, norder=2, df=None, spar=0, method=1):

    import numpy as np

    if w is None:
        w = np.ones(len(x))
    if df is None:
        df = norder + 2

    n = len(x)
    if y.ndim == 2:
        nvar = y.shape[1] # Number of columns
    else:
        nvar = 1
    y = np.array(y)

    if len(w) == 1:
        w = w * np.ones(n)
    if y.shape[0] != n | len(w) != n:
        raise ValueError('Argument arrays of wrong length')
    if method != 1 | method !=2 | method != 3 | method != 4:
        raise ValueError('Wrong value for argument "method"')
    if norder <= 1 | norder >= 19:
        raise ValueError('Wrong value for argument "norder"')

    yhat = np.array((0, n, nvar, 0, 0))
    nworksiz = (n - norder) * (4 * norder + 3) + n
    work = np.zeros(nworksiz)
    lev = np.zeros(n)
    gcv = 0.0
    cv = 0.0
    dfmax = n * 1.0
    ier = 0.0
    irerun = 0
    lam = 0

    result = fit_pspline(norder, x, w, y, yhat, lev, gcv, cv, df, lam, dfmax,
            work, method, irerun, ier, n=n, nvar=nvar)

    return result

import numpy as np

x = np.arange(0,5,1)
y = x**2 + np.random.randint(5)

pspline(x, y, method=3)




