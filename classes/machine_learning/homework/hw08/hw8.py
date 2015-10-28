#!/usr/bin/python

import numpy as np
def gauss(x, A, sigma, x0):
    y = A * np.exp(-(x - x0)**2 / (2 * sigma))
    return y

def main():

    ncomp = 3
    x = np.arange(0, 100, 1)
    y_true = np.zeros((100,))

    Alist = [0.5, 0.2, 1.0]
    sigmalist = [10, 20, 5]
    for comp in xrange(ncomp):
        x0 = np.random.randint(0, 100)
        y_true += gauss(x, Alist[comp], sigmalist[comp], x0)

    noise = 0.05 * np.random.normal(0,1,100)
    y_noisey = y_true + noise

    pspline(x=x, y=y_noisey, w=(1.0,)*len(x), method=3, norder=2)


if __name__ == '__main__':
    main()
