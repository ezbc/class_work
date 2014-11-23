#!/usr/bin/python

''' Script to test generalized cross validation on same test set as in Lindner
et al. (2014)
http://adsabs.harvard.edu/abs/2014arXiv1409.2840L

'''

def gauss(x, sigma, x0, A):

    import numpy as np

    return A * np.exp(-(x - x0)**2 / (2.0 * sigma**2))

def gauss_1st_deriv(x, sigma, x0, A):

    ''' Returns the first derivative of a Gaussian. See gauss_4thderiv.nb
    mathematica code.

    '''

    import numpy as np

    expo = np.exp(-(x - x0)**2 / (2 * sigma**2))

    result = - A * expo * (x - x0) / sigma**2

    return result

def gauss_4th_deriv(x, sigma, x0, A):

    ''' Returns the fourth derivative of a Gaussian. See gauss_4thderiv.nb
    mathematica code.

    '''

    import numpy as np

    expo = np.exp(-(x - x0)**2 / (2 * sigma**2))

    result = A * expo * (x - x0)**4 / sigma**8 - \
             6 * A * expo * (x - x0)**2 / sigma**6 + \
             3 * A * expo / sigma**4

    return result

def test_prep_spectrum():

    import pickle
    import pspline

    test_data = pickle.load(open('HT2003_data_test100.pickle'))
    #for key in test_data: print key

    # Grab the first spectrum from the data
    x = test_data['x_values'][0]
    y = test_data['amplitudes'][0]

    pspline.prep_spectrum(x, y)

def test_fit_spline():

    import numpy as np
    import pspline
    import matplotlib.pyplot as plt
    from scipy.integrate import simps


    sigma = 1e-10
    x = np.linspace(-30, 30, 100)
    y = gauss(x, 5, -10, 10) + gauss(x, 40, 10, 4)
    y = gauss(x, 40, 10, 4) + np.random.normal(0, 0.1, len(y))
    y = gauss(x, 5, -10, 10) \
        + np.random.normal(0, sigma, len(y))
        #+ gauss(x, 110, 10, 4) \

    # Define range of chi values to fit with
    chis = [1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1, 1e2, 1e4]
    chis = np.logspace(-10, -1, 10)
    #chis = (1e-4,)

    A_C, h, lam_C, Vs = pspline.fit_spline(x, y, chis=chis)

    #N_k = lam_C.shape[0]
    #Delta = lam_C[1,0] - lam_C[0,0]
    #B = pspline.construct_B(N_k, N_k, Delta)

    plt.clf(); plt.close()
    if 0:
        plt.plot(chis, Vs)
        plt.xscale('log')
        plt.yscale('log')
        plt.show()

    #deriv_3 = B*h
    #deriv_2 = B*deriv_3
    #deriv_1 = B*deriv_2
    #deriv_0 = B*deriv_1

    plt.plot(lam_C, h)
    #plt.plot(lam_C, deriv_0)
    plt.plot(x, A_C)
    plt.plot(x, y)
    plt.savefig('splines.png')
    #plt.show()

def test_B():

    import pspline
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.integrate import simps

    Delta = 2

    #B_calc = pspline.construct_B(4, 4, Delta)

    #print B_calc.T

    B_true = np.array((
                       (0, 0.5*Delta, 0.5*Delta, 0.5*Delta),
                       (0, 0.5*Delta, Delta,     Delta,),
                       (0, 0,         0.5*Delta, Delta),
                       ))

    #assert np.array_equal(B_calc, B_true)

    # Test integration of fourth derivative of gaussians
    ngauss = 1
    x0s = (0, 10)
    sigmas = (5, 10)
    As = (10, 5)

    x = np.linspace(-30, 30, 100)
    x = np.arange(-30, 30, 0.5)
    if ngauss == 2:
        y = gauss(x, sigmas[0], x0s[0], As[0]) \
            + gauss(x, sigmas[1], x0s[1], As[1])
        y_4d = gauss_4th_deriv(x, sigmas[0], x0s[0], As[0]) \
               + gauss_4th_deriv(x, sigmas[1], x0s[1], As[1])
        y_1d = gauss_1st_deriv(x, sigmas[0], x0s[0], As[0]) \
               + gauss_1st_deriv(x, sigmas[1], x0s[1], As[1])
    elif ngauss == 1:
        y = gauss(x, sigmas[0], x0s[0], As[0])
        y_4d = gauss_4th_deriv(x, sigmas[0], x0s[0], As[0])
        y_1d = gauss_1st_deriv(x, sigmas[0], x0s[0], As[0])

    y_4d = np.matrix(y_4d).T
    y_1d = np.matrix(y_1d).T

    A_M, lam_M, lam_C, N_k = pspline.prep_spectrum(x, y)

    lam_0 = lam_M[0, 0]
    Delta = x[1] - x[0]
    N_D = len(x)
    print 'lam_0', lam_0
    print 'Delta', Delta

    B = pspline.construct_B(N_D, N_D, Delta, lam_C, lam_0).T

    print B.shape, y_4d.shape
    print B

    # Integrate to original function
    y_calc = B * y_4d

    y = np.squeeze(np.asarray(y))
    y_1d = np.squeeze(np.asarray(y_1d))
    y_calc = np.squeeze(np.asarray(y_calc))

    # Plot
    scale = np.max(y) / np.max(y_4d)
    plt.clf(); plt.close()
    plt.plot(x, y, label='f(x)',)
    plt.plot(x, y_4d * scale, label="f''''(x) x " + str(scale))
    plt.plot(x, y_calc, label="Integrated f''''(x)")
    plt.plot(x, y_calc - y, label="Integrated f''''(x) - f(x)")
    plt.xlim(-30, 30)
    plt.ylim(-15, 15)
    plt.legend(loc='best')
    plt.savefig('integration_test.png')

def test_beta():

    import pspline
    import numpy as np

    lam_C = np.linspace(0, 1, 5)

    beta_calc = pspline.construct_beta(lam_C)

    beta_true = np.array((
                       (1, -2, 1,  0,  0,  0, 0, 0, 0),
                       (0, 1,  -2, 1,  0,  0, 0, 0, 0),
                       (0, 0,  1,  -2, 1,  0, 0, 0, 0),
                       ))

    assert np.array_equal(beta_calc, beta_true)

def main():

    #test_prep_spectrum()
    test_B()
    #test_beta()
    #test_fit_spline()

    import pickle
    import pspline

    test_data = pickle.load(open('HT2003_data_test100.pickle'))
    #for key in test_data: print key

    # Grab the first spectrum from the data
    x = test_data['x_values'][0]
    y = test_data['data_list'][0]

    import csv
    with open('spectrum0.csv', 'wb') as csvfile:
        csv_file = csv.writer(csvfile, delimiter=' ',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for i in xrange(0, len(x)):
            csv_file.writerow((x[i], y[i]))

if __name__ == '__main__':
    main()



