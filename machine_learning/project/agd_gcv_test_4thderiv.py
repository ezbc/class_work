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

    expo = np.exp(-(x - x0)**2 / (2.0 * sigma**2))

    result = A * expo * (x - x0)**4.0 / sigma**8.0 - \
             6.0 * A * expo * (x - x0)**2.0 / sigma**6.0 + \
             3.0 * A * expo / sigma**4.0

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

    sigma = 1
    x = np.linspace(-30, 30, 71)
    y = gauss(x, 5, -10, 10) + gauss(x, 40, 10, 4)
    y = gauss(x, 40, 10, 4) + np.random.normal(0, 0.1, len(y))
    y = gauss(x, 5, -10, 10) \
        + np.random.normal(0, sigma, len(y))
        #+ gauss(x, 110, 10, 4) \
    y_1d = gauss_1st_deriv(x, 5, -10, 10)

    # Define range of chi values to fit with
    chis = [1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1, 1e2, 1e4]
    chis = np.logspace(-14, 5, 10)
    #chis = np.logspace(-4, 0, 10)
    #chis = np.logspace(-2, -1, 10)
    #chis = (1e-4,)

    #A_C, h, lam_C, Vs, derivs = pspline.fit_spline(x, y, chis=chis)
    A_C, h, derivs, lam_C = pspline.fit_spline(x, y)

    #N_k = lam_C.shape[0]
    #Delta = lam_C[1,0] - lam_C[0,0]
    #B = pspline.construct_B(N_k, N_k, Delta)

    if 0:
        plt.clf(); plt.close()
        plt.plot(chis, Vs)
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel(r'$\chi$')
        plt.xlabel(r'V($\chi$)')
        plt.savefig('figures/v_vs_chi.png')

    #deriv_3 = B*h
    #deriv_2 = B*deriv_3
    #deriv_1 = B*deriv_2
    #deriv_0 = B*deriv_1
    A_C = np.squeeze(np.asarray(A_C))

    plt.clf(); plt.close()
    plt.plot(lam_C, derivs[2], label=r"f'($\lambda$)")
    plt.plot(x, y_1d, label=r"True f'($\lambda$)")
    plt.legend(loc='best')
    plt.savefig('figures/spline_1stderiv.png')

    plt.clf(); plt.close()
    scale = np.max(y) / np.max(h)
    plt.plot(lam_C, h*scale, label=r"f''''($\lambda$)")
    plt.plot(lam_C, derivs[0], label=r"f'''($\lambda$)")
    plt.plot(lam_C, derivs[1], label=r"f''($\lambda$)")
    plt.plot(lam_C, derivs[2], label=r"f'($\lambda$)")
    plt.plot(lam_C, derivs[3], label=r"f($\lambda$)")
    #plt.plot(lam_C, deriv_0)
    plt.plot(x, A_C, label=r"Integ f''''($\lambda$)")
    plt.plot(x, y, label=r"f($\lambda$)")
    plt.plot(x, y - A_C, label=r"f''''($\lambda$) - Integ f''''($\lambda$)")
    plt.xlabel(r'$\lambda$')
    #plt.legend(loc='best')
    plt.savefig('figures/splines.png')
    #plt.show()

def test_B():

    import pspline
    import numpy as np
    import matplotlib.pyplot as plt

    lam_M = np.matrix((0, 2, 4))
    lam_C = np.matrix((0, 1, 2, 3, 4))
    A_M = 1/120.0 * np.power(lam_M, 5)
    A_M_4d = np.power(lam_C, 1) #* 1/24.

    N_k = lam_C.shape[1]
    N_D = lam_M.shape[1]

    # Prep arrays to be of correct shape
    A_M, lam_M, lam_C, N_k = pspline.prep_spectrum(lam_M, A_M, N_k=N_k)

    lam_0 = lam_M[0, 0]
    Delta = lam_M[1, 0] - lam_M[0, 0]

    B = pspline.construct_B(lam_C, lam_M)

    B_true = Delta / 12.0 * np.array((
                                      (0, 0, 0, 0, 0),
                                      (0, 8, 0, 0, 0),
                                      (0, 16, 64, 0, 0),
                                      ))

    from pprint import pprint

    if 0:
        print 'B calculated'
        print B / np.min(B[B != 0])
        print 'B True'
        print B_true / Delta * 12.0

    #assert np.array_equal(B, B_true)

    # Integrate to original function
    A_C = B* A_M_4d.T

    A_M = np.squeeze(np.asarray(A_M))
    A_M_4d = np.squeeze(np.asarray(A_M_4d))
    A_C = np.squeeze(np.asarray(A_C))

    if 0:
        print 'A_C = '
        print A_C
        print ''
        print 'A_M = '
        print A_M
        print ''
    if 1:
        print 'A_C / A_M = '
        print A_C / A_M
        print ''

    x_C = np.squeeze(np.array(x_C))
    #y_4d = np.squeeze(np.array(y_4d_calc[:-4]))
    y_4d = np.squeeze(np.array(y_4d))
    y_3d = cumtrapz(y_4d, x_C, initial=0)
    y_2d = cumtrapz(y_3d, x_C, initial=0)
    y_1d = cumtrapz(y_2d, x_C, initial=0)
    y_0d = cumtrapz(y_1d, x_C, initial=0)

    # Plot
    plt.clf(); plt.close()
    plt.plot(lam_M, A_M, label='f(x)',)
    plt.plot(lam_C, A_M_4d, label="f''''(x)")
    plt.plot(lam_M, A_C, label="Integrated f''''(x)")
    #plt.plot(lam_M, A_C - A_M, label="Integrated f''''(x) - f(x)")
    #plt.xlim(-30, 30)
    #plt.ylim(-15, max(As)*1.1)
    plt.legend(loc='best')
    plt.savefig('figures/linear_integration_test.png')

def test_gauss_integration(ngauss=1):

    import pspline
    import numpy as np
    from scipy import linalg
    import matplotlib.pyplot as plt

    # Test integration of fourth derivative of gaussians
    Delta = 0.5
    x0s = (0, 15)
    sigmas = (5, 10)
    As = (20, 10)
    x0s = (-5, 0)
    sigmas = (5, 5)
    As = (10, 20)

    x = np.arange(-30, 30, Delta)
    x = np.linspace(-30, 30, 71)

    if ngauss == 2:
        y = gauss(x, sigmas[0], x0s[0], As[0]) \
            + gauss(x, sigmas[1], x0s[1], As[1]) \
            + np.random.normal(0,1, len(x))
    elif ngauss == 1:
        y = gauss(x, sigmas[0], x0s[0], As[0]) \
            + np.random.normal(0,0.1, len(x))

    A_M, lam_M, lam_C, N_k = pspline.prep_spectrum(x, y)

    x_C = np.squeeze(np.asarray(lam_C))

    if ngauss == 2:
        y_4d = gauss_4th_deriv(x_C, sigmas[0], x0s[0], As[0]) \
               + gauss_4th_deriv(x_C, sigmas[1], x0s[1], As[1])
    elif ngauss == 1:
        y_4d = gauss_4th_deriv(x_C, sigmas[0], x0s[0], As[0])

    y_4d = np.matrix(y_4d).T
    #y_1d = np.matrix(y_1d).T

    lam_0 = lam_M[0, 0]
    Delta = x[1] - x[0]
    N_D = len(x)
    #print 'lam_0', lam_0
    #print 'Delta', Delta

    B = pspline.construct_B(lam_C, lam_C)

    trap_integ = pspline.construct_trap_integ(lam_C, lam_C)

    #print B.shape, y_4d.shape
    #print B

    # Integrate to original function
    y_calc = B[:-2,:] * y_4d[:, 0]

    y = np.squeeze(np.asarray(y))
    #y_1d = np.squeeze(np.asarray(y_1d))
    y_calc = np.squeeze(np.asarray(y_calc))

    lam_0 = lam_M[0, 0]
    Delta = abs(lam_C[-1, 0] - lam_C[0, 0]) / (N_k - 1)
    ones = np.matrix(np.ones(A_M.shape[0])).T
    lam_M0 = lam_M - ones * lam_0
    B = pspline.construct_B(lam_C, lam_M)
    B_prime = np.hstack((B,
                         ones,
                         lam_M0,
                         np.power(lam_M0, 2) / 2.0,
                         np.power(lam_M0, 3) / 6.0))
    beta = pspline.construct_beta(lam_C)
    chi = 10
    y_4d_calc = linalg.inv(B_prime.T * B_prime + chi / Delta**4 * \
            beta.T*beta) * B_prime.T * A_M

    #print y_calc.shape
    #print y_4d.shape

    #x_C = x_[1:-1]

    # Integrate 4th derive to 0th deriv
    from scipy.integrate import cumtrapz

    x_C = np.squeeze(np.array(x_C))
    y_4d = np.squeeze(np.array(y_4d_calc[:-4]))
    #y_4d = np.squeeze(np.array(y_4d))
    if 1:
        y_3d = cumtrapz(y_4d, x_C, initial=y_4d_calc[-4, 0])
        y_2d = cumtrapz(y_3d, x_C, initial=y_4d_calc[-3, 0])
        y_1d = cumtrapz(y_2d, x_C, initial=y_4d_calc[-2, 0])
        y_0d = cumtrapz(y_1d, x_C, initial=y_4d_calc[-1, 0])
    else:
        y_3d = cumtrapz(y_4d, x_C, initial=0)
        y_2d = cumtrapz(y_3d, x_C, initial=0)
        y_1d = cumtrapz(y_2d, x_C, initial=0)
        y_0d = cumtrapz(y_1d, x_C, initial=0)

    # Calculate spectrum to compare with computed
    if ngauss == 2:
        y_C = gauss(x_C, sigmas[0], x0s[0], As[0]) \
            + gauss(x_C, sigmas[1], x0s[1], As[1])
    elif ngauss == 1:
        y_C = gauss(x_C, sigmas[0], x0s[0], As[0])

    norm = np.sum((y_0d - y_C)**2)
    print('L2 norm of calc y and true y = {0:.2f}'.format(norm))

    print('\ny_0d[0] =', y_0d[0])

    # Plot
    scale = np.max(y) / np.max(y_4d)
    plt.clf(); plt.close()
    plt.plot(x, y, label='f(x)',)
    plt.plot(x_C, y_4d * scale, label="f''''(x) x " + str(scale))
    plt.plot(x_C, y_3d, label="Calc. f'''(x)")
    plt.plot(x_C, y_2d, label="Calc. f''(x)")
    plt.plot(x_C, y_1d, label="Calc. f'(x)")
    plt.plot(x_C, y_0d, label="Calc. f(x)",)#marker='o')
    plt.plot(x_C[:-2], y_calc, label="Integrated f''''(x)",) #marker='+')
    #plt.plot(x, y_calc - y, label="Integrated f''''(x) - f(x)")
    plt.xlim(-30, 30)
    #plt.ylim(-15, max(As)*1.1)
    plt.legend(loc='best')
    plt.savefig('figures/integration_test_ngauss' + str(ngauss) + '.png')
    plt.show()

    plt.clf(); plt.close()
    plt.plot(x_C, y_4d, label="f''''(x)")
    plt.plot(x_C, y_4d_calc[:-4], label="f''''(x) Calc")
    #plt.plot(x_C, y_3d, label="Calc. f'''(x)")
    #plt.plot(x_C, y_2d, label="Calc. f''(x)")
    #plt.plot(x_C, y_1d, label="Calc. f'(x)")
    #plt.plot(x_C, y_0d, label="Calc. f(x)")
    #plt.plot(x, y_calc - y, label="Integrated f''''(x) - f(x)")
    plt.xlim(-30, 30)
    #plt.ylim(-15, max(As)*1.1)
    plt.legend(loc='best')
    plt.savefig('figures/integration_test_derivs_ngauss_' + str(ngauss) + '.png')

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

def plot_splines(results, show=False):

    # Import external modules
    import numpy as np
    from mpl_toolkits.axes_grid1 import ImageGrid
    import matplotlib.pyplot as plt
    import matplotlib

    # Write data
    x_data = results['x_data']
    x_fit = results['x_fit']

    for i in xrange(len(results['y_data_list'])):
        # Set up plot aesthetics
        plt.clf()
        plt.close()
        plt.rcdefaults()
        colormap = plt.cm.gist_ncar
        font_scale = 10
        params = {#'backend': .pdf',
                  'axes.labelsize': font_scale,
                  'axes.titlesize': font_scale,
                  'text.fontsize': font_scale,
                  'legend.fontsize': font_scale * 3 / 4.0,
                  'xtick.labelsize': font_scale,
                  'ytick.labelsize': font_scale,
                  'font.weight': 500,
                  'axes.labelweight': 500,
                  'text.usetex': False,
                  #'figure.figsize': (8, 8 * y_scaling),
                  #'axes.color_cycle': color_cycle # colors of different plots
                 }
        plt.rcParams.update(params)

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9,4))

        # Plot spline fit to data
        # -----------------------
        # Write data
        y_data = results['y_data_list'][i]
        y_fit = results['y_fit_list'][i]

        ax1.plot(x_fit, y_fit,
                 label=r"Spline Fit",
                 linestyle='-',
                 color='k')
        ax1.plot(x_data, y_data,
                 label=r"HT03 Data",
                 linestyle='-',
                 color='r',
                 #marker='o',
                 )
        ax1.set_xlabel(r'Velocity (km/s)')
        ax1.set_ylabel(r'T$_{\rm obs}$ (K km/s)')
        ax1.legend(loc='lower left')
        ax1.set_xlim(-40, 40)

        # Plot derivatives
        # -----------------------
        # Write data
        y_4d = results['y_4d_list'][i]
        y_3d = results['y_3d_list'][i]
        y_2d = results['y_2d_list'][i]
        y_1d = results['y_1d_list'][i]

        ax2.plot(x_fit, y_4d / y_4d.max(),
                 label=r"Spline Fit 4th Deriv.",
                 linestyle='-',
                 color='r')
        ax2.plot(x_fit, y_3d / y_3d.max(),
                 label=r"Spline Fit 3rd Deriv.",
                 linestyle='-',
                 color='c')
        ax2.plot(x_fit, y_2d / y_2d.max(),
                 label=r"Spline Fit 2nd Deriv.",
                 linestyle='-',
                 color='g')
        ax2.plot(x_fit, y_1d / y_1d.max(),
                 label=r"Spline Fit 1st Deriv.",
                 linestyle='-',
                 color='b')
        ax2.plot(x_fit, y_fit / y_fit.max(),
                 label=r"Spline Fit",
                 linestyle='-',
                 color='k')
        ax2.set_xlabel(r'Velocity (km/s)')
        ax2.set_ylabel(r'T$_{\rm obs}$ / max(T$_{\rm obs}$)')
        ax2.legend(loc='lower right')
        ax2.set_xlim(-100, 100)
        ax2.set_ylim(-2, 1.1)

        plt.savefig('figures/ht03_fits/ht03_spline{0:.0f}'.format(i) + '.png',
                    bbox_inches='tight')
        if 0:
            plt.clf(); plt.close()
            plt.plot(x_fit, y_4d / y_4d.max(),
                     label=r"Spline Fit 4th Deriv.",
                     linestyle='-',
                     color='r')
            plt.show()

        if show:
            plt.show()

def main():

    #test_prep_spectrum()
    #test_B()
    #test_gauss_integration(ngauss=1)
    #test_gauss_integration(ngauss=2)
    #test_beta()
    #test_fit_spline()

    import pickle
    import pspline
    import numpy as np

    test_data = pickle.load(open('data/HT2003_data_test100.pickle'))

    # load data instead of fitting?
    load_data = 0

    # Grab the first spectrum from the data
    x = test_data['x_values'][0]
    y_list = test_data['data_list'][0:2]

    if not load_data:
        A_C_list, h_list, derivs_list, lam_C = \
                pspline.fit_spline(x, y_list, N_k=len(x), init_guess=0.0052)

        y_3d_list = []
        y_2d_list = []
        y_1d_list = []
        for deriv_list in derivs_list:
        	y_3d_list.append(deriv_list[0])
        	y_2d_list.append(deriv_list[1])
        	y_1d_list.append(deriv_list[2])

        if len(y_list) > 1:
            results = {}
            results['y_fit_list'] = A_C_list
            results['y_4d_list'] = h_list
            results['y_3d_list'] = y_3d_list
            results['y_2d_list'] = y_2d_list
            results['y_1d_list'] = y_1d_list
            results['x_fit'] = lam_C
            results['y_data_list'] = y_list
            results['x_data'] = x
        elif len(y_list) == 1:
            results = {}
            results['y_fit_list'] = (A_C_list,)
            results['y_4d_list'] = (h_list,)
            results['y_3d_list'] = y_3d_list
            results['y_2d_list'] = y_2d_list
            results['y_1d_list'] = y_1d_list
            results['x_fit'] = lam_C
            results['y_data_list'] = (np.array(y_list[0]),)
            results['x_data'] = x

        with open('data/spline_fits.pickle', 'w') as f:
            pickle.dump(results, f)
    elif load_data:
        results = pickle.load(open('data/spline_fits.pickle'))

    plot_splines(results, show=0)

    #import csv
    #with open('spectrum0.csv', 'wb') as csvfile:
    #    csv_file = csv.writer(csvfile, delimiter=' ',
    #                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
    #    for i in xrange(0, len(x)):
    #        csv_file.writerow((x[i], y[i]))

if __name__ == '__main__':
    main()



