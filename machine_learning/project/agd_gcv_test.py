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
    x0s = (3, 14, 16)
    sigmas = (1.7, 2, 0.4)
    As = (0.6, 1.055, 0.1)

    x = np.arange(-30, 30, Delta)
    x = np.linspace(-30, 30, 71)
    x = np.arange(-98, 98, 0.20010214)

    if ngauss == 2:
        y = gauss(x, sigmas[0], x0s[0], As[0]) \
            + gauss(x, sigmas[1], x0s[1], As[1]) \
            + gauss(x, sigmas[2], x0s[2], As[2])
    elif ngauss == 1:
        y = gauss(x, sigmas[0], x0s[0], As[0]) \
            + np.random.normal(0,1, len(x))

    # Add noise
    uniform = 0
    if uniform:
        y += np.random.normal(0,0.008007, len(x))
    else:
        noise = np.random.normal(0,0.008007, len(x))
        indices = np.where((x > -1) & (x < 17))
        noise[indices] = np.random.normal(0,0.008007/3.0, len(indices[0]))
        y += noise

    y_list = [y,]

    if 0:
        import matplotlib.pyplot as plt
        plt.plot(x, y_list[0])
        plt.show()

    print('sim ylist length', len(y_list))

    temp_results = pspline.fit_spline(x, y_list, N_k=len(x), init_guess=20)

    y_fit = temp_results['reg 1st deriv']['y_fit_list'][0]
    y_4d = temp_results['reg 4th deriv']['y4d_list'][0]
    y_3d = temp_results['reg 3rd deriv']['y3d_list'][0]
    y_2d = temp_results['reg 2nd deriv']['y2d_list'][0]
    y_1d = temp_results['reg 1st deriv']['y1d_list'][0]
    x_fit = temp_results['x_fit']
    y_calc = temp_results['reg 1st deriv']['y_fit_list'][0]
    x_data = temp_results['x_data']

    #norm = np.sum(y_calc - y)
    #print('Sum of diff. of calc y and true y = {0:.2f}'.format(norm))
    norm = np.sum((y_calc - y)**2)
    print('L2 norm of calc y and true y = {0:.2f}'.format(norm))

    # Plot
    scale = np.max(y) / np.max(y_4d)
    plt.clf(); plt.close()
    plt.plot(x, y, label='f(x)',)
    plt.plot(x_fit, y_4d * scale, label="f''''(x) x " + str(scale))
    plt.plot(x_fit, y_3d, label="Calc. f'''(x)")
    plt.plot(x_fit, y_2d, label="Calc. f''(x)")
    plt.plot(x_fit, y_1d, label="Calc. f'(x)")
    plt.plot(x_fit, y_calc, label="Calc. f(x)",)#marker='o')
    #plt.plot(x, y_calc - y, label="Integrated f''''(x) - f(x)")
    plt.xlim(-30, 30)
    #plt.ylim(-15, max(As)*1.1)
    plt.legend(loc='best')
    plt.savefig('figures/integration_test_ngauss' + str(ngauss) + '.png')
    plt.show()

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

def plot_splines(results, show=False, filename=None, wide=False):

    # Import external modules
    import numpy as np
    from mpl_toolkits.axes_grid1 import ImageGrid
    import matplotlib.pyplot as plt
    import matplotlib

    # Write data
    x_data = results['x_data']
    x_fit = results['x_fit']
    if wide:
        filename += '_wide'

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

        if not wide:
            fig = plt.figure(figsize=(5,8))
        else:
            fig = plt.figure(figsize=(10, 5))

        # Plot spline fit to data
        # -----------------------
        # Write data
        y_data = results['y_data_list'][i]
        y_fit = results['y_fit_list'][i]

        if not wide:
            ax1 = fig.add_subplot(211)
        else:
            ax1 = fig.add_subplot(121)

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
        ax1.legend(loc='upper left')
        ax1.set_xlim(-40, 40)


        # Plot derivatives
        # -----------------------
        nrows_ncols=(2,2)
        ngrids=4
        if not wide:
            subplot = (2,1,2)
        else:
            subplot = (1,2,2)

        imagegrid = ImageGrid(fig, subplot,
                     nrows_ncols=nrows_ncols,
                     ngrids=ngrids,
                     axes_pad=0.1,
                     aspect=False,
                     label_mode='L',
                     share_all=True)

        # Write data
        y_4d = results['y_4d_list'][i]
        y_3d = results['y_3d_list'][i]
        y_2d = results['y_2d_list'][i]
        y_1d = results['y_1d_list'][i]

        ax2 = imagegrid[0]
        ax2.plot(x_fit, y_4d / y_4d.max(),
                 label=r"Spline Fit 4th Deriv.",
                 linestyle='-',
                 color='r')
        ax2.legend(loc='lower right')

        ax3 = imagegrid[1]
        ax3.plot(x_fit, y_3d / y_3d.max(),
                 label=r"Spline Fit 3rd Deriv.",
                 linestyle='-',
                 color='c')
        ax3.legend(loc='lower right')


        ax4 = imagegrid[2]
        ax4.plot(x_fit, y_2d / y_2d.max(),
                 label=r"Spline Fit 2nd Deriv.",
                 linestyle='-',
                 color='g')
        ax4.legend(loc='lower right')

        ax5 = imagegrid[3]
        ax5.plot(x_fit, y_1d / y_1d.max(),
                 label=r"Spline Fit 1st Deriv.",
                 linestyle='-',
                 color='b')
        ax5.legend(loc='lower right')

        ax2.set_xlabel(r'Velocity (km/s)')
        ax2.set_ylabel(r'T$_{\rm obs}$ / max(T$_{\rm obs}$)')
        ax3.set_xlabel(r'Velocity (km/s)')
        ax3.set_ylabel(r'T$_{\rm obs}$ / max(T$_{\rm obs}$)')
        ax4.set_xlabel(r'Velocity (km/s)')
        ax4.set_ylabel(r'T$_{\rm obs}$ / max(T$_{\rm obs}$)')
        ax5.set_xlabel(r'Velocity (km/s)')
        ax5.set_ylabel(r'T$_{\rm obs}$ / max(T$_{\rm obs}$)')
        ax5.set_xlim(-39, 40)
        ax5.set_ylim(-2, 1.1)

        if filename is not None:
            if len(results['y_data_list']) > 1:

                plt.savefig(filename + '{:03d}'.format(i) + '.png',
                            bbox_inches='tight',
                            dpi=400)
            else:
                print('Writing figure ' + filename + '.png')
                plt.savefig(filename + '.png',
                            bbox_inches='tight',
                            dpi=400)

        if 0:
            plt.clf(); plt.close()
            plt.plot(x_fit, y_4d / y_4d.max(),
                     label=r"Spline Fit 4th Deriv.",
                     linestyle='-',
                     color='r')
            plt.show()

        if show:
            plt.show()

def run_ht03_sim(uniform_noise=True, load_data=False):

    import pspline
    import pickle
    import numpy as np
    from scipy import linalg
    import matplotlib.pyplot as plt

    if not load_data:
        # Test integration of fourth derivative of gaussians
        x0s = (3, 14, 16)
        sigmas = (1.7, 2, 0.4)
        As = (0.6, 1.055, 0.1)

        x = np.arange(-98, 98, 0.20010214)

        print len(x)

        y = gauss(x, sigmas[0], x0s[0], As[0]) \
            + gauss(x, sigmas[1], x0s[1], As[1]) \
            + gauss(x, sigmas[2], x0s[2], As[2])

        # Add noise
        if uniform_noise:
            y += np.random.normal(0,0.008007, len(x))
        else:
            noise = np.random.normal(0,0.008007, len(x))
            indices = np.where((x > -1) & (x < 5))
            noise[indices] = np.random.normal(0.01,0.008007/3.0,
                                              len(indices[0]))
            indices = np.where((x > 5) & (x < 10))
            noise[indices] = np.random.normal(-0.01,0.008007/3.0,
                                              len(indices[0]))
            indices = np.where((x > 10) & (x < 17))
            noise[indices] = np.random.normal(0,0.008007/1.0, len(indices[0]))
            indices = np.where((x > -20) & (x < -4))
            noise[indices] = np.random.normal(0.01,0.008007, len(indices[0]))
            indices = np.where((x > -30) & (x < -20))
            noise[indices] = np.random.normal(-0.01,0.008007/1.3, len(indices[0]))
            indices = np.where((x > 18) & (x < 25))
            noise[indices] = np.random.normal(0.01,0.008007*1.3, len(indices[0]))
            indices = np.where((x > 25) & (x < 40))
            noise[indices] = np.random.normal(-0.01,0.008007, len(indices[0]))
            y += noise

        y_list = [y,]

        temp_results = pspline.fit_spline(x, y_list, N_k=len(x),
                init_guess=1e-1)

        results = {}
        results['y_fit_list'] = temp_results['reg 4th deriv']['y_fit_list']
        results['y_4d_list'] = temp_results['reg 4th deriv']['y4d_list']
        results['y_3d_list'] = temp_results['reg 3rd deriv']['y3d_list']
        results['y_2d_list'] = temp_results['reg 2nd deriv']['y2d_list']
        results['y_1d_list'] = temp_results['reg 1st deriv']['y1d_list']
        results['x_fit'] = temp_results['x_fit']
        results['y_data_list'] = y_list
        results['x_data'] = temp_results['x_data']

        if uniform_noise:
            data_filename = 'ht03_sim_uniform_noise.pickle'
        else:
            data_filename = 'ht03_sim_nonuniform_noise.pickle'

        with open('data/' + data_filename, 'w') as f:
            pickle.dump(results, f)
    elif load_data:
        if uniform_noise:
            data_filename = 'ht03_sim_uniform_noise.pickle'
        else:
            data_filename = 'ht03_sim_nonuniform_noise.pickle'

        results = pickle.load(open('data/' + data_filename))

    #norm = np.sum(y_calc - y)
    #print('Sum of diff. of calc y and true y = {0:.2f}'.format(norm))
    norm = np.sum((results['y_data_list'][0] - results['y_fit_list'][0])**2)
    print('L2 norm of calc y and true y = {0:.2f}'.format(norm))

    # filename
    if uniform_noise:
        filename = 'ht03_sim_uniform_noise'
    else:
        filename = 'ht03_sim_nonuniform_noise'

    plot_splines(results, filename='figures/' + filename)

def run_ht03_fits(chi=None, load_data=False):

    import pickle
    import pspline
    import numpy as np
    import os

    test_data = pickle.load(open('data/HT2003_data_test100.pickle'))

    # Grab the first spectrum from the data
    x = test_data['x_values'][0]
    y_list = test_data['data_list'][0:]

    smooth = 0
    if smooth:
        x = np.linspace(x[0], x[-1], len(x) / 2)
        for i, y in enumerate(y_list):
            y_new = np.zeros(len(x))
            for j in xrange(len(x)):
                y_new[j] = np.average(y[2*j:2*j+1])
            y_list[i] = y_new


    # Do a fourier transform
    if 1:
        import matplotlib.pyplot as plt

        ft = np.fft.fft(y_list[0])
        freq = np.fft.fftfreq(y_list[0].shape[0])
        plt.plot(freq, ft)
        plt.show()
        return None


    if 0:
        print('last and first =', x[-1], x[0])
        print('noise =', np.std(y_list[0][0:50]))
        print('max = ', np.max(y_list[0]))
        print('snr = ', np.max(y_list[0])/ np.std(y_list[0][0:50]))
        print('Delta =', x[1] - x[0])
        print('ht03 ylist length', len(y_list))

    if 0:
        import matplotlib.pyplot as plt
        plt.plot(x, y_list[0])
        plt.show()

    if not load_data:

        temp_results = pspline.fit_spline(x, y_list, N_k=len(x),
                #init_guess=0.0001525,
                init_guess=1,
                chi=chi,
                )

        results = {}
        results['y_fit_list'] = temp_results['reg 4th deriv']['y_fit_list']
        results['y_4d_list'] = temp_results['reg 4th deriv']['y4d_list']
        results['y_3d_list'] = temp_results['reg 3rd deriv']['y3d_list']
        results['y_2d_list'] = temp_results['reg 2nd deriv']['y2d_list']
        results['y_1d_list'] = temp_results['reg 1st deriv']['y1d_list']
        results['x_fit'] = temp_results['x_fit']
        results['y_data_list'] = y_list
        results['x_data'] = temp_results['x_data']

        if chi is None:
            data_filename = 'spline_fits.pickle'
        else:
            data_filename = 'spline_fits_chi{0:.2f}.pickle'.format(chi)

        with open('data/' + data_filename, 'w') as f:
            pickle.dump(results, f)

    elif load_data:
        results = pickle.load(open('data/spline_fits.pickle'))

    # Plot
    if chi is None:
        filename = 'figures/ht03_fits/ht03_spline_'
    else:
        filename = 'figures/ht03_fits/ht03_spline_chi{0:.2f}_'.format(chi)

    plot_splines(results, show=0, filename=filename)
    plot_splines(results, show=0, filename=filename, wide=True)

    # Make a movie with the plots
    if chi is None:
        os.system('avconv -r 2 -i figures/ht03_fits/ht03_spline_%3d.png ' + \
                  'figures/ht03_fits.mp4')
    else:
        os.system('avconv -r 2 -i figures/ht03_fits/' + \
                  'ht03_spline_chi{0:.2f}_%3d.png '.format(chi) + \
                  'figures/ht03_fits_chi{0:.2f}.mp4'.format(chi))

def main():

    #test_prep_spectrum()
    #test_B()
    #test_gauss_integration(ngauss=1)
    #test_gauss_integration(ngauss=2)
    #test_beta()
    #test_fit_spline()

    # Simulations
    if 0:
        load_data = 0
        run_ht03_sim(uniform_noise=False, load_data=load_data)
        run_ht03_sim(uniform_noise=True, load_data=load_data)

    # Data fits
    if 1:
        load_data = 1
        run_ht03_fits(load_data=load_data)
        run_ht03_fits(chi=0.01, load_data=load_data)

if __name__ == '__main__':
    main()



