#!/usr/bin/python

def prep_spectrum(x, y, N_k=None):

    '''

    Parameters
    ----------
    x, y : array-like

    N_k : int
        Number of samples to create model derivatives. If not set, N_k = 100 *
        length(x).


    '''

    import numpy as np

    # Convert spectral array data into vectors
    A_M = np.matrix(y).T
    lam_M = np.matrix(x).T

    # Initialize model vectors, amps + wavelengths
    # --------------------------------------------
    # Define sampling size of model vectors
    if N_k is None:
        N_k = 10 * len(x)

    # wavelengths
    lam_C_array = np.linspace(lam_M[0, 0], lam_M[-1, 0], N_k)
    lam_C = np.matrix(lam_C_array).T

    return A_M, lam_M, lam_C, N_k

def fit_spline(x, y, N_k=None, chis=None):

    import numpy as np

    V_list = []
    h_list = []
    A_C_list = []
    coeffs_list = []
    for chi in chis:
        A_C, h, coeffs, lam_C, V = calc_V(x, y, chi, N_k=N_k)

        V_list.append(V)
        h_list.append(h)
        A_C_list.append(A_C)
        coeffs_list.append(coeffs)

    import matplotlib.pyplot as plt
    plt.clf(); plt.close()
    fig = plt.figure(figsize=(12,8))
    n = np.ceil(np.sqrt(len(h_list)))
    lam_C = np.squeeze(np.asarray(lam_C))
    for i in xrange(0, len(h_list)):
        h_hat = np.squeeze(np.asarray(h_list[i]))
        ax = fig.add_subplot(5, 5, i+1)
        ax.plot(lam_C, h_hat)
    plt.tight_layout()
    fig.savefig('figures/4th_derivs.png')

    # Get best-fit 4th deriv and associated spline
    h_hat = h_list[np.where(V_list == min(V_list))[0]]
    A_C_hat = A_C_list[np.where(V_list == min(V_list))[0]]
    coeffs_hat = coeffs_list[np.where(V_list == min(V_list))[0]]

    # Integrate up to original function
    from scipy.integrate import cumtrapz
    x_C = np.squeeze(np.array(lam_C))
    h_hat = np.squeeze(np.array(h_hat))
    y_3d = cumtrapz(h_hat, x_C, initial=coeffs_hat[0])
    y_2d = cumtrapz(y_3d, x_C, initial=coeffs_hat[-2])
    y_1d = cumtrapz(y_2d, x_C, initial=coeffs_hat[-3])
    y_0d = cumtrapz(y_1d, x_C, initial=coeffs_hat[-4])
    derivatives = (y_3d, y_2d, y_1d, y_0d)

    return A_C_hat, h_hat, lam_C, V_list, derivatives

def fit_spline(x, y, N_k=None, init_guess=None, verbose=True):

    '''
    Fits a spline to a one dimensional dataset using Tikhonov regularization of
    the 4th derivative.

    Parameters
    ----------
    x : array-like
        Spectral axis of spectra.
    y : array-like, list
        Amplitude as a function of x. If multiple arrays provided in a list,
        splines will be fitted to each of the y arrays. This reduces initial
        setup computation.

    '''


    import numpy as np
    from scipy.optimize import minimize
    from scipy.integrate import cumtrapz

    # Construct matrices calculating residual sum of squares, V(chi)
    # Subscript C = computed spectrum, sample length N_D
    # Subscript M = measured spectrum, sample length N_k
    # --------------------------------------------------------------

    if verbose:
        print('\nPrepping matrices...')

    N_D = len(x)

    # Check if y consists of just one spectrum or multiple spectra
    # ------------------------------------------------------------
    y = np.array(y)
    if y.ndim == 1:
    	A_M_list = (np.matrix(y).T,)
    elif y.ndim == 2:
        A_M_list = []
        y_list = y.tolist()
        for y in y_list:
        	A_M_list.append(np.matrix(y).T)
    else:
    	raise ValueError('y must be a one dimensional array or a list of' + \
    	                 'dimensional arrays of len(x)')

    # Initialize model vectors, wavelengths
    # --------------------------------------------
    # Define sampling size of model vectors
    if N_k is None:
        N_k = 10 * len(x)

    # Convert spectral array data into vectors
    lam_M = np.matrix(x).T

    # wavelengths
    lam_C_array = np.linspace(lam_M[0, 0], lam_M[-1, 0], N_k)
    lam_C = np.matrix(lam_C_array).T

    #assert A_M.shape == (N_D, 1)
    assert lam_M.shape == (N_D, 1)
    assert lam_C.shape == (N_k, 1)

    # Reference wavelength
    lam_0 = lam_M[0, 0]
    Delta = abs(lam_C[-1, 0] - lam_C[0, 0]) / (N_k - 1)

    # Ones is a column vector of ones of shape N_D x 1
    ones = np.matrix(np.ones(N_D)).T
    lam_M0 = lam_M - ones * lam_0
    assert ones.shape == lam_M0.shape

    B = construct_B(lam_C, lam_M)
    assert B.shape == (N_D, N_k)

    # See Eq 2b of Yeow et al. 2005

    # Define regularization parameters
    # --------------------------------
    # Condition (i) ensures that the computed
    # spectrum approximates the measured spectrum closely
    # Condition (ii) ensures that the fourth derivative spectrum does not show
    # spurious fluctuations.

    # S_1 = np.sum((A_M - A_C)**2)
    # S_2 = np.sum(np.diff(h)**2)

    # Tikhonov regularization parameter
    # R = S_1 + chi * S_2

    # add columns to B
    B_prime = np.hstack((B,
                         ones,
                         lam_M0,
                         np.power(lam_M0, 2) / 2.0,
                         np.power(lam_M0, 3) / 6.0))
    assert B_prime.shape == (N_D, N_k + 4)

    # Construct beta
    beta = construct_beta(lam_C)
    assert beta.shape[1] == N_k + 4

    # Fit each spectrum provided
    # --------------------------
    A_C_hat_list = []
    h_hat_list = []
    derivs_list = []
    chi_hat = None

    for i, A_M in enumerate(A_M_list):

        # Begin minimization of V(chi)
        # ----------------------------
        if verbose:
            if len(A_M_list) > 1:
                print('\nPerforming spline fit of spectrum {0:.0f}'.format(i))
            else:
                print('\nPerforming spline fit')

        if init_guess is None:
            init_guess = 1

        # Change initial guess to best guess of previous spectrum
        if chi_hat is not None:
        	init_guess = chi_hat

        args = (N_D, A_M, lam_M, lam_C, N_k, B_prime, beta, lam_0, Delta)

        result = minimize(calc_V, x0=init_guess, args=args,
                          method='Nelder-Mead', tol=1e-6)
        chi_hat = result.x[0]
        #chi_hat = 0.005

        verbose = True
        if verbose:
            print('Chi minimum = ' + str(chi_hat))

        # Get best-fit model 4th deriv and spectrum
        A_C_hat, h_hat, coeffs = calc_4th_deriv(chi_hat,
                                                B_prime,
                                                Delta,
                                                beta,
                                                A_M,
                                                N_k)

        # Integrate up to original function
        x_C = np.squeeze(np.array(lam_C))
        h_hat = np.squeeze(np.array(h_hat))
        if 1:
            y_3d = cumtrapz(h_hat, x_C, initial=coeffs[3])
            y_2d = cumtrapz(y_3d, x_C, initial=coeffs[2])
            y_1d = cumtrapz(y_2d, x_C, initial=coeffs[1])
            y_0d = cumtrapz(y_1d, x_C, initial=coeffs[0])
        if 0:
            y_3d = cumtrapz(h_hat, x_C, initial=0) + coeffs[3]
            y_2d = cumtrapz(y_3d, x_C, initial=0) + coeffs[2]
            y_1d = cumtrapz(y_2d, x_C, initial=0) + coeffs[1]
            y_0d = cumtrapz(y_1d, x_C, initial=0) + coeffs[0]
        derivs = (y_3d, y_2d, y_1d, y_0d)

        # Convert matrices to arrays, since matrices are annoying
        A_C_hat = np.squeeze(np.array(A_C_hat))
        h_hat = np.squeeze(np.array(h_hat))
        lam_C = np.squeeze(np.array(lam_C))

        A_C_hat_list.append(A_C_hat)
        h_hat_list.append(h_hat)
        derivs_list.append(derivs)

    if len(A_M_list) == 1:
        return A_C_hat_list[0], h_hat_list[0], derivs_list[0], lam_C
    else:
        return A_C_hat_list, h_hat_list, derivs_list, lam_C

def calc_V(chi, N_D=None, A_M=None, lam_M=None, lam_C=None, N_k=None,
        B_prime=None, beta=None, lam_0=None, Delta=None):

    '''

    A_C = computed spectrum [N_D x 1]
    B = coefficients [N_D x N_k]
    ones = column vector of ones [N_D x 1]
    lambda_M = measured wavelengths [N_D x 1]
    lambda_C = computed wavelengths [N_k x 1]

    '''

    import numpy as np
    import numpy
    import scipy
    from scipy import linalg

    assert B_prime.shape == (N_D, N_k + 4)
    assert beta.shape[1] == N_k + 4

    # Chi is a list from optimize
    if type(chi) is list or type(chi) is numpy.ndarray:
        chi = chi[0]

    # Chi cannot be zero, else finite differencing and no need for this at all!
    if chi == 0:
        return np.Inf

    # h_E_matrix is computation heavy, calculate for h_prime and E both
    h_E_matrix = calc_h_E_matrix(chi, B_prime, Delta, beta)

    # Calculate the fourth derivative and coefficients
    A_C, h, coeffs = calc_4th_deriv(chi, B_prime, Delta, beta, A_M, N_k,
                                    h_E_matrix=h_E_matrix)

    # Eq 8
    E = B_prime * h_E_matrix

    # The optimum chi comes from minimizing V(chi)
    V = ((A_M - A_C).T * (A_M - A_C) / N_D) \
        / (1 - np.trace(E) / N_D)**2

    return V[0,0]

def calc_h_E_matrix(chi, B_prime, Delta, beta):

    from scipy import linalg

    # Write matrix shared by h_prime and E
    h_E_matrix = linalg.inv(B_prime.T * B_prime + chi / Delta**4 * \
                            beta.T * beta) * B_prime.T

    return h_E_matrix

def calc_4th_deriv(chi, B_prime, Delta, beta, A_M, N_k, h_E_matrix=None):

    ''' Calculates 4th deriv and estimated function

    If h_E_matrix is provide, used to calculate h_prime

    '''

    import numpy as np

    # B_prime^T * B_prime --> N_k+4 x N_D * N_D x N_k+4 --> N_k+4 x N_k+4
    # beta^T * beta --> N_k+4 x N_k-2 * N_k-2 x N_k+4 --> N_k+4 x N_k+4

    # Solve for h_prime, eq 4
    if h_E_matrix is None:
        h_prime = calc_h_E_matrix(chi, B_prime, Delta, beta) * A_M
    else:
        h_prime = h_E_matrix * A_M
    assert h_prime.shape == (N_k+4, 1)

    coeffs = np.squeeze(np.array(h_prime[-4:])).tolist()
    h_hat = h_prime[0:N_k]

    # Eq 6, N_D x N_k+4 * N_k+4 x 1= 1 x N_D
    A_C = B_prime * h_prime
    assert A_C.shape == A_M.shape

    return A_C, h_hat, coeffs

def construct_trap_integ(lam_C, lam_M):

    import numpy as np

    # Prep Original variables
    #lam_M = np.reshape(lam_M, (N_D, 1))
    #lam_C = np.reshape(lam_C, (N_k, 1))
    lam_M = np.squeeze(np.asarray(lam_M))
    lam_C = np.squeeze(np.asarray(lam_C))
    lam_0 = lam_M[0]
    N_D = lam_M.size
    N_k = lam_C.size

    Delta = (lam_C[-1] - lam_0) / (N_k - 1.0)

    # Prep scalings of axes
    dxk = 1.0 / (N_k - 1.0)
    dxu = (lam_M - lam_0) / (lam_M[-1] - lam_0)
    xu = lam_M
    xknot = np.linspace(0,1,N_k)
    dxl = np.zeros((N_D))
    xl = dxl
    xu = dxu
    trap_integ = np.matrix(np.zeros((N_D, N_k)))

    testing = False

    if testing:
        print 'dxk', dxk
        print 'dxu', dxu
        print 'xu', xu
        print 'xknot', xknot
        print 'dxl', dxl
        print 'xl', xl

    if testing:
        assert dxk == 0.0025
        np.testing.assert_array_almost_equal(dxl,
            np.array(
            (0.,0.,0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
            0.,0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
            0.,0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
            0.,0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
            0.,0., 0.)),
            decimal=4)
        np.testing.assert_array_almost_equal(dxu,
            np.array(
            (0., 0.0142857, 0.0285714, 0.0428571, 0.0571429, 0.0714286,
            0.0857143, 0.1, 0.114286, 0.128571, 0.142857, 0.157143, 0.171429,
            0.185714, 0.2, 0.214286, 0.228571, 0.242857, 0.257143, 0.271429,
            0.285714, 0.3, 0.314286, 0.328571, 0.342857, 0.357143, 0.371429,
            0.385714, 0.4, 0.414286, 0.428571, 0.442857, 0.457143, 0.471429,
            0.485714, 0.5, 0.514286, 0.528571, 0.542857, 0.557143, 0.571429,
            0.585714, 0.6, 0.614286, 0.628571, 0.642857, 0.657143, 0.671429,
            0.685714, 0.7, 0.714286, 0.728571, 0.742857, 0.757143, 0.771429,
            0.785714, 0.8, 0.814286, 0.828571, 0.842857, 0.857143, 0.871429,
            0.885714, 0.9, 0.914286, 0.928571, 0.942857, 0.957143, 0.971429,
            0.985714, 1.)),
            decimal=4)
        np.testing.assert_array_almost_equal(xknot,
            np.array(
            (0.,0.0025,0.005,0.0075,0.01, 0.0125, 0.015, 0.0175, 0.02, 0.0225, \
            0.025, 0.0275, 0.03, 0.0325, 0.035, 0.0375, 0.04, 0.0425, 0.045, \
            0.0475, 0.05, 0.0525, 0.055, 0.0575, 0.06, 0.0625, 0.065, 0.0675, \
            0.07, 0.0725, 0.075, 0.0775, 0.08, 0.0825, 0.085, 0.0875, 0.09, \
            0.0925, 0.095, 0.0975, 0.1, 0.1025, 0.105, 0.1075, 0.11, 0.1125, \
            0.115, 0.1175, 0.12, 0.1225, 0.125, 0.1275, 0.13, 0.1325, 0.135, \
            0.1375, 0.14, 0.1425, 0.145, 0.1475, 0.15, 0.1525, 0.155, 0.1575, \
            0.16, 0.1625, 0.165, 0.1675, 0.17, 0.1725, 0.175, 0.1775, 0.18, \
            0.1825, 0.185, 0.1875, 0.19, 0.1925, 0.195, 0.1975, 0.2, 0.2025, \
            0.205, 0.2075, 0.21, 0.2125, 0.215, 0.2175, 0.22, 0.2225, 0.225, \
            0.2275, 0.23, 0.2325, 0.235, 0.2375, 0.24, 0.2425, 0.245, 0.2475, \
            0.25, 0.2525, 0.255, 0.2575, 0.26, 0.2625, 0.265, 0.2675, 0.27, \
            0.2725, 0.275, 0.2775, 0.28, 0.2825, 0.285, 0.2875, 0.29, 0.2925, \
            0.295, 0.2975, 0.3, 0.3025, 0.305, 0.3075, 0.31, 0.3125, 0.315, \
            0.3175, 0.32, 0.3225, 0.325, 0.3275, 0.33, 0.3325, 0.335, 0.3375, \
            0.34, 0.3425, 0.345, 0.3475, 0.35, 0.3525, 0.355, 0.3575, 0.36, \
            0.3625, 0.365, 0.3675, 0.37, 0.3725, 0.375, 0.3775, 0.38, 0.3825, \
            0.385, 0.3875, 0.39, 0.3925, 0.395, 0.3975, 0.4, 0.4025, 0.405, \
            0.4075, 0.41, 0.4125, 0.415, 0.4175, 0.42, 0.4225, 0.425, 0.4275, \
            0.43, 0.4325, 0.435, 0.4375, 0.44, 0.4425, 0.445, 0.4475, 0.45, \
            0.4525, 0.455, 0.4575, 0.46, 0.4625, 0.465, 0.4675, 0.47, 0.4725, \
            0.475, 0.4775, 0.48, 0.4825, 0.485, 0.4875, 0.49, 0.4925, 0.495, \
            0.4975, 0.5, 0.5025, 0.505, 0.5075, 0.51, 0.5125, 0.515, 0.5175, \
            0.52, 0.5225, 0.525, 0.5275, 0.53, 0.5325, 0.535, 0.5375, 0.54, \
            0.5425, 0.545, 0.5475, 0.55, 0.5525, 0.555, 0.5575, 0.56, 0.5625, \
            0.565, 0.5675, 0.57, 0.5725, 0.575, 0.5775, 0.58, 0.5825, 0.585, \
            0.5875, 0.59, 0.5925, 0.595, 0.5975, 0.6, 0.6025, 0.605, 0.6075, \
            0.61, 0.6125, 0.615, 0.6175, 0.62, 0.6225, 0.625, 0.6275, 0.63, \
            0.6325, 0.635, 0.6375, 0.64, 0.6425, 0.645, 0.6475, 0.65, 0.6525, \
            0.655, 0.6575, 0.66, 0.6625, 0.665, 0.6675, 0.67, 0.6725, 0.675, \
            0.6775, 0.68, 0.6825, 0.685, 0.6875, 0.69, 0.6925, 0.695, 0.6975, \
            0.7, 0.7025, 0.705, 0.7075, 0.71, 0.7125, 0.715, 0.7175, 0.72, \
            0.7225, 0.725, 0.7275, 0.73, 0.7325, 0.735, 0.7375, 0.74, 0.7425, \
            0.745, 0.7475, 0.75, 0.7525, 0.755, 0.7575, 0.76, 0.7625, 0.765, \
            0.7675, 0.77, 0.7725, 0.775, 0.7775, 0.78, 0.7825, 0.785, 0.7875, \
            0.79, 0.7925, 0.795, 0.7975, 0.8, 0.8025, 0.805, 0.8075, 0.81, \
            0.8125, 0.815, 0.8175, 0.82, 0.8225, 0.825, 0.8275, 0.83, 0.8325, \
            0.835, 0.8375, 0.84, 0.8425, 0.845, 0.8475, 0.85, 0.8525, 0.855, \
            0.8575, 0.86, 0.8625, 0.865, 0.8675, 0.87, 0.8725, 0.875, 0.8775, \
            0.88, 0.8825, 0.885, 0.8875, 0.89, 0.8925, 0.895, 0.8975, 0.9, \
            0.9025, 0.905, 0.9075, 0.91, 0.9125, 0.915, 0.9175, 0.92, 0.9225, \
            0.925, 0.9275, 0.93, 0.9325, 0.935, 0.9375, 0.94, 0.9425, 0.945, \
            0.9475, 0.95, 0.9525, 0.955, 0.9575, 0.96, 0.9625, 0.965, 0.9675, \
            0.97, 0.9725, 0.975, 0.9775, 0.98, 0.9825, 0.985, 0.9875, 0.99, \
            0.9925, 0.995, 0.9975, 1.)),
            decimal=4)

    for n_k in xrange(0, N_k - 1):
        for n_D in xrange(0, N_D - 1):

            # The following cases are for the edges of integration
            # See the yeow_code.nb for the original
            # Odd and evens should be switched with Mathematica code, python
            # has 0 indexing, mathematica has 1 indexing

            # 1
            if ((n_k % 2 == 1) &
                (xknot[n_k] - xl[n_D] >= 0.9999 * dxk) &
                (xu[n_D] - xknot[n_k] >= 0.9999 * dxk)):
                a = 4.0
            # 2
            elif ((n_k % 2 == 0) &
                  (xknot[n_k] - xl[n_D] > 2 * dxk) &
                  (xu[n_D] - xknot[n_k] > 2 * dxk)):
                a = 2.0
            # 3
            elif ((n_k % 2 == 1) &
                  (xknot[n_k] - xu[n_D] >= dxk)):
                a = 0.0
            # 4
            elif ((n_k % 2 == 0) &
                  (xknot[n_k] - xu[n_D] >= 2 * dxk)):
                a = 0.0
            # 5
            elif ((n_k % 2 == 1) &
                  (xl[n_D] - xknot[n_k] >= dxk)):
                print 'yes'
                a = 0.0
            # 6
            elif ((n_k % 2 == 0) &
                  (xl[n_D] - xknot[n_k] >= 2 * dxk)):
                a = 0.0
            # 7
            elif ((n_k % 2 == 1) &
                  (xu[n_D] - xknot[n_k - 1] > 0) &
                  (xknot[n_k + 1] - xu[n_D] <= 2.0 * dxk) &
                  (xknot[n_k] - xl[n_D] >= dxk)):
                a = 3.0 * coeff_2((xu[n_D] - xknot[n_k - 1]) / dxk)
            # 8
            elif ((n_k % 2 == 0) &
                  (xu[n_D] - xknot[n_k] >= 0) &
                  (xu[n_D] - xknot[n_k] <= 2 * dxk) &
                  (xknot[n_k] - xl[n_D] >= 2 * dxk)):
                a = 1.0 + 3.0 * coeff_1((xu[n_D] - xknot[n_k]) / dxk)
            # 9
            elif ((n_k % 2 == 0) &
                  (xknot[n_k] - xu[n_D] >= 0) &
                  (xknot[n_k] - xu[n_D] <= 2 * dxk) &
                  (xknot[n_k] - xl[n_D] >= 2 * dxk)):
                a = 3.0 * coeff_3((xu[n_D] - xknot[n_k - 1]) / dxk + 1)
            # 10
            elif ((n_k % 2 == 1) &
                  (xl[n_D] - xknot[n_k - 1] >= 0) &
                  (xl[n_D] - xknot[n_k + 1] <= 0) &
                  (xu[n_D] - xknot[n_k] >= dxk)):
                a = 4.0 - 3.0 * coeff_2((xl[n_D] - xknot[n_k - 1]) / dxk)
            # 11
            elif ((n_k % 2 == 1) &
                  (xu[n_D] - xknot[n_k - 1] <= 2.00001 * dxk) &
                  (xl[n_D] - xknot[n_k - 1] <= 2.00001 * dxk) &
                  (xl[n_D] - xknot[n_k - 1] >= -0.00001)):
                a = 3.0 * (coeff_2((xu[n_D] - xknot[n_k - 1]) / dxk) \
                           - coeff_2((xl[n_D] - xknot[n_k - 1]) / dxk))
            # 12
            elif ((n_k % 2 == 0) and
                  (xu[n_D] - xknot[n_k] < 2 * dxk) and
                  (xu[n_D] - xknot[n_k] > 0) and
                  (xl[n_D] - xknot[n_k])):
                a = 3.0 * (coeff_1((xu[n_D] - xknot[n_k]) / dxk) \
                           - coeff_1((xl[n_D] - xknot[n_k]) / dxk))
            # 13
            elif ((n_k % 2 == 0) &
                  (xknot[n_k] - xu[n_D] >= -0.000001) &
                  (xknot[n_k] - xu[n_D] < 2 * dxk) &
                  (xknot[n_k] - xl[n_D] <= 2 * dxk)):
                a = 3.0 * (coeff_3((xu[n_D] - xknot[n_k - 2]) / dxk) \
                           - coeff_3((xl[n_D] - xknot[n_k - 2]) / dxk))
            # 14
            elif ((n_k % 2 == 0) &
                  (xu[n_D] - xknot[n_k] > 0) &
                  (xu[n_D] - xknot[n_k] < 2 * dxk) &
                  (xknot[n_k] - xl[n_D] < 2 * dxk)):
                a = 3.0 * coeff_1((xu[n_D] - xknot[n_k]) / dxk) + \
                    3.0 * coeff_3((xl[n_D] - xknot[n_k - 2]) / dxk) + 1
            # 15
            elif ((n_k % 2 == 0) &
                  (xknot[n_k] - xl[n_D] > 0) &
                  (xknot[n_k] - xl[n_D] <= 2 * dxk) &
                  (xu[n_D] - xknot[n_k] >= 2 * dxk)):
                a = 1 + 1 + 3.0 * coeff_3(2 - (xknot[n_k] - xl[n_D]) / dxk)
            # 16
            elif ((n_k % 2 == 0) &
                  (xl[n_D] - xknot[n_k] >= 0) &
                  (xl[n_D] - xknot[n_k] <= 2 * dxk) &
                  (xu[n_D] - xknot[n_k] >= 2 * dxk)):
                a = 1 - 3.0 * coeff_1((xl[n_D] - xknot[n_k]) / dxk)

            trap_integ[n_D, n_k] = a

    if testing:
        test_trap_integ(trap_integ)

    trap_integ *= Delta

    return trap_integ

def test_trap_integ(trap_integ):

    import numpy as np

    np.testing.assert_array_almost_equal(trap_integ[0:3, :],
        np.array(((0., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (1., 4,
  2., 4, 2.04956, 3.77843, 0.314869, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0), (1., 4, 2., 4, 2, 4, 2, 4, 2, 4, 2.1516, 3.207, -0.0728863, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0))),
        decimal=4)

    np.testing.assert_array_almost_equal(trap_integ[-4:, :],
        np.array(
        ((1., 4, 2., 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, \
    2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, \
    4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, \
    2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, \
    4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, \
    2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, \
    4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, \
    2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, \
    4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, \
    2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, \
    4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, \
    2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, \
    4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, \
    2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, \
    4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, \
    2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, \
    4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, \
    2.07289, 0.793003, -0.151603, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (1., 4, \
    2., 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, \
    4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, \
    2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, \
    4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, \
    2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, \
    4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, \
    2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, \
    4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, \
    2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2,
      4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, \
    4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, \
    2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, \
    4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, \
    2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, \
    4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, \
    2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, \
    4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, \
    2, 4, 1.68513, 0.221574, -0.0495627, 0, 0, 0, 0), (1., 4, 2., 4, 2, \
    4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, \
    2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, \
    4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, \
    2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, \
    4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, \
    2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, \
    4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, \
    2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, \
    4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, \
    2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, \
    4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, \
    2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, \
    4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, \
    2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, \
    4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, \
    2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, \
    4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, \
    2, 4, 2., 4, 1.))
            ),
        decimal=4)

def coeff_1(al):

    return al - 3.0 * al**2 / 4.0 + al**3 / 6.0

def coeff_2(al):

    return al**2 - al**3 / 3.0

def coeff_3(al):

    return - al**2 / 4.0 + al**3 / 6.0

def construct_B(lam_C, lam_M):

    ''' Constructs an integral operator matrix.

    http://mathfaculty.fullerton.edu/mathews/n2003/SimpsonsRule2DMod.html

    '''

    import numpy as np

    N_D = lam_M.size
    N_k = lam_C.size
    lam_M = np.reshape(lam_M, (N_D, 1))
    lam_C = np.reshape(lam_C, (N_k, 1))
    Delta = lam_M[1, 0] - lam_M[0, 0]
    Delta = lam_C[1, 0] - lam_C[0, 0]
    lam_0 = lam_M[0, 0]

    if 0:
        B = np.matrix(np.zeros((N_D,N_k)))

        for i in xrange(1, N_D):
            B[i, :] += B[i-1, :]
            B[i, i-1] += (lam_M[i - 1, 0] - lam_0)**3
            B[i, i] += (lam_M[i, 0] - lam_0)**3
            #B[i, i-1] += (Delta)**3
            #B[i, i] += (Delta)**3

        B *= 0.5 / 6.0 * Delta
    elif 1:

        '''
        # Create Trapezoidal integration matrix
        trap_integ = np.matrix(np.zeros((N_D, N_k)))
        scale = (N_k) / (N_D / 1.0)

        for i in xrange(1, N_D):
            lam_pos = i * scale
            trap_integ[i, :lam_pos] = np.ones(lam_pos)

        print lam_pos, N_k
        trap_integ[:, 1:lam_pos:2] *= 4
        trap_integ[:, 2:lam_pos:2] *= 2
        print trap_integ[-10:-1]
        trap_integ *= Delta
        '''

        # Construct variables
        lam_M = np.squeeze(np.asarray(lam_M))
        lam_C = np.squeeze(np.asarray(lam_C))
        Delta = lam_C[1] - lam_C[0]
        lam_0 = lam_M[0]
        N_D = lam_M.size
        N_k = lam_C.size

        # Prep scalings of axes
        dxk = 1.0 / (N_k - 1.0)
        dxu = (lam_M - lam_0) / (lam_M[-1] - lam_0)
        xu = lam_M
        xknot = np.linspace(0,1,N_k)
        dxl = np.zeros((N_D))
        xl = dxl
        xu = dxu

        trap_integ = construct_trap_integ(lam_C, lam_M)

        C = np.matrix(np.zeros((N_D,N_k)))

        for i in xrange(0, N_D):
            for j in xrange(0, N_k):
                if dxu[i] >= xknot[j]:
                    C[i, j] = (dxu[i] - xknot[j])**3 / 6.0

        testing = False
        if testing:
            assert dxk == 0.0025
            #np.testing.assert_array_almost_equal(dxl,
            #    np.array(
            #    )
            #    decimal=4)

        lam_M = np.reshape(lam_M, (N_D, 1))
        lam_C = np.reshape(lam_C, (N_k, 1))

        for i in xrange(0, N_D):
            for j in xrange(0, N_k):
                if lam_M[i, 0] >= lam_C[j, 0]:
                    C[i, j] = (lam_M[i, 0] - lam_C[j, 0])**3 / 6.0 / 3.0

        B = np.multiply(trap_integ, C)

    return B

def construct_beta(lam_C):

    import numpy as np

    offset = 2

    ones_array = np.ones(lam_C.shape[0] + offset)
    beta_diag = np.diag(ones_array[:], 0) \
                - 2 * np.diag(ones_array[:-1], 1) \
                + np.diag(ones_array[:-2], 2)
    beta_diag = beta_diag[:lam_C.shape[0] - 2, :-2]
    beta_zeros = np.zeros((ones_array.shape[0] - offset - 2, 4))

    beta = np.matrix(np.hstack((beta_diag, beta_zeros)))
    bega_diag, beta_zeros = None, None

    testing = False
    if testing:
        test_beta(beta.T * beta)

    return beta

def test_beta(beta_inprod):


    import numpy as np

    np.testing.assert_array_almost_equal(beta_inprod[0:4, :],
        np.array(
    ((1, -2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (-2, 5, -4, 1, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0), (1, -4, 6, -4, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, \
    1, -4, 6, -4, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            ),
        decimal=4)

    np.testing.assert_array_almost_equal(beta_inprod[-5,:],
        np.array(
    ((0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 1, -2, 1, 0, 0, 0, 0),)
            ),
        decimal=4)

    assert beta_inprod.shape == (405, 405)


