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
    A_M = np.matrix(y)
    lam_M = np.matrix(x)

    # Initialize model vectors, amps + wavelengths
    # --------------------------------------------
    # Define sampling size of model vectors
    if N_k is None:
        N_k = 100 * len(x)

    # wavelengths
    lam_C = np.matrix(np.ones((N_k)))

    # fourth derivative
    h = np.matrix(np.ones((N_k)))

    # amps
    A_C = np.matrix(np.ones((N_k)))

    return A_M, lam_M, lam_C, h, A_C

def fit_spline(x, y):

    import numpy as np
    from scipy.optimize import minimize

    result = minimize(calc_V, init_params, method='BFGS', options={'disp':True})

def derive_minimizer(x, y, h, r_0, f_0, g_0):

    import numpy as np
    import scipy
    from scipy import linalg

    def calc_V(params):
        h, r_0, f_0, g_0, chi = params

        # initialize matrices
        A_M, lam_M, lam_C, h, A_C = prep_spectrum(x, y)

        # Reference wavelength
        lam_0 = lam_M[0]
        A_0 = A_M[0]
        Delta = abs(lam_C[1] - lam_C[0])

        # Ones is a column vector of ones of shape N_D x 1
        ones = np.matrix(np.ones(len(A_M))).T

        B = integ_h(h, lam_M, lam_0)

        # See Eq 2b of Yeow et al. 2005
        A_C = 1/6.0 * B * h \
              + ones * A_0 \
              + (lam_M - ones * lam_0) * r_0 \
              + (lam_M - ones * lam_0)**2 / 2.0 * f_0 \
              + (lam_M - ones * lam_0)**3 / 6.0 * g_0

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

        #
        B_prime = np.hstack((B, A_0, r_0, f_0, g_0))
        h_prime = np.hstack((h, A_0, r_0, f_0, g_0))

        # Construct beta
        beta_diag = np.diag(ones, -1) - 2 * np.diag(ones,0) + np.diag(ones, 1)
        beta_zeros = np.zeros((len(ones), 4))
        beta = np.matrix(np.hstack((beta_diag, beta_zeros)))
        bega_diag, beta_zeros = None, None

        E = linalg.inv(B_prime.T * B_prime + chi / Delta**4 * beta.T * beta) \
                  * beta.T * A_M

        # The optimum chi comes from minimizing V(chi)
        V = ((A_M - A_C).T * (A_M - A_C) / N_D \
            / (1 - linalg.trace(E) / N_D)**2

        return V

    return calc_V

def integ_h(h, lam_M, lam_0):

    import numpy as np
    from scipy.integrate import quad

    def integrand(h, lam, lam_prime):
        return (lam - lam_prime)**3 * h

    integral = np.matrix(np.zeros((len(lam_M), len(h))))

    for i, lam in enumerate(lam_M):
        for j, lam_prime in enumerate(h):
            integral[i, j] = quad(integrand,
                                  lam_0,
                                  lam,
                                  args=(lam, lam_prime))[0]

    return integral






