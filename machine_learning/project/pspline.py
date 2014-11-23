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
    for chi in chis:
        A_C, h, coeffs, lam_C, V = calc_V(x, y, chi, N_k=N_k)

        V_list.append(V)
        h_list.append(h)
        A_C_list.append(A_C)


    import matplotlib.pyplot as plt
    plt.clf(); plt.close()
    fig = plt.figure(figsize=(12,8))
    n = np.ceil(np.sqrt(len(h_list)))
    lam_C = np.squeeze(np.asarray(lam_C))
    for i in xrange(0, len(h_list)):
        h_hat = np.squeeze(np.asarray(h_list[i]))
        ax = fig.add_subplot(3, 4, i)
        ax.plot(lam_C, h_hat)

    plt.tight_layout()
    fig.savefig('4th_derivs.png')

    h_hat = h_list[np.where(V_list == min(V_list))[0]]
    A_C_hat = A_C_list[np.where(V_list == min(V_list))[0]]

    return A_C_hat, h_hat, lam_C, V_list

def calc_V(x, y, chi, N_k=None):

    '''

    A_C = computed spectrum [N_D x 1]
    B = coefficients [N_D x N_k]
    ones = column vector of ones [N_D x 1]
    lambda_M = measured wavelengths [N_D x 1]
    lambda_C = computed wavelengths [N_k x 1]

    '''

    import numpy as np
    import scipy
    from scipy import linalg

    N_D = len(x)

    # Subscript C = computed spectrum, sample length N_D
    # Subscript M = measured spectrum, sample length N_k

    # initialize matrices
    A_M, lam_M, lam_C, N_k = prep_spectrum(x, y, N_k=N_k)

    assert A_M.shape == (N_D, 1)
    assert lam_M.shape == (N_D, 1)
    assert lam_C.shape == (N_k, 1)

    # Reference wavelength
    lam_0 = lam_M[0, 0]
    Delta = abs(lam_C[-1, 0] - lam_C[0, 0]) / (N_k - 1)

    # Ones is a column vector of ones of shape N_D x 1
    ones = np.matrix(np.ones(A_M.shape[0])).T
    lam_M0 = lam_M - ones * lam_0
    assert ones.shape == lam_M0.shape

    B = construct_B(N_D, N_k, Delta, lam_C, lam_0)
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
    assert beta.shape == (N_k - 2, N_k + 4)

    # B_prime^T * B_prime --> N_k+4 x N_D * N_D x N_k+4 --> N_k+4 x N_k+4
    # beta^T * beta --> N_k+4 x N_k-2 * N_k-2 x N_k+4 --> N_k+4 x N_k+4

    # Write matrix shared by h_prime and E
    h_E_matrix = linalg.inv(B_prime.T * B_prime + chi / Delta**4 * \
                            beta.T * beta) * B_prime.T

    # Solve for h_prime, eq 4
    h_prime = h_E_matrix * A_M
    assert h_prime.shape == (N_k+4, 1)

    # Eq 6, N_D x N_k+4 * N_k+4 x 1= 1 x N_D
    A_C = B_prime * h_prime
    assert A_C.shape == A_M.shape

    # Eq 8
    E = B_prime * h_E_matrix

    # The optimum chi comes from minimizing V(chi)
    V = ((A_M - A_C).T * (A_M - A_C) / N_D) \
        / (1 - np.trace(E) / N_D)**2

    coeffs = h_prime[-4:]
    h_hat = h_prime[0:N_k]

    return A_C, h_hat, coeffs, lam_C, V[0,0]

def construct_B(N_D, N_k, Delta, lam_C, lam_0):

    ''' Constructs an integral operator matrix.

    http://mathfaculty.fullerton.edu/mathews/n2003/SimpsonsRule2DMod.html

    '''

    import numpy as np

    B = np.matrix(np.zeros((N_D,N_k)))
    if 1:
        for i in xrange(1, N_D):
            B[:, i] += B[:, i-1]
            B[i-1, i] += 0.5*Delta * (lam_C[i-1] - lam_0)**3 / 6.0
            B[i, i] += 0.5*Delta * (lam_C[i] - lam_0)**3 / 6.0
        #B[0, :] = 0.5*Delta
        #B[1, :-1] *= 2
        #B[0,0] = 0.5*Delta
    else:
        B[0, :] = np.ones((N_k))
        B[0, 1:-1] *= 2.0
        B[-1, :] = B[0, :]

        for j in xrange(1, N_D- 1):
            B[j, :] = 2 * B[0, :]

        indices = np.tril_indices(N_D)
        B[indices] = 0.0

        B *= Delta

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

    return beta



