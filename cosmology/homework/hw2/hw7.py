def t(z):

    import numpy as np

    omega_lam = 0.7
    omega_mat = 0.3
    H_0 = 2.28e-18 # s, #70.4 km/s Mpc

    t = 2/3. * 1 / (H_0 * omega_lam**0.5) *\
            np.log(((omega_lam/omega_mat) * (1+z)**-3)**0.5 + \
                   (1 + (omega_lam / omega_mat) * (1 + z)**-3)**0.5)

    return t

ratio =  t(10.0) / t(0.0)
print 'ratio = ', ratio
print 'time = ', ratio * 13.7, 'Gyr'

