#!/usr/bin/python


import numpy as np
from scipy.optimize import curve_fit



def problem_1():

    q = 4/9.
    q_0 = 0.2

    i_rad = np.arccos(((q**2 - q_0**2) / (1 - q_0**2))**0.5)

    i_deg = np.degrees(i_rad)

    print('i = %.2f deg' % i_deg)

    u_obs = np.array([20.1, 21.1, 22.7])
    r_arcmin = np.array([1, 2, 8]) # arcmin

    distance = 9e3 # kpc

    r_kpc = r_arcmin / (206265 / 1e3) * distance
    u_real = u_obs - 2.5*np.log10(9/4.)

    sigma = (10**(4.74 + 21.572 - u_real) / 2.5) / (1e3**2)

    def profile(r, rd, sigma_0):
        return sigma_0 * np.exp(-r / rd)

    rd, sigma_0 = curve_fit(profile, r_kpc, sigma,
            maxfev=1000, p0=(100, 10))[0]

    print('rd = %.2f kpc' % rd)
    print('sigma_0 = %.2f L_r,sun kpc**-2' % sigma_0)

def main():
    problem_1()

if __name__ == '__main__':
    main()

