def prob_1():

    import numpy as np
    from constants import cgs

    F = 10**((-29 + 0.03) / 2.5)

    F = 3640 * 10**(-0.4 * (29 - 0.3))

    photon_flux = F * 5500 / (cgs.h * cgs.c)

    print('F = {0:.2g} erg s^-1 cm^-2 Hz^-1'.format(F))
    print('F_gamma = {0:.3g}'.format(photon_flux))


def main():
    prob_1()

if __name__ == '__main__':
    main()
