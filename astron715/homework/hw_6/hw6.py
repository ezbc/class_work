#!/usr/bin/python

def prob1a():

    import numpy as np
    from scipy.integrate import simps as integrate

    print('\n---------------')
    print('1a')
    print('---------------\n')

    A = 2.73e-7
    A = 5.297e3
    Gamma = 1.3

    phi = lambda m: A * m**-(Gamma + 1.0)

    tau_ms = lambda m: m**-2.5 * 10**10

    masses = np.logspace(np.log10(0.5), np.log10(40), 20)

    N_m = phi(masses) #* tau_ms(masses)

    for mass in masses:
        print mass

    for i, n_m in enumerate(N_m):
        print('{1:.1f}\t{0:.0f}'.format(n_m, masses[i]))

    print('\nTotal # = {0:.0f}'.format(integrate(N_m, masses)))

def prob1b():

    import numpy as np
    import pandas as pd
    from scipy.integrate import simps as integrate
    from scipy.integrate import cumtrapz

    print('\n---------------')
    print('1b')
    print('---------------\n')

    df = pd.DataFrame.from_csv('cluster_data.csv', index_col=None)

    Ms = df['Mass [Msun]']
    Ls = df['Luminosity [Lsun]']
    Ts = df['Teff [K]']

    A = 5.297e3
    Gamma = 1.3
    phi = lambda m: A * m**-(Gamma + 1.0)
    N_m = phi(Ms) #* tau_ms(masses)

    L_bol = integrate(N_m * Ls, Ms)

    print('\nBolometric luminosity = {0:e} Lsun'.format(L_bol))


    L_cdf = cumtrapz(N_m * Ls, Ms)

    M_half = np.interp(0.5 * L_bol, L_cdf, Ms[1:])

    if 0:
        import matplotlib.pyplot as plt
        plt.close(); plt.clf()
        plt.plot(Ms[1:], L_cdf/L_cdf.max())
        plt.show()

    print('\nHalf light mass = {0:.2f} Msun'.format(M_half))

def prob1c():

    import numpy as np
    import pandas as pd
    from scipy.integrate import simps as integrate
    from scipy.integrate import cumtrapz

    print('\n---------------')
    print('1c')
    print('---------------\n')

    df = pd.DataFrame.from_csv('cluster_data.csv', index_col=None)

    Ms = df['Mass [Msun]']
    Ls = df['Luminosity [Lsun]']
    Ts = df['Teff [K]']

    A = 5.297e3
    Gamma = 1.3
    phi = lambda m: A * m**-(Gamma + 1.0)
    N_m = phi(Ms) #* tau_ms(masses)

    df_colors = pd.DataFrame.from_csv('colors.txt', index_col=None)

    color_table = df_colors['B-V']
    Teff_table = df_colors['Teff']

    color_table[color_table == '...'] = np.nan
    color_table = color_table.astype(float)

    print color_table

    colors = np.interp(Ts, Teff_table, color_table)

    print colors


def main():

    #prob1a()
    #prob1b()
    prob1c()

if __name__ == '__main__':
    main()
