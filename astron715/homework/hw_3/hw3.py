#!/usr/bin/python

import warnings
warnings.filterwarnings('ignore')

''' Plotting Functions
'''

def plot_data(data, L_fit=None, L_broken_fit=None, R_fit=None,
        R_broken_fit=None, T_fit=None, T_broken_fit=None, filename=None,
        show=False):

    # Import external modules
    import matplotlib.pyplot as plt
    import matplotlib
    import numpy as np

    # Set up plot aesthetics
    # ----------------------
    plt.close;plt.clf()
    plt.rcdefaults()

    # Color map
    cmap = plt.cm.gnuplot

    # Color cycle, grabs colors from cmap
    color_cycle = [cmap(i) for i in np.linspace(0, 0.8, 3)]
    font_scale = 10
    line_weight = 800
    font_weight = 800
    params = {#'backend': .pdf',
              'axes.labelsize': font_scale,
              'axes.titlesize': font_scale,
              'axes.weight': line_weight,
              'text.fontsize': font_scale,
              'legend.fontsize': font_scale*3/4,
              'xtick.labelsize': font_scale,
              'xtick.weight': line_weight,
              'ytick.labelsize': font_scale,
              'ytick.weight': line_weight,
              'font.weight': font_weight,
              'axes.labelweight': font_weight,
              'text.usetex': True,
              #'font.family': 'sans-serif',
              'figure.figsize': (7, 7),
              'figure.titlesize': font_scale,
              'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)

    # Create figure instance
    fig = plt.figure()

    # L vs M
    # ------
    ax = fig.add_subplot(221)

    ax.plot(data['mass'], data['luminosity'],
            linewidth=0,
            marker='+',
            markersize=4,
            color='k',
            label='Model',
            )

    ax.plot(data['mass'], L_fit,
            #color='r',
            label='Power-law fit',
            alpha=0.5,
            )

    ax.plot(data['mass'], L_broken_fit,
            #color='g',
            label='Broken power-law fit',
            alpha=0.5,
            )

    ax.set_xlabel('Mass [$M_\odot$]')
    ax.set_ylabel('Luminosity [$L_\odot$]')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(loc='best')
    ax.set_xlim([10**-0.5, 10**2])
    ax.set_ylim([10**-1.8, 10**1])

    # R vs M
    # ------
    ax = fig.add_subplot(222)

    ax.plot(data['mass'], data['radius'],
            linewidth=0,
            marker='+',
            markersize=4,
            color='k',
            label='Model',
            )

    ax.plot(data['mass'], R_fit,
            #color='r',
            label='Power-law fit',
            alpha=0.5,
            )

    ax.plot(data['mass'], R_broken_fit,
            #color='g',
            label='Broken power-law fit',
            alpha=0.5,
            )

    ax.set_xlabel('Mass [$M_\odot$]')
    ax.set_ylabel('Radius [$R_\odot$]')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(loc='best')
    ax.set_xlim([10**-0.5, 10**2])
    ax.set_ylim([10**-1.8, 10**0.1])

    # Tc vs M
    # ------
    ax = fig.add_subplot(223)

    ax.plot(data['mass'], 10**data['T_center'],
            linewidth=0,
            marker='+',
            markersize=4,
            color='k',
            label='Model',
            )

    # fit power law
    calc_y = lambda M, alpha, a, b: b * M**alpha + a

    ax.plot(data['mass'], T_fit,
            #color='r',
            label='Power-law fit',
            alpha=0.5,
            )

    ax.plot(data['mass'], T_broken_fit,
            #color='g',
            label='Broken power-law fit',
            alpha=0.5,
            )

    ax.set_xlabel('Mass [$M_\odot$]')
    ax.set_ylabel('T$_c$ [K]')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(loc='best')
    ax.set_xlim([10**-0.5, 10**2])
    ax.set_ylim([10**6.5, 10**7.8])

    plt.tight_layout()

    if filename is not None:
        plt.savefig(filename, bbox_inches='tight')
    if show:
        plt.show()

def plot_pressure(x, y, filename=None, show=False):

    # Import external modules
    import matplotlib.pyplot as plt
    import matplotlib
    import numpy as np

    # Set up plot aesthetics
    # ----------------------
    plt.close;plt.clf()
    plt.rcdefaults()

    # Color map
    cmap = plt.cm.gnuplot

    # Color cycle, grabs colors from cmap
    color_cycle = [cmap(i) for i in np.linspace(0, 0.8, 3)]
    font_scale = 10
    line_weight = 800
    font_weight = 800
    params = {#'backend': .pdf',
              'axes.labelsize': font_scale,
              'axes.titlesize': font_scale,
              'axes.weight': line_weight,
              'text.fontsize': font_scale,
              'legend.fontsize': font_scale*3/4,
              'xtick.labelsize': font_scale,
              'xtick.weight': line_weight,
              'ytick.labelsize': font_scale,
              'ytick.weight': line_weight,
              'font.weight': font_weight,
              'axes.labelweight': font_weight,
              'text.usetex': True,
              #'font.family': 'sans-serif',
              'figure.figsize': (3, 3),
              'figure.titlesize': font_scale,
              'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)

    # Create figure instance
    fig = plt.figure()

    # L vs M
    # ------
    ax = fig.add_subplot(111)

    ax.plot(x, y,
            linewidth=0,
            marker='+',
            markersize=4,
            color='k',
            label='Model',
            )

    if 0:
        ax.plot((10^-10, 10^-10), (10^100, 10^100),
                color='k',
                label='Linear relationship',
                linestyle='-',
                linewidth=2)

    if 0:
        ax.plot(x, y_fit,
                #color='r',
                label='Power-law fit',
                alpha=0.5,
                )

    ax.set_xlabel('M$^2$ / R$^4$ [erg cm$^{-3}$]')
    ax.set_ylabel('P$_c$ [erg cm$^{-3}$]')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(loc='best')
    ax.set_xlim([10**22.5, 10**24.5])
    ax.set_ylim([10**16, 10**18])

    plt.tight_layout()

    if filename is not None:
        plt.savefig(filename, bbox_inches='tight')
    if show:
        plt.show()

def plot_opacities(mass, k_bf, k_ff, k_e, filename=None, show=False):

    # Import external modules
    import matplotlib.pyplot as plt
    import matplotlib
    import numpy as np

    # Set up plot aesthetics
    # ----------------------
    plt.close;plt.clf()
    plt.rcdefaults()

    # Color map
    cmap = plt.cm.gnuplot

    # Color cycle, grabs colors from cmap
    color_cycle = [cmap(i) for i in np.linspace(0, 0.8, 3)]
    font_scale = 10
    line_weight = 800
    font_weight = 800
    params = {#'backend': .pdf',
              'axes.labelsize': font_scale,
              'axes.titlesize': font_scale,
              'axes.weight': line_weight,
              'text.fontsize': font_scale,
              'legend.fontsize': font_scale*3/4,
              'xtick.labelsize': font_scale,
              'xtick.weight': line_weight,
              'ytick.labelsize': font_scale,
              'ytick.weight': line_weight,
              'font.weight': font_weight,
              'axes.labelweight': font_weight,
              'text.usetex': True,
              #'font.family': 'sans-serif',
              'figure.figsize': (3, 3),
              'figure.titlesize': font_scale,
              'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)

    # Create figure instance
    fig = plt.figure()

    # L vs M
    # ------
    ax = fig.add_subplot(111)

    ax.plot(mass, k_bf,
            linewidth=0,
            marker='+',
            markersize=3,
            label='Bound-free',
            )
    ax.plot(mass, k_ff,
            linewidth=0,
            marker='+',
            markersize=3,
            label='Free-free',
            )
    ax.plot(mass, k_e,
            linewidth=0,
            marker='+',
            markersize=3,
            label='e$^-$',
            )

    ax.set_xlabel('Mass [M$_\odot$]')
    ax.set_ylabel('Opacity [cm$^2$ g$^{-1}$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(loc='best')
    #ax.set_xlim([10**-0.5, 10**2])
    #ax.set_ylim([10**-1.8, 10**1])

    plt.tight_layout()

    if filename is not None:
        plt.savefig(filename, bbox_inches='tight')
    if show:
        plt.show()

def fit_power_law(x, y, p0=None):

    from scipy.optimize import curve_fit

    # fit power law
    calc_y = lambda x, alpha, a, b: b * x**alpha + a

    fit_params = curve_fit(calc_y,
                           x,
                           y,
                           p0=p0,
                           maxfev=100000)[0]

    fit_y = calc_y(x, *fit_params)

    return fit_y

def calc_broke_y(x, alpha1, a1, beta1, alpha2, a2, beta2, x0):

    import numpy as np

    y = np.zeros(len(x))
    y[np.where(x < x0)[0]] = beta1 * x[np.where(x < x0)[0]]**alpha1 + a1
    y[np.where(x >= x0)[0]] = beta2 * x[np.where(x >= x0)[0]]**alpha2 + a2

    return y

def fit_broken_power_law(x, y, p0=None):

    from scipy.optimize import curve_fit

    # fit power law

    fit_params = curve_fit(calc_broke_y,
                           x,
                           y,
                           p0=p0,
                           maxfev=100000)[0]

    fit_y = calc_broke_y(x, *fit_params)

    return fit_y

def fit_slope(x, y, p0=None):

    from scipy.optimize import curve_fit

    # fit power law
    calc_y = lambda x, m: m * x

    fit_params = curve_fit(calc_y,
                           x,
                           y,
                           p0=p0,
                           maxfev=100000)[0]

    print('P = {0:e} M^2 / R^4 erg cm^-3'.format(fit_params[0]))

    fit_y = calc_y(x, *fit_params)

    return fit_y

def prob1a():

    ''' Executes script.
    '''

    # import external modules
    import numpy as np
    import os
    import pandas as pd
    import h5py

    # Script parameters
    # -----------------


    # Analysis
    # --------
    df = pd.read_csv('stars.csv')

    cols = df.keys()
    new_cols = {}
    for col in cols:
        new_cols[col] = col[0:col.find('(')]

    df.rename(columns=new_cols, inplace=True)

    # Fit data
    L_fit = fit_power_law(df['mass'],
                          df['luminosity'],
                          p0=[4.51275795e-05,
                              -3.61314583e+04,
                              3.61313727e+04],
                          )

    L_broken_fit = fit_broken_power_law(df['mass'],
                                        df['luminosity'],
                                        p0=[-6.62937152e-05,
                                            2.69226500e+04,
                                            -2.69227557e+04,
                                            3.09578562e-04,
                                            -4.24308204e+03,
                                            4.24377163e+03,
                                            4.00000000e+00],
                                        )
    R_fit = fit_power_law(df['mass'],
                          df['radius'],
                          )

    R_broken_fit = fit_broken_power_law(df['mass'],
                                        df['radius'],
                                        )

    T_fit = fit_power_law(df['mass'],
                          df['T_center'],
                          )

    T_broken_fit = fit_broken_power_law(df['mass'],
                                        df['T_center'],
                                        )

    plot_data(df,
              filename='fig_1a',
              L_fit=L_fit,
              L_broken_fit=L_broken_fit,
              R_fit=R_fit,
              R_broken_fit=R_broken_fit,
              T_fit=T_fit,
              T_broken_fit=T_broken_fit,
              )

def prob1b():

    ''' Executes script.
    '''

    # import external modules
    import numpy as np
    import os
    import pandas as pd
    import h5py

    # Script parameters
    # -----------------


    # Analysis
    # --------
    df = pd.read_csv('stars.csv')

    cols = df.keys()
    new_cols = {}
    for col in cols:
        new_cols[col] = col[0:col.find('(')]

    df.rename(columns=new_cols, inplace=True)

    mass = df['mass'] * 1.9e33 # g
    radius = 10**df['radius'] * 6.955e10 # cm

    Pc_fit = fit_slope(mass**2 / radius**4,
                       10**df['pressure'],
                       p0=1)

    plot_pressure(mass**2 / radius**4, 10**df['pressure'],
                  filename='fig_1b.png')

    G = 6.67e-11
    Msun = 1.9e30
    Lsun = 3.846e26
    Rsun = 6.9e8

    print('1b constant = {0:e}'.format(-3*G/(2*np.pi)))

def prob1c():

    # import external modules
    import numpy as np
    import os
    import pandas as pd
    import h5py

    # Script parameters
    # -----------------


    # Analysis
    # --------
    df = pd.read_csv('stars.csv')

    cols = df.keys()
    new_cols = {}
    for col in cols:
        if '(' in col:
            new_cols[col] = col[0:col.find('(')]
        else:
            new_cols[col] = col

    df.rename(columns=new_cols, inplace=True)

    Z = 0.02
    X = df['X']
    Y = df['Y']
    mass = df['mass'] * 1.9e33
    radius = df['radius'] * 6.955e8
    density = df['density']
    T = df['T_center']

    k_bf = 4e25 * Z * (1 + X) * density * T**-3.5 # cm^2 g^-1

    k_ff = 4e22 * (X + Y) * (1 + X) * density * T**-3.5 # cm^2 g^-1

    k_e = 0.2 * (1 + X)

    plot_opacities(df['mass'], k_bf, k_ff, k_e,
                  filename='fig_1c.png')

def prob2():

    G = 6.67e-11
    Msun = 1.9e30
    Lsun = 3.846e26
    Rsun = 6.9e8

    t = 3 * G * Msun**2 / Lsun / Rsun

    print(t / (365.0 * 24 * 3600.0))

def main():

    #prob1a()
    prob1b()
    prob1c()
    prob2()

if __name__ == '__main__':

    main()

