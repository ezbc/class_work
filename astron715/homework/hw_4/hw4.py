#!/usr/bin/python

import warnings
warnings.filterwarnings('ignore')

''' Plotting Functions
'''

def plot_r_vs_m(df_m1, df_m5, filename=None, show=False):

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
    params = {#'backend': .png',
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

    # Temp gradient difference
    ax = fig.add_subplot(221)
    temp_diff_m1 = np.abs(df_m1['Radiative temperature gradie']) - \
                   np.abs(df_m1['Adiabatic temperature gradie'])
    temp_diff_m5 = np.abs(df_m5['Radiative temperature gradie']) - \
                   np.abs(df_m5['Adiabatic temperature gradie'])

    ax.plot(df_m1['Radius coordinate'], df_m1['Radiative temperature gradie'],
            linewidth=2,
            #marker='+',
            #markersize=4,
            color='r',
            label=r'Radiative',
            )
    ax.plot(df_m1['Radius coordinate'], df_m1['Adiabatic temperature gradie'],
            linewidth=2,
            linestyle='--',
            #marker='+',
            #markersize=4,
            color='k',
            label=r'Adiabatic',
            )

    ax.set_xlabel('$r$ [R$_\odot$]')
    ax.set_ylabel(r'$\nabla$T')
    #ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(loc='best')

    ax = fig.add_subplot(223)

    if 0:
        ax.plot(df_m5['Radius coordinate'], temp_diff_m5,
                linewidth=2,
                linestyle='--',
                #marker='+',
                #markersize=4,
                color='k',
                label=r'M = 5 M$_\odot$',
                )

    ax.plot(df_m5['Radius coordinate'], df_m5['Radiative temperature gradie'],
            linewidth=2,
            #marker='+',
            #markersize=4,
            color='r',
            label=r'Radiative',
            )

    ax.plot(df_m5['Radius coordinate'], df_m5['Adiabatic temperature gradie'],
            linewidth=2,
            linestyle='--',
            #marker='+',
            #markersize=4,
            color='k',
            label=r'Adiabatic',
            )

    ax.set_xlabel('$r$ [R$_\odot$]')
    ax.set_ylabel(r'$\nabla$T')
    #ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(loc='best')

    # r vs M
    # ------
    ax = fig.add_subplot(222)
    ax.plot(df_m1['Radius coordinate'], df_m1['Lagrangian mass coordinate'],
            linewidth=2,
            #marker='+',
            #markersize=4,
            color='r',
            label=r'M = 1 M$_\odot$',
            )

    ax.plot(df_m5['Radius coordinate'], df_m5['Lagrangian mass coordinate'],
            linewidth=2,
            linestyle='--',
            #marker='+',
            #markersize=4,
            color='k',
            label=r'M = 5 M$_\odot$',
            )

    ax.set_xlabel('$r$ [R$_\odot$]')
    ax.set_ylabel('$M_r$ [M$_\odot$]')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(loc='best')
    #ax.set_xlim([10**22.5, 10**24.5])
    #ax.set_ylim([10**16, 10**18])

    plt.tight_layout()

    if filename is not None:
        plt.savefig(filename, bbox_inches='tight', dpi=800)
    if show:
        plt.show()

def plot_pressure_vs_density(df_m1, df_m5, m1_fit=None, m5_fit=None,
        filename=None, show=False):

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
    color_cycle = [cmap(i) for i in np.linspace(0, 0.8, 5)]
    font_scale = 10
    line_weight = 800
    font_weight = 800
    params = {#'backend': .png',
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
              'figure.figsize': (4,4),
              'figure.titlesize': font_scale,
              'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)

    # Create figure instance
    fig = plt.figure()

    # r vs M
    # ------
    ax = fig.add_subplot(111)

    ax.plot(df_m1['Density'], df_m1['Total pressure'],
            linewidth=2,
            linestyle='--',
            #marker='+',
            #markersize=4,
            #color='r',
            label=r'M = 1 M$_\odot$',
            )

    ax.plot(df_m5['Density'], df_m5['Total pressure'],
            linewidth=2,
            linestyle='--',
            #marker='+',
            #markersize=4,
            #color='k',
            label=r'M = 5 M$_\odot$',
            )

    ax.plot(df_m1['Density'], m1_fit,
                #color='r',
                label='1 M$_\odot$ fit',
                linewidth=2,
                alpha=0.5,
                )
    ax.plot(df_m5['Density'], m5_fit,
                #color='r',
                label='5 M$_\odot$ fit',
                linewidth=2,
                alpha=0.5,
                )

    ax.set_xlabel('Density [g cm$^{-3}$]')
    ax.set_ylabel('Pressure [dyn cm$^{-2}$]')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(loc='best')
    #ax.set_xlim([10**22.5, 10**24.5])
    #ax.set_ylim([10**16, 10**18])

    plt.tight_layout()

    if filename is not None:
        plt.savefig(filename, bbox_inches='tight', dpi=800)
    if show:
        plt.show()

def plot_pressure_vs_density_norm(df_m1, df_m5, p_fit_m1=None, p_fit_m5=None,
        filename=None, show=False):

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
    color_cycle = [cmap(i) for i in np.linspace(0, 0.8, 5)]
    font_scale = 10
    line_weight = 800
    font_weight = 800
    params = {#'backend': .png',
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
              'figure.figsize': (4,4),
              'figure.titlesize': font_scale,
              'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)

    # Create figure instance
    fig = plt.figure()

    # r vs M
    # ------
    ax = fig.add_subplot(111)

    ax.plot(df_m1['Density'],
            df_m1['Total pressure'],# / np.max(df_m1['Total pressure']),
            linewidth=2,
            linestyle='--',
            #marker='+',
            #markersize=4,
            #color='r',
            label=r'M = 1 M$_\odot$',
            )

    ax.plot(df_m5['Density'],
            df_m5['Total pressure'],# / np.max(df_m5['Total pressure']),
            linewidth=2,
            linestyle='--',
            #marker='+',
            #markersize=4,
            #color='k',
            label=r'M = 5 M$_\odot$',
            )

    ax.plot(df_m1['Density'],
            p_fit_m1,
                #color='r',
                label='1 M$_\odot$ Adiabatic',
                linewidth=2,
                alpha=0.5,
                )
    ax.plot(df_m5['Density'], p_fit_m5,
                #color='r',
                label='5 M$_\odot$ Adiabatic',
                linewidth=2,
                alpha=0.5,
                )

    ax.set_xlabel(r'Normalized Density')
    ax.set_ylabel(r'Normalized Pressure')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(loc='best')
    #ax.set_xlim([10**22.5, 10**24.5])
    #ax.set_ylim([10**16, 10**18])

    plt.tight_layout()

    if filename is not None:
        plt.savefig(filename, bbox_inches='tight', dpi=800)
    if show:
        plt.show()

def plot_gamma(df_m1, df_m5, n_m1=None, n_m5=None,
        filename=None, show=False):

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
    color_cycle = [cmap(i) for i in np.linspace(0, 0.8, 5)]
    font_scale = 10
    line_weight = 800
    font_weight = 800
    params = {#'backend': .png',
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
              'figure.figsize': (4,4),
              'figure.titlesize': font_scale,
              'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)

    # Create figure instance
    fig = plt.figure()

    # r vs M
    # ------
    ax = fig.add_subplot(111)

    gamma_m1 = df_m1['First adiabatic expone']
    gamma_m5 = df_m5['First adiabatic expone']

    ax.plot(df_m5['Radius coordinate'],
            gamma_m1 - (n_m1 + 1) / n_m1,# / np.max(df_m1['Total pressure']),
            linewidth=2,
            linestyle='-',
            #marker='+',
            #markersize=4,
            #color='r',
            label=r'M = 1 M$_\odot$',
            )

    ax.plot(df_m5['Radius coordinate'],
            gamma_m5 - (n_m5 + 1) / n_m5,# / np.max(df_m1['Total pressure']),
            linewidth=2,
            linestyle='--',
            #marker='+',
            #markersize=4,
            #color='r',
            label=r'M = 5 M$_\odot$',
            )

    ax.set_xlabel(r'Radius [M$_\odot$]')
    ax.set_ylabel(r'$\Gamma_1$ - (n + 1) / n')
    ax.set_xscale('log')
    #ax.set_yscale('log')
    ax.legend(loc='best')
    #ax.set_xlim([10**0.5, 10**-2])
    #ax.set_ylim([10**16, 10**18])

    plt.tight_layout()

    if filename is not None:
        plt.savefig(filename, bbox_inches='tight', dpi=800)
    if show:
        plt.show()

def fit_pressure(density, pressure, return_params=False, p0=None):

    from scipy.optimize import curve_fit
    import numpy as np

    # fit power law
    calc_y = lambda x, K, n: K * (1 + n) / n * np.log(density) + np.log(K)

    if 0:
        fit_params = curve_fit(calc_y,
                               density,
                               pressure,
                               p0=p0,
                               maxfev=100000)[0]

    results = np.polyfit(np.log(density), np.log(pressure), deg=1)

    fit_y = np.exp(results[1])*density**results[0]

    K = np.exp(results[1])

    n = 1.0 / (results[0] - 1.0)

    print('n = {0:e}'.format(n))
    print('K = {0:e}'.format(K))

    if return_params:
        return fit_y, n, K
    else:
        return fit_y

def calc_broke_y(x, alpha1, a1, beta1, alpha2, a2, beta2, x0):

    import numpy as np

    y = np.zeros(len(x))
    y[np.where(x < x0)[0]] = beta1 * x[np.where(x < x0)[0]]**alpha1 + a1
    y[np.where(x >= x0)[0]] = beta2 * x[np.where(x >= x0)[0]]**alpha2 + a2

    return y

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
    df_m1 = pd.read_csv('m1_structure.txt', delimiter='\t')
    df_m5 = pd.read_csv('m5_structure.txt', delimiter='\t')

    cols = df_m1.keys()
    new_cols = {}
    for col in cols:
        new_cols[col] = col[0:col.find('(')-1]

    df_m1.rename(columns=new_cols, inplace=True)
    df_m5.rename(columns=new_cols, inplace=True)

    plot_r_vs_m(df_m1, df_m5,
                filename='fig_1a.png',
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
    df_m1 = pd.read_csv('m1_structure.txt', delimiter='\t')
    df_m5 = pd.read_csv('m5_structure.txt', delimiter='\t')

    cols = df_m1.keys()
    new_cols = {}
    for col in cols:
        new_cols[col] = col[0:col.find('(')-1]

    df_m1.rename(columns=new_cols, inplace=True)
    df_m5.rename(columns=new_cols, inplace=True)

    print('\nMsun = 1')
    m1_fit = fit_pressure(df_m1['Density'], df_m1['Total pressure'])

    print('\nMsun = 5')
    m5_fit = fit_pressure(df_m5['Density'], df_m5['Total pressure'])

    plot_pressure_vs_density(df_m1,
                             df_m5,
                             m1_fit=m1_fit,
                             m5_fit=m5_fit,
                             filename='fig_1b.png',
                             )

def prob1c():

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
    df_m1 = pd.read_csv('m1_structure.txt', delimiter='\t')
    df_m5 = pd.read_csv('m5_structure.txt', delimiter='\t')

    cols = df_m1.keys()
    new_cols = {}
    for col in cols:
        new_cols[col] = col[0:col.find('(')-1]

    df_m1.rename(columns=new_cols, inplace=True)
    df_m5.rename(columns=new_cols, inplace=True)

    n = 1.0 / (1.1 - 1)

    #p_fit_m1 = df_m1['Density']**df_m1['First adiabatic expone']
    #p_fit_m5 = df_m5['Density']**df_m5['First adiabatic expone']
    p_fit_m1 = df_m1['Density']*np.exp(df_m1['First adiabatic expone'])
    p_fit_m5 = df_m5['Density']*np.exp(df_m5['First adiabatic expone'])
    #p_fit_m1 /= np.max(p_fit_m1)
    #p_fit_m5 /= np.max(p_fit_m5)

    print('\nMsun = 1')
    m1_fit, n_m1, k = fit_pressure(df_m1['Density'], df_m1['Total pressure'],
                          return_params=True)

    print('\nMsun = 5')
    m5_fit, n_m5, k = fit_pressure(df_m5['Density'], df_m5['Total pressure'],
                          return_params=True)

    if 0:
        plot_pressure_vs_density_norm(df_m1,
                                 df_m5,
                                 p_fit_m1=p_fit_m1,
                                 p_fit_m5=p_fit_m5,
                                 filename='fig_1c',
                                 )

    plot_gamma(df_m1,
               df_m5,
               n_m1,
               n_m5,
               filename='fig_1c.png',
               )

def prob2():

    G = 6.67e-11
    Msun = 1.9e30
    Lsun = 3.846e26
    Rsun = 6.9e8

    t = 3 * G * Msun**2 / Lsun / Rsun

    print(t / (365.0 * 24 * 3600.0))

def main():

    prob1a()
    prob1b()
    prob1c()
    #prob2()

if __name__ == '__main__':

    main()

