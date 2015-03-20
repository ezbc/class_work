#!/usr/bin/python

import warnings
warnings.filterwarnings('ignore')

def plot_L_vs_T(T, L, masses=None, filename=None, show=False, limits=None,
        scale=['linear','linear']):

    # Import external modules
    import numpy as np
    import math
    import pyfits as pf
    import matplotlib.pyplot as plt
    import matplotlib
    from mpl_toolkits.axes_grid1 import ImageGrid
    from astroML.plotting import scatter_contour

    # Set up plot aesthetics
    # ----------------------
    plt.close;plt.clf()
    plt.rcdefaults()

    # Color map
    cmap = plt.cm.gnuplot

    # Color cycle, grabs colors from cmap
    color_cycle = [cmap(i) for i in np.linspace(0, 0.8, 3)]
    font_scale = 13
    line_weight = 600
    font_weight = 600

    params = {
              'axes.color_cycle': color_cycle, # colors of different plots
              'axes.labelsize': font_scale,
              'axes.titlesize': font_scale,
              #'axes.weight': line_weight,
              'axes.linewidth': 1.2,
              'axes.labelweight': font_weight,
              'legend.fontsize': font_scale*3/4,
              'xtick.labelsize': font_scale,
              'ytick.labelsize': font_scale,
              'font.weight': font_weight,
              'font.serif': 'computer modern roman',
              'text.fontsize': font_scale,
              'text.usetex': True,
              #'text.latex.preamble': r'\usepackage[T1]{fontenc}',
              #'font.family': 'sans-serif',
              'figure.figsize': (5, 5),
              'figure.dpi': 600,
              'backend' : 'pdf',
              #'figure.titlesize': font_scale,
             }
    plt.rcParams.update(params)

    # Create figure instance
    fig = plt.figure()

    ax = fig.add_subplot(111)


    bins = (T[::-1], L)

    bins = 15
    L = np.log10(L)
    scale = 'linear', 'linear'
    limits=[40000,500,-3,6]

    l1 = scatter_contour(T, L,
                         threshold=1,
                         log_counts=True,
                         levels=6,
                         ax=ax,
                         histogram2d_args=dict(bins=bins,
                                range=((limits[1], limits[0]),
                                       (limits[2], limits[3])),
                                ),
                         plot_args=dict(marker='o',
                                        linestyle='none',
                                        color='black',
                                        alpha=0.7,
                                        markersize=2),
                         contour_args=dict(
                                           cmap=plt.cm.gray_r,
                                           #cmap=cmap,
                                           ),
                                 )

    if limits is not None:
        ax.set_xlim(limits[0],limits[1])
        ax.set_ylim(limits[2],limits[3])

    ax.set_xscale(scale[0])
    ax.set_yscale(scale[1])


    # Adjust asthetics
    ax.set_xlabel(r'$T_{\rm eff}$ [K]')
    ax.set_ylabel(r'log$_{10}( L $ [L$_\odot$])')

    ax.locator_params(nbins=5, axis='x')

    if filename is not None:
        plt.savefig(filename, bbox_inches='tight')
    if show:
        fig.show()

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

    Teff_table = Teff_table[~np.isnan(color_table)]
    color_table = color_table[~np.isnan(color_table)]

    colors = np.interp(Ts, np.sort(Teff_table), np.sort(color_table))

    L_tot = integrate(N_m * Ls, Ms)
    N_tot = integrate(N_m, Ms)

    print('Ntot', N_tot)
    print colors
    print N_m

    color_dm = np.log10(2.512**(colors) * Ls/L_tot * N_m/N_tot) \
               / np.log10(2.512)

    color_dm = 2.512**(colors) * Ls/L_tot * N_m/N_tot

    avg_color = integrate(color_dm, Ms)

    print avg_color

    print('\nAverage color = {0:e} mag'.format(avg_color))

def prob1d():

    import numpy as np
    import pandas as pd
    from scipy.integrate import simps as integrate
    from scipy.integrate import cumtrapz

    print('\n---------------')
    print('1d')
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

    Teff_table = Teff_table[~np.isnan(color_table)]
    color_table = color_table[~np.isnan(color_table)]

    colors = np.interp(Ts, np.sort(Teff_table), np.sort(color_table))

    Ts_plot = []
    Ls_plot = []
    for i, T in enumerate(Ts):
        for j in xrange(int(N_m[i])):
            Ts_plot.append(T)
            Ls_plot.append(Ls[i])

    Ts_plot = np.array(Ts_plot)
    Ls_plot = np.array(Ls_plot)

    plot_L_vs_T(Ts_plot, Ls_plot,
                masses=Ms,
                limits=[40e3, 2e3, 0.001, 1e5],
                scale=['linear', 'log'],
                filename='fig1c.png')

    plot_L_vs_T(Ts_plot, Ls_plot,
                masses=Ms,
                limits=[40e3, 2e3, 0.001, 1e5],
                scale=['linear', 'log'],
                filename='fig1c.pdf')

def main():

    prob1a()
    prob1b()
    prob1c()
    prob1d()

if __name__ == '__main__':
    main()
