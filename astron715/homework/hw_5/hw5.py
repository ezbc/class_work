#!/usr/bin/python

import warnings
warnings.filterwarnings('ignore')

''' Plotting Functions
'''

def plot_T_vs_density(density, T, filename=None, show=False):

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

    ax.plot(density, T,
            linewidth=2,
            linestyle='-',
            #marker='+',
            #markersize=4,
            #color='r',
            label=r'M = 1 M$_\odot$',
            )

    ax.annotate('Degenerate',
                xy=(0.5,0.1),
                xycoords='axes fraction',
                fontsize=font_scale)
    ax.annotate('Non-Degenerate',
                xy=(0.1,0.75),
                xycoords='axes fraction',
                fontsize=font_scale)

    ax.set_xlabel(r'$\rho / \mu_n$ [cm$^{-3}$]')
    ax.set_ylabel('T [K]')
    ax.set_xscale('log')
    ax.set_yscale('log')
    #ax.legend(loc='best')
    #ax.set_xlim([10**22.5, 10**24.5])
    #ax.set_ylim([10**16, 10**18])

    plt.tight_layout()

    if filename is not None:
        plt.savefig(filename, bbox_inches='tight', dpi=800)
    if show:
        plt.show()

def prob1():

    ''' Executes script.
    '''

    # import external modules
    import numpy as np
    import os
    import pandas as pd
    import h5py
    from constants import cgs

    # Script parameters
    # -----------------

    # Analysis
    # --------
    density_transition = 10e6 * cgs.me / cgs.mp
    density = np.logspace(0, 18, 100)
    T = np.zeros(density.shape)
    T_nonrel = (density * cgs.mp / cgs.me / 6.0e-9)**(1/1.5)
    T_rel = (density * cgs.mp / cgs.me / 4.6e-24)**(1/3.0)

    T[density < density_transition] = T_nonrel[density < density_transition]
    T[density > density_transition] = T_rel[density > density_transition]

    #density_nonrel = cgs.me / cgs.mp * 6.0e-9 * T**1.5 # g cm^-3
    #density_rel = cgs.me / cgs.mp * 4.6e-24 * T**3 # g cm^-3

    plot_T_vs_density(density, T,
                filename='fig_1.png',
                )

def prob2c():

    ''' Executes script.
    '''

    # import external modules
    import numpy as np
    import os
    import pandas as pd
    import h5py
    from constants import cgs

    # Script parameters
    # -----------------

    # Analysis
    # --------
    df = pd.read_csv('m1_structure.txt', delimiter='\t')
    cols = df.keys()
    new_cols = {}
    for col in cols:
        new_cols[col] = col[0:col.find('(')-1]

    df.rename(columns=new_cols, inplace=True)

    radius = np.argmin(np.abs(df['Radius coordinate'] - 0.62))
    X = df['Hydrogen mass fraction'][radius]

    n = df['Density'][radius] / cgs.mp * X
    T = df['Temperature'][radius]

    print('Density at 0.62 Rsun = {0:e} cm^-3'.format(n))
    print('Temp at 0.62 Rsun = {0:e} K'.format(T))

    NA = 6.0221413e+23
    T6 = T / 10**6
    sigma_v = 7.2*10**10 * T6**(-2/3.0) * np.exp(-84.7*T6**(-1/3.0)) / NA
    t = np.log(140.0)*(sigma_v * n)**-1

    t /= (365*24*3600.0)

    print('Time of depletion = {0:e} yr'.format(t))


def main():

    #prob1()
    prob2c()

if __name__ == '__main__':

    main()

