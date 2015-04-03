#!/usr/bin/python

import warnings
warnings.filterwarnings('ignore')

def plot_L_vs_T(df, filename=None, show=False, limits=None,
        scale=['linear','linear']):

    # Import external modules
    import numpy as np
    import math
    import pyfits as pf
    import matplotlib.pyplot as plt
    import matplotlib
    from mpl_toolkits.axes_grid1 import ImageGrid
    from astroML.plotting import scatter_contour
    from itertools import cycle
    from matplotlib.ticker import FormatStrFormatter
    from matplotlib.ticker import ScalarFormatter
    import myplotting as myplt

    # Set up plot aesthetics
    # ----------------------
    plt.close;plt.clf()
    plt.rcdefaults()

    # Color map
    cmap = plt.cm.gnuplot

    # Color cycle, grabs colors from cmap
    color_cycle = [cmap(i) for i in np.linspace(0, 0.8, 5)]
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

    lines = ["-","--","-.",":"]
    linecycler = cycle(lines)

    for i, Z in enumerate(df):
        ax.plot(df[Z]['Teff [K]'],
                df[Z]['Luminosity [Lsun]'],
                #linestyle=next(linecycler),
                label=r'Z = {0:s} Z$_\odot$'.format(Z)
                )

    if limits is not None:
        ax.set_xlim(limits[0],limits[1])
        ax.set_ylim(limits[2],limits[3])


    # set scales, linear or log
    ax.set_xscale(scale[0])
    ax.set_yscale(scale[1])

    # change labels to be integer scalars
    fig.canvas.draw()
    if 1:
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.xaxis.set_minor_formatter(ScalarFormatter())

    #ax.xaxis.set_minor_formatter(FormatStrFormatter('%.1f'))

    # Make sure no tick labels are overplotting one another
    ax = myplt.delete_overlapping_xlabels(fig, ax)

    # Adjust asthetics
    ax.set_xlabel(r'$T_{\rm eff}$ [K]')
    ax.set_ylabel(r'L [L$_\odot$]')
    ax.legend(loc='best')

    #ax.locator_params(nbins=5, axis='x')

    if filename is not None:
        plt.tight_layout()
        plt.savefig(filename, bbox_inches='tight')
    if show:
        fig.show()

def prob1a():

    import numpy as np
    import pandas as pd
    from scipy.integrate import simps as integrate
    from scipy.integrate import cumtrapz
    import pickle

    print('\n---------------')
    print('1d')
    print('---------------\n')

    with open('evol_data.pkl', 'rb') as handle:
        df = pickle.load(handle)

    plot_L_vs_T(df,
                limits=[15e3, 2e3, 0.1, 5e4],
                scale=['log', 'log'],
                filename='fig1a.png')
    if 0:
        plot_L_vs_T(df,
                    limits=[150e3, 2e3, 0.1, 5e4],
                    scale=['log', 'log'],
                    filename='fig1a.png')

    if 0:
        plot_L_vs_T(df,
                    limits=[15e3, 2e3, 0.1, 5e4],
                    scale=['log', 'log'],
                    filename='fig1a.pdf')


def main():

    prob1a()
    #prob1b()
    #prob1c()
    #prob1d()

if __name__ == '__main__':
    main()
