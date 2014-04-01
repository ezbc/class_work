#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import pyfits as pf

def plot_spectrum(wavelengths = None, magnitudes = None, limits = None, savedir
        = './', filename = None, show = True, title = ''):

        ''' Plots

        amp_functions = tuple of functions
        '''

        # Import external modules
        import numpy as np
        import pyfits as pf
        import matplotlib.pyplot as plt
        import matplotlib
        from mpl_toolkits.axes_grid1 import ImageGrid

        # Set up plot aesthetics
        plt.clf()
        plt.rcdefaults()
        fontScale = 10
        params = {#'backend': .pdf',
                  'axes.labelsize': fontScale,
                  'axes.titlesize': fontScale,
                  'text.fontsize': fontScale,
                  'legend.fontsize': fontScale*3/4,
                  'xtick.labelsize': fontScale,
                  'ytick.labelsize': fontScale,
                  'font.weight': 500,
                  'axes.labelweight': 500,
                  'text.usetex': False,
                  'figure.figsize': (10, 10),
                 }
        plt.rcParams.update(params)

        size = len(amp_functions)

        # Create figure
        fig = plt.figure()
        grid = ImageGrid(fig, (1,1,1),
                      nrows_ncols=(size,1),
                      ngrids = size,
                      direction='row',
                      axes_pad=1,
                      aspect=False,
                      share_all=True,
                      label_mode='All')

        colors = ['k','b','g','r','c']
        linestyles = ['-','--','-.','-','-']
        letters = ['a','b']

        for i, in range(1):
            ax = grid[i]

            ax.plot(wavelengths, magnitudes,
                    color = colors[j],
                    #label = 'G = %s' % G,
                    linestyle = '-',
                    )

            if limits is not None:
                ax.set_xlim(limits[0],limits[1])
                ax.set_ylim(limits[2],limits[3])

            # Adjust asthetics
            ax.set_xlabel(r'$\lambda$',)
            ax.set_ylabel(r'Magnitude')
            if size > 1:
                ax.annotate('(%s)' % letters[i],
                        xy = (0.9, 0.1),
                        xycoords='axes fraction',
                        textcoords='axes fraction')
            ax.grid(True)
            ax.legend(loc='lower left')

        if filename is not None:
            plt.savefig(savedir + filename,bbox_inches='tight', dpi=600)
        if show:
            fig.show()

def calc_flux2mag(flux):

    # Flux in erg s^-1 cm^-2 Hz^-1

    mag = -2.5 * log10(flux / 3.6308e-20)

    return mag

def problem_1():

    # import the data
    z004 = pf.open('bc03_padova1994_chab_z004_ssp.fit')
    z05 = pf.open('bc03_padova1994_chab_z05_ssp.fit')
    z02 = pf.open('bc03_padova1994_chab_z02_ssp.fit')

    z004_cols = z004[1].columns
    z05_cols = z05[1].columns
    z02_cols = z02[1].columns

    z004_data = z004[1].data
    z05_data = z05[1].data
    z02_data = z02[1].data

    names = z02_cols.names

    magnitudes = calc_flux2mag(z02_


def main():
    problem_1()

if __name__ == '__main__':
    main()




