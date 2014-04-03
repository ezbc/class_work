#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import pyfits as pf

def plot_mags(wavelengths = None, mag_list = None, ages = None, limits = None,
        savedir = './', filename = None, show = True, title = '',):

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
        fontScale = 12
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
                  'figure.figsize': (8, 8),
                 }
        plt.rcParams.update(params)

        # Create figure
        fig = plt.figure()
        grid = ImageGrid(fig, (1,1,1),
                      nrows_ncols=(1,1),
                      ngrids = 1,
                      direction='row',
                      axes_pad=1,
                      aspect=False,
                      share_all=True,
                      label_mode='All')

        colors = ['k','b','g','r','c']
        linestyles = ['-','--','-.','-','-']
        letters = ['a','b']

        for i in range(1):
            ax = grid[i]

            for i, mags in enumerate(mag_list):
                ax.plot(wavelengths, mags,
                        label = 'Age = %s Gyr' % ages[i],
                        )

            if limits is not None:
                ax.set_xlim(limits[0],limits[1])
                ax.set_ylim(limits[2],limits[3])

            # Adjust asthetics
            ax.set_xlabel(r'$\lambda$ (Angstrom)',)
            ax.set_ylabel(r'$M_{AB}$ d$\lambda$')
            ax.grid(True)
            ax.legend(loc='lower right')
            ax.set_title(title)

        if filename is not None:
            plt.savefig(savedir + filename,bbox_inches='tight', dpi=600)
        if show:
            fig.show()

def plot_mass2light(ages = None, m2l = None, limits =
        None, savedir = './', filename = None, show = True, title = '',):

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
        fontScale = 12
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
                  'figure.figsize': (8, 8),
                 }
        plt.rcParams.update(params)

        # Create figure
        fig = plt.figure()
        grid = ImageGrid(fig, (1,1,1),
                      nrows_ncols=(1,1),
                      ngrids = 1,
                      direction='row',
                      axes_pad=1,
                      aspect=False,
                      share_all=True,
                      label_mode='All')

        colors = ['k','b','g','r','c']
        linestyles = ['-','--','-.','-','-']
        letters = ['a','b']

        for i in range(1):
            ax = grid[i]

            ax.plot(ages, m2l,
                    #label = 'Age = %s Gyr' % ages[i],
                    )

            if limits is not None:
                ax.set_xlim(limits[0],limits[1])
                ax.set_ylim(limits[2],limits[3])

            # Adjust asthetics
            #ax.set_xscale('log')
            ax.set_xlabel(r'Age (Gyr)',)
            ax.set_ylabel(r'$M / L (M_\odot / L_\odot$)')
            ax.grid(True)
            ax.legend(loc='lower right')
            ax.set_title(title)

        if filename is not None:
            plt.savefig(savedir + filename,bbox_inches='tight', dpi=600)
        if show:
            fig.show()

def plot_fluxes(wavelengths = None, flux_list = None, ages = None, limits =
        None, savedir = './', filename = None, show = True, title = '',):

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
        fontScale = 12
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
                  'figure.figsize': (8, 8),
                 }
        plt.rcParams.update(params)

        # Create figure
        fig = plt.figure()
        grid = ImageGrid(fig, (1,1,1),
                      nrows_ncols=(1,1),
                      ngrids = 1,
                      direction='row',
                      axes_pad=1,
                      aspect=False,
                      share_all=True,
                      label_mode='All')

        colors = ['k','b','g','r','c']
        linestyles = ['-','--','-.','-','-']
        letters = ['a','b']

        for i in range(1):
            ax = grid[i]

            for i, fluxes in enumerate(flux_list):
                ax.plot(wavelengths, fluxes,
                        label = 'Age = %s Gyr' % ages[i],
                        )

            if limits is not None:
                ax.set_xlim(limits[0],limits[1])
                ax.set_ylim(limits[2],limits[3])

            # Adjust asthetics
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_xlabel(r'$\lambda (\AA$)',)
            ax.set_ylabel(r'$f_\lambda / f_{5500 \AA}$')
            ax.grid(True)
            ax.legend(loc='lower right')
            ax.set_title(title)

        if filename is not None:
            plt.savefig(savedir + filename,bbox_inches='tight', dpi=600)
        if show:
            fig.show()

def calc_flux2mag(flux):

    # Flux in erg s^-1 cm^-2 Hz^-1

    mag = -2.5 * np.log10(flux / 3.6308e-20)

    return mag

def calc_m2M(mag, distance):
    distance *= 3.24077929e-19 # cm to parsecs

    return mag - 5 * np.log10(distance / 10.)

def calc_z2dist(z):

    # c*z = H_0 * D where H_0 = 70 km / s / Mpc = 0.07 m / s / m

    c = 2.99e8 # m/s
    H_0 = 0.07 # m / s / m

    dist = c*z / H_0

    return dist

def problem_1(section = 'a'):
    #########
    # 1a
    #########

    columns = ['AGE', 'WAVE', 'FLUX']

    f = pf.open('bc03_padova1994_chab_%s_ssp.fit' % 'z02')

    f_data = f[1].data

    ages = np.squeeze(f_data[columns[0]])
    wavelengths = np.squeeze(f_data[columns[1]])
    fluxes = np.squeeze(f_data[columns[2]])

    # assume a distance of 1 Mpc
    distance = 1e6 / 3.24077929e-19 # cm

    fluxes_cgs = fluxes / (4 * np.pi * distance**2) * 3.839e33

    mag_ab = -2.5 * np.log10(fluxes_cgs / 3.6308e-20)

    Mags_ab = calc_m2M(mag_ab, distance)

    age_list = (0.001, 0.01, 0.1, 0.3, 0.5, 1, 5, 10)

    Mag_list = []
    for age in age_list:
        age_indices = np.argmin(np.abs(ages - age*1e9))
        Mags = Mags_ab[age_indices]
        Mag_list.append(Mags)

    if section == 'a':
        plot_mags(wavelengths = wavelengths, mag_list = Mag_list,
                ages = age_list,
                limits = (900, 25000, -15, -40),
                savedir = './', filename = 'p_1a.png',
                show = False, title = 'Problem 1a: SSPs')

    #########
    # 1b
    #########

    wavelength_index = np.argmin(np.abs(wavelengths - 5500))

    m2l = 1 / (fluxes[:, wavelength_index] * 5500.)

    if section == 'b':
        plot_mass2light(m2l = m2l,
                ages = ages / 1e9,
                #limits = (10**-2, 10**0, 0, 10**23),
                savedir = './', filename = 'p_1b.png',
                show = False, title = 'Problem 1b: SSP Mass to Light Ratios')

    #########
    # 1c
    #########

    flux_norm_list = []

    for age in age_list:
        age_indices = np.argmin(np.abs(ages - age*1e9))
        flux_norm = fluxes[age_indices] / fluxes[age_indices, wavelength_index]
        flux_norm_list.append(flux_norm)

    if section == 'c':
        plot_fluxes(flux_list = flux_norm_list,
                wavelengths = wavelengths,
                ages = age_list,
                limits = (900, 25000, 10**-4, 500),
                savedir = './', filename = 'p_1c.png',
                show = False, title = 'Problem 1c: SSP Normalized Spectra')

def problem_3():

    data_names = ['z004', 'z05', 'z02']
    redshifts = [0.004, 0.05, 0.02]

    columns = ['AGE', 'WAVE', 'FLUX']

    # columns = ['Age', 'Wavelength', 'Flux']

    data = {}

    for i, name in enumerate(data_names):
        f = pf.open('bc03_padova1994_chab_%s_ssp.fit' % name)

        f_data = f[1].data

        for column in columns:
            data[name] = {'redshift': redshifts[i],
                          'age': np.squeeze(f_data[column]),
                          'wavelength': np.squeeze(f_data[column]),
                          'flux': np.squeeze(f_data[column]),
                          }
        w = data[name]['wavelength']
        print w.max(), w.min()

    for redshift in data:
    	datum = data[redshift]
        datum['distance'] = calc_z2dist(datum['redshift'])

        datum['mag_ab'] = calc_flux2mag(datum['flux'])

        #datum['Mag_ab'] = -2.5 * np.log10(datum['flux'])
        #datum['Mag_ab'] = calc_lum2mag(datum['flux'])
        datum['Mag_ab'] = calc_m2M(datum['mag_ab'], datum['distance'])

    shape = data['z02']['Mag_ab'].shape
    Mags = np.zeros((shape[0], shape[1], 3))
    wavelengths = np.copy(Mags)

    for i, key in enumerate(data):
        Mags[:,:,i] = data[key]['Mag_ab']
        wavelengths[:,:,i] = data[key]['wavelength']

    Mags = np.ravel(Mags)
    wavelengths = np.ravel(wavelengths)

    plot_mags(wavelengths = wavelengths, magnitudes = Mags, limits =
            (900, 25000, 10**-3, 10**2), savedir = './', filename = 'p_1a.png',
            show = True, title = 'Problem 1a: SSPs')

def main():
    #problem_1(section = 'a')
    #problem_1(section = 'b')
    problem_1(section = 'c')
    #problem_1(section = 'd')

if __name__ == '__main__':
    main()




