#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import pyfits as pf

def plot_mass_vs_age(mass_list = None, age_list = None, list_names = None,
        limits = None, savedir = './', filename = None, show = True, title =
        '', redshifts = None):

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
        colormap = plt.cm.gist_ncar
        color_cycle = [colormap(i) for i in np.linspace(0, 0.9, len(mass_list))]
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
                  'figure.figsize': (6, 6),
                  'axes.color_cycle': color_cycle # colors of different plots
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

            for i, mass in enumerate(mass_list):
                ax.plot(age_list[i], mass,
                        label = '%s Gyr' % list_names[i],
                        marker = 's'
                        )
                if i < 2:
                    for j, z in enumerate(redshifts):
                        ax.annotate('z=%s' % z,
                                xy = (age_list[i][j], mass_list[i][j]),
                                textcoords = 'offset points',
                                xytext = (2,3),
                                size = fontScale * 0.75
                                )

            if limits is not None:
                ax.set_xlim(limits[0],limits[1])
                ax.set_ylim(limits[2],limits[3])

            # Adjust asthetics
            ax.set_xlabel(r'Age (Gyr)',)
            ax.set_ylabel(r'Mass ($M_\odot$)')
            ax.grid(True)
            ax.legend(loc='lower right')
            ax.set_title(title)

        if filename is not None:
            plt.savefig(savedir + filename,bbox_inches='tight', dpi=600)
        if show:
            fig.show()

def plot_color_index(g_i_list = None, i_list = None, list_names = None,
        limits = None, savedir = './', filename = None, show = True, title =
        '', redshifts = None):

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
        colormap = plt.cm.gist_ncar
        color_cycle = [colormap(i) for i in np.linspace(0, 0.9, len(i_list))]
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
                  'figure.figsize': (6, 6),
                  'axes.color_cycle': color_cycle # colors of different plots
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

            for i in range(len(i_list)):
                ax.plot(i_list[i], g_i_list[i],
                        label = '%s Gyr' % list_names[i],
                        marker = 's'
                        )
                if i < 2:
                    for j, z in enumerate(redshifts):
                        ax.annotate('z=%s' % z,
                                xy = (i_list[i][j], g_i_list[i][j]),
                                textcoords = 'offset points',
                                xytext = (2,3),
                                size = fontScale * 0.75
                                )

            if limits is not None:
                ax.set_xlim(limits[0],limits[1])
                ax.set_ylim(limits[2],limits[3])

            # Adjust asthetics
            ax.set_xlabel(r'$M_i$ (mag)',)
            ax.set_ylabel(r'$M_g - M_i$ (mag)')
            ax.grid(True)
            ax.legend(loc='bottom right')
            ax.set_title(title)

        if filename is not None:
            plt.savefig(savedir + filename,bbox_inches='tight', dpi=600)
        if show:
            fig.show()

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
        colormap = plt.cm.gist_ncar
        color_cycle = [colormap(i) for i in np.linspace(0, 0.9, len(mag_list))]
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
                  'figure.figsize': (6, 6),
                  'axes.color_cycle': color_cycle # colors of different plots
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

            # filters
            filters = ['U', 'B', 'V', 'R', 'I', 'J', 'H', 'K', 'NUV', 'FUV']
            filter_centers = [3630, 4450, 5510, 6580, 8060, 12200, 16300,
                    21900, 2274, 1542,]
            for j in range(len(filters)):
                ax.axvline(x = filter_centers[j],
                           ymin = 0, ymax = 1e10,
                           color = 'k')
                ax.annotate(filters[j],
                        xy = (filter_centers[j], -15),
                        textcoords = 'offset points',
                        xytext = (2,3)
                        )

            if limits is not None:
                ax.set_xlim(limits[0],limits[1])
                ax.set_ylim(limits[2],limits[3])

            # Adjust asthetics
            ax.set_xscale('log')
            ax.set_xlabel(r'$\lambda$ ($\AA$)',)
            ax.set_ylabel(r'$M_{AB}\ d\lambda$ (mag / $\AA$)')
            #ax.grid(True)
            ax.legend(loc='upper right')
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
        colormap = plt.cm.gist_ncar
        color_cycle = [colormap(i) for i in np.linspace(0, 0.9, len(ages))]
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
                  'figure.figsize': (6, 6),
                  'axes.color_cycle': color_cycle # colors of different plots
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
                    marker = 's',
                    color = 'k',
                    markersize = 3
                    #label = 'Age = %s Gyr' % ages[i],
                    )
            ax.axhline(y = 4.83 / 4.64,
                       xmin = -1, xmax = 100,
                       color = 'k')
            ax.annotate(r'$M_\odot / L_\odot$',
                    xy = (10, 4.83 / 4.64),
                    textcoords = 'offset points',
                    xytext = (2,3)
                    )

            if limits is not None:
                ax.set_xlim(limits[0],limits[1])
                ax.set_ylim(limits[2],limits[3])

            # Adjust asthetics
            #ax.set_xscale('log')
            ax.set_xlabel(r'Age (Gyr)',)
            ax.set_ylabel(r'$M / L (M_\odot / L_\odot$)')
            ax.grid(True)
            ax.legend(loc='upper right')
            ax.set_title(title)

        if filename is not None:
            plt.savefig(savedir + filename,bbox_inches='tight', dpi=600)
        if show:
            fig.show()

def plot_fluxes(wavelengths = None, flux_list = None, ages = None, metals =
        None, limits = None, savedir = './', filename = None, show = True,
        title = '', log_scale = (1,1), normalized = True, attenuations = None,
        age_unit = 'Gyr', balmer_line = False):

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
        colormap = plt.cm.gist_ncar
        color_cycle = [colormap(i) for i in np.linspace(0, 0.9, len(flux_list))]
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
                  'figure.figsize': (6, 6),
                  'axes.color_cycle': color_cycle # colors of different plots
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
            if balmer_line:
                lines = [6560,4861,4341,4102,3970,
                        3889, 3835,3646]
                for line in lines:
                    ax.axvline(x = line,
                           ymin = 0, ymax = 1e10,
                           color = 'k',
                           alpha = 0.5)

            for i, fluxes in enumerate(flux_list):
                if ages is not None and metals is None:
                    ax.plot(wavelengths, fluxes,
                            label = 'Age = %.1f %s' % (ages[i], age_unit)
                            )
                elif metals is not None and ages is None:
                    ax.plot(wavelengths, fluxes,
                            label = r'Z = %s ' % metals[i],
                            )
                elif ages is not None and metals is not None:
                    ax.plot(wavelengths, fluxes,
                            label = 'Age = %.1f %s, Z = %s ' % \
                                    (ages[i], age_unit, metals[i]),
                            )
                elif attenuations is not None:
                    ax.plot(wavelengths, fluxes,
                            label = r'$A_V = $ %s' % attenuations[i],
                            )

            if limits is not None:
                ax.set_xlim(limits[0],limits[1])
                ax.set_ylim(limits[2],limits[3])

            # Adjust asthetics
            if log_scale[0]:
                ax.set_xscale('log')
            if log_scale[1]:
                ax.set_yscale('log')
            ax.set_xlabel(r'$\lambda (\AA$)',)
            if normalized:
                ax.set_ylabel(r'$f_\lambda / f_{5500 \AA}$')
            else:
                ax.set_ylabel(r'$f_\lambda d\lambda$')

            ax.grid(True)
            ax.legend(loc='upper right')
            ax.set_title(title)

        if filename is not None:
            plt.savefig(savedir + filename,bbox_inches='tight', dpi=600)
        if show:
            fig.show()

def attenuate_flux(fluxes, wavelengths, Av):

    tau_v = Av / 1.086

    tau_d = tau_v * (wavelengths / 5500)**-0.7

    fluxes_attenuated = fluxes * np.exp(-tau_d)

    return fluxes_attenuated

def calc_flux2mag(flux):

    # Flux in erg s^-1 cm^-2 Hz^-1

    mag = -2.5 * np.log10(flux / 3.6308e-20)

    return mag

def calc_m2M(mag, distance):
    distance *= 3.24077929e-19 # cm to parsecs

    return mag - 5 * np.log10(distance / 10.)

def calc_mag2lum(Mag):

    lum = 10**(Mag / -2.5)

    return lum

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

    f = pf.open('manga_catalog_v0.1.fit')

    f_data = f[1].data

    ages = np.squeeze(f_data[columns[0]])
    wavelengths = np.squeeze(f_data[columns[1]])
    fluxes = np.squeeze(f_data[columns[2]])

    # assume a distance of 10 pc
    distance = 10 / 3.24077929e-19 # cm

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
                limits = (900, 25000, -15, -37),
                savedir = './', filename = 'p_1a.pdf',
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
                savedir = './', filename = 'p_1b.pdf',
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
                balmer_line = True,
                log_scale = (0,0),
                limits = (3500, 7000, 10**-4, 5),
                savedir = './', filename = 'p_1c.pdf',
                show = False, title = 'Problem 1c: SSP Normalized Spectra')

def problem_2(section = 'a'):

    #########
    # 2a
    #########

    data_names = ['z004', 'z05', 'z02']
    metallicities = [0.004, 0.05, 0.02]

    columns = ['AGE', 'WAVE', 'FLUX']

    # columns = ['Age', 'Wavelength', 'Flux']

    data = {}

    for i, name in enumerate(data_names):
        f = pf.open('bc03_padova1994_chab_%s_ssp.fit' % name)

        f_data = f[1].data

        data[name] = {'Z': metallicities[i],
                      'ages': np.squeeze(f_data[columns[0]]),
                      'wavelengths': np.squeeze(f_data[columns[1]]),
                      'fluxes': np.squeeze(f_data[columns[2]]),
                      }

    flux_norm_list = []
    wavelengths = data['z004']['wavelengths']
    ages = data['z004']['ages']
    wavelength_index = np.argmin(np.abs(wavelengths - 5500))

    age_index = np.argmin(np.abs(ages - 5*1e6))

    for i, metal in enumerate(data):
    	flux_norm = data[metal]['fluxes'][age_index]
    	flux_norm /= flux_norm[wavelength_index]

    	flux_norm_list.append(flux_norm)

    if section == 'a':
        plot_fluxes(flux_list = flux_norm_list,
                wavelengths = wavelengths,
                metals = metallicities,
                limits = (3500, 7000, 10**-1, 4),
                log_scale = (0,0),
                savedir = './', filename = 'p_2a.pdf',
                show = False, title = 'Problem 2a: SSP Normalized Spectra')

    #########
    # 2b
    #########

    ''' Find the age at which the high metallicity model best matches the 5 Gry
    old solar Z model.  Compare models normalized at 5500 and evaluate the fit
    over 3500 - 7000 . Do the same A A thing for the low metallicity model and
    plot then all three models together from 3500-7000 .  A On average, how
    large are the residuals? (1%, 10%, 20%?) '''

    age_index_solar = np.argmin(np.abs(ages - 5*1e9))
    wavelength_indices = np.where((wavelengths > 3500) & (wavelengths < 7000))
    wavelength_index = np.argmin(np.abs(wavelengths[wavelength_indices] - 5500))
    flux_z02 = (data['z02']['fluxes'][age_index_solar, wavelength_indices])[0]
    flux_norm_z02 = flux_z02 / flux_z02[wavelength_index]

    flux_norm_z02 = np.ma.array(flux_norm_z02,
            mask=(flux_norm_z02 != flux_norm_z02))

    residuals_array_z05 = np.zeros(len(ages))
    residuals_array_z004 = np.zeros(len(ages))

    for i in range(len(ages)):
    	flux_z05 = (data['z05']['fluxes'][i, wavelength_indices])[0]
    	flux_z004 = (data['z004']['fluxes'][i, wavelength_indices])[0]

        flux_norm_z05 = flux_z05 / flux_z05[wavelength_index]
        flux_norm_z004 = flux_z004 / flux_z004[wavelength_index]

        # mask NaNs
        flux_norm_z05 = np.ma.array(flux_norm_z05,
                mask=(flux_norm_z05 != flux_norm_z05))
        flux_norm_z004 = np.ma.array(flux_norm_z004,
                mask=(flux_norm_z004 != flux_norm_z004))

        residual = np.sum(((flux_norm_z05)**2 - (flux_norm_z02)**2)**0.5)
        residuals_array_z05[i] = residual

        residual = np.sum(((flux_norm_z004)**2 - (flux_norm_z02)**2)**0.5)
        residuals_array_z004[i] = residual

    residuals_array_z05 = np.ma.array(residuals_array_z05,
            mask=(residuals_array_z05 != residuals_array_z05))
    residuals_array_z004 = np.ma.array(residuals_array_z004,
            mask=(residuals_array_z004 != residuals_array_z004))

    if section == 'b':
        # print size of residuals
        resid_size_z05 = np.mean(residuals_array_z05) / np.sum(flux_norm_z02)
        resid_size_z004 = np.mean(residuals_array_z004) / np.sum(flux_norm_z02)
        print('Average residual size for z=0.05: %s %%' %  resid_size_z05)
        print('Average residual size for z=0.004: %s %%' %  resid_size_z004)

    age_index_z05 = np.where(residuals_array_z05 == residuals_array_z05.min())
    age_index_z004 = np.where(residuals_array_z004 == \
            residuals_array_z004.min())

    age_indices = [age_index_solar, age_index_z05[0], age_index_z004[0]]

    age_list = [5, ages[age_index_z05][0]/1e9, ages[age_index_z004][0]/1e9]

    metals = ['z02', 'z05', 'z004']
    metals_list = [0.02, 0.05, 0.004]
    flux_norm_list = []

    for i, metal in enumerate(metals):
        flux_norm = np.squeeze(data[metal]['fluxes'][age_indices[i],:])
    	flux_norm /= flux_norm[wavelength_index]
    	flux_norm_list.append(flux_norm)

    if section == 'b':
        plot_fluxes(flux_list = flux_norm_list,
                wavelengths = wavelengths,
                ages = age_list,
                metals = metals_list,
                log_scale = (0, 0),
                limits = (3500, 7000, 0, 2),
                savedir = './', filename = 'p_2b.pdf',
                show = False, title = 'Problem 2b: SSP Normalized Spectra')

    #########
    # 2c
    #########

    if section == 'c':
        plot_fluxes(flux_list = flux_norm_list,
                wavelengths = wavelengths,
                ages = age_list,
                metals = metals_list,
                log_scale = (0, 0),
                limits = (900, 30000, 10**-3, 1.5),
                savedir = './', filename = 'p_2c.pdf',
                show = False, title = 'Problem 2c: SSP Normalized Spectra')

def problem_3(section = 'a'):
    #########
    # 3a
    #########

    data_names = ['z004', 'z05', 'z02']
    metallicities = [0.004, 0.05, 0.02]

    columns = ['AGE', 'WAVE', 'FLUX']

    # columns = ['Age', 'Wavelength', 'Flux']

    data = {}

    for i, name in enumerate(data_names):
        f = pf.open('bc03_padova1994_chab_%s_ssp.fit' % name)

        f_data = f[1].data

        data[name] = {'Z': metallicities[i],
                      'ages': np.squeeze(f_data[columns[0]]),
                      'wavelengths': np.squeeze(f_data[columns[1]]),
                      'fluxes': np.squeeze(f_data[columns[2]]),
                      }

    flux_norm_list = []
    wavelengths = data['z02']['wavelengths']
    ages = data['z02']['ages']
    fluxes = data['z02']['fluxes']
    wavelength_index = np.argmin(np.abs(wavelengths - 5500))

    flux_list = []
    Av_list = [0.2, 0.5, 1]
    age_index = np.argmin(np.abs(ages - 100*1e6))

    for Av in Av_list:
        fluxes_attenuated = attenuate_flux(fluxes[age_index], wavelengths, Av)
        flux_list.append(fluxes_attenuated)

    if section == 'a':
        plot_fluxes(flux_list = flux_list,
                wavelengths = wavelengths,
                attenuations = Av_list,
                normalized = False,
                log_scale = (0, 1),
                limits = (900, 30000, 10**-6, 10**-2),
                savedir = './', filename = 'p_3a.pdf',
                show = False, title = 'Problem 3a: Attenuated SSP Spectra')


    #########
    # 3b
    #########

    Av = 0.5
    fluxes_solar = flux_list[1]
    age_index_solar = age_index

    wavelength_indices = np.where((wavelengths > 3500) & (wavelengths < 7000))
    flux_z02 = (data['z02']['fluxes'][age_index_solar, wavelength_indices])[0]
    wavelength_index = np.argmin(np.abs(wavelengths - 5500))

    flux_z02 = attenuate_flux(flux_z02[age_index_solar],
            wavelengths[wavelength_indices], Av)

    #flux_z02 /= flux_z02[wavelength_index]

    flux_z02 = np.ma.array(flux_z02,
            mask=(flux_z02 != flux_z02))

    residuals_array_z02 = np.zeros(len(ages))

    for i in range(len(ages)):
    	flux_z02_age = (data['z02']['fluxes'][i, wavelength_indices])[0]

        #flux_z02_age /= flux_z02_age[wavelength_index]

        # mask NaNs
        flux_z02_age = np.ma.array(flux_z02_age,
                mask=(flux_z02_age != flux_z02_age))

        residual = np.sum(((flux_z02_age)**2 - (flux_z02)**2)**0.5)
        residuals_array_z02[i] = residual

    residuals_array_z02 = np.ma.array(residuals_array_z02,
            mask=(residuals_array_z02 != residuals_array_z02))

    age_index_z02 = np.where(residuals_array_z02 == residuals_array_z02.min())

    age_indices = [age_index_solar, age_index_z02[0][0]]

    age_list = [100, ages[age_index_z02][0]/1e6]

    # redefine fluxes across entire wavelength range
    flux_z02 = (data['z02']['fluxes'][age_index_solar, :])
    flux_z02_age = (data['z02']['fluxes'][age_index_z02[0], :])[0]

    flux_z02 /= flux_z02[wavelength_index]
    flux_z02_age /= flux_z02_age[wavelength_index]

    flux_list = [flux_z02, flux_z02_age]

    if section == 'b':
        plot_fluxes(flux_list = flux_list,
                wavelengths = wavelengths,
                ages = age_list,
                age_unit = 'Myr',
                log_scale = (0, 1),
                limits = (900, 30000, 10**-2, 10**1),
                savedir = './', filename = 'p_3b.pdf',
                show = False,
                title = 'Problem 3b: Age vs. Attenuation SSP Spectra')

def problem_4(section = 'a'):
    #########
    # 4a
    #########

    bc_10 = np.loadtxt('bc03_tau10.txt')
    bc_1 = np.loadtxt('bc03_tau1.txt')
    m_10 = np.loadtxt('m05_tau10.txt')
    m_1 = np.loadtxt('m05_tau1.txt')

    data_list = [bc_10, bc_1, m_10, m_1]
    data_names = [r'BC03 $\tau = 10$',
                  r'BC03 $\tau = 1$',
                  r'M05 $\tau = 10$',
                  r'M05 $\tau = 1$',]
    columns = ['z', 'age', 'mass', 'sloan_g', 'sloan_i']

    mass_list = []
    age_list = []
    for data in data_list:
        age_list.append(data[:, 1])
        mass_list.append(data[:, 2])

    if section == 'a':
        plot_mass_vs_age(mass_list = mass_list,
                age_list = age_list,
                list_names = data_names,
                redshifts = data_list[0][:, 0],
                #limits = (900, 30000, 10**-6, 1),
                savedir = './',
                filename = 'p_4a.pdf',
                show = False,
                title = 'Problem 4a: SSP Mass with Age')
    #########
    # 4b
    #########

    g_i_list = []
    i_list = []
    for data in data_list:
        i_mag = data[:, 4]
        g_mag = data[:, 3]
        g_i_list.append(g_mag - i_mag)
        i_list.append(i_mag)

    if section == 'b':
        plot_color_index(g_i_list = g_i_list,
                i_list = i_list,
                list_names = data_names,
                redshifts = data_list[0][:, 0],
                limits = (6.5, 3.5, 1.4, -0.2),
                savedir = './',
                filename = 'p_4b.pdf',
                show = False,
                title = 'Problem 4b: Sloan g-i Mags with Time')
    #########
    # 4c
    #########

    # assume the ages of the two mergins galaxies are 9.86 Gyr old
    data_late = np.loadtxt('bc03_tau10.txt')
    data_early = np.loadtxt('bc03_tau1.txt')

    g_lum_early_10 = calc_mag2lum(data_early[-2, 3])
    i_lum_early_10 = calc_mag2lum(data_early[-2, 4])


    g_lum_early_5 = calc_mag2lum(data_early[5, 3])
    i_lum_early_5 = calc_mag2lum(data_early[5, 4])

    i_lum_late_10 = calc_mag2lum(data_late[-2, 4])
    g_lum_late_10 = calc_mag2lum(data_late[-2, 3])

    g_i_lum_early_10 = 1/g_lum_early_10 - 1/i_lum_early_10
    g_i_lum_late_10 = 1/g_lum_late_10 - 1/i_lum_late_10
    g_i_lum_early_5 = 1/g_lum_early_5 - 1/i_lum_early_5

    alpha = (1 / g_i_lum_early_5 - 1 / g_i_lum_late_10) \
            / (1 / g_i_lum_early_10 - 1 / g_i_lum_early_5)

    print('Mass_early / mass_late = %.2f' % alpha)

def main():
    #problem_1(section = 'a')
    #problem_1(section = 'b')
    #problem_1(section = 'c')
    #problem_1(section = 'd')

    #problem_2(section = 'a')
    #problem_2(section = 'b')
    #problem_2(section = 'c')

    #problem_3(section = 'a')
    #problem_3(section = 'b')

    problem_4(section = 'a')
    problem_4(section = 'b')
    problem_4(section = 'c')

if __name__ == '__main__':
    main()




