#!/usr/bin/python

'''

Writes plots for mach. learning presentation

'''
def gauss(x, sigma, x0, A):

    import numpy as np

    return A * np.exp(-(x - x0)**2 / (2.0 * sigma**2))

def gauss_1st_deriv(x, sigma, x0, A):

    ''' Returns the first derivative of a Gaussian. See gauss_4thderiv.nb
    mathematica code.

    '''

    import numpy as np

    expo = np.exp(-(x - x0)**2 / (2 * sigma**2))

    result = - A * expo * (x - x0) / sigma**2

    return result

def gauss_2nd_deriv(x, sigma, x0, A):

    ''' Returns the fourth derivative of a Gaussian. See gauss_4thderiv.nb
    mathematica code.

    '''

    import numpy as np

    expo = np.exp(-(x - x0)**2 / (2.0 * sigma**2))

    result = A * expo * (x - x0)**2.0 / sigma**4.0 - \
             A * expo * (x - x0)**0.0 / sigma**2.0

    return result

def gauss_3rd_deriv(x, sigma, x0, A):

    ''' Returns the fourth derivative of a Gaussian. See gauss_4thderiv.nb
    mathematica code.

    '''

    import numpy as np

    expo = np.exp(-(x - x0)**2 / (2.0 * sigma**2))

    result = -A * expo * (x - x0)**3.0 / sigma**6.0 + \
             3.0 * A * expo * (x - x0)**1.0 / sigma**4.0

    return result

def gauss_4th_deriv(x, sigma, x0, A):

    ''' Returns the fourth derivative of a Gaussian. See gauss_4thderiv.nb
    mathematica code.

    '''

    import numpy as np

    expo = np.exp(-(x - x0)**2 / (2.0 * sigma**2))

    result = A * expo * (x - x0)**4.0 / sigma**8.0 - \
             6.0 * A * expo * (x - x0)**2.0 / sigma**6.0 + \
             3.0 * A * expo / sigma**4.0

    return result

def gauss(x, sigma, x0, A):

    import numpy as np

    return A * np.exp(-(x - x0)**2 / (2.0 * sigma**2))

def plot_temps():

    import numpy as np
    import matplotlib.pyplot as plt

    x = np.arange(0, 100, 0.1)

    gauss1 = gauss(x, 5, 50, 1)
    gauss2 = gauss(x, 10, 50, 1)
    gauss3 = gauss(x, 15, 50, 1)


    gauss1 /= np.sum(gauss1)
    gauss2 /= np.sum(gauss2)
    gauss3 /= np.sum(gauss3)

    scale = np.max(np.array((gauss1, gauss2, gauss3)))

    gauss1 /= scale
    gauss2 /= scale
    gauss3 /= scale

    # Set up plot aesthetics
    plt.clf()
    plt.close()
    plt.rcdefaults()
    colormap = plt.cm.gist_ncar
    font_scale = 20
    params = {#'backend': .pdf',
              'axes.labelsize': font_scale,
              'axes.titlesize': font_scale,
              'text.fontsize': font_scale,
              'legend.fontsize': font_scale * 4.0 / 4.0,
              'xtick.labelsize': font_scale,
              'ytick.labelsize': font_scale,
              'font.weight': 500,
              'axes.labelweight': 500,
              'text.usetex': False,
              #'figure.figsize': (8, 8 * y_scaling),
              #'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)

    fig, ax = plt.subplots(1, 1, figsize=(7, 7))

    ax.plot(x, gauss1, label='Cold', linewidth=3)
    ax.plot(x, gauss2, label='Warm', linewidth=3)
    ax.plot(x, gauss3, label='Hot', linewidth=3)

    ax.set_xlabel('Velocity (km/s)')
    ax.set_ylabel('Intensity')
    ax.legend()

    plt.savefig('gauss_temps.png')

def plot_gauss():

    import numpy as np
    import matplotlib.pyplot as plt

    x = np.arange(0, 100, 0.1)

    gauss1 = gauss(x, 15, 50, 1)

    gauss1 /= np.sum(gauss1)

    scale = np.max(np.array((gauss1,)))

    gauss1 /= scale

    # Set up plot aesthetics
    plt.clf()
    plt.close()
    plt.rcdefaults()
    colormap = plt.cm.gist_ncar
    font_scale = 20
    params = {#'backend': .pdf',
              'axes.labelsize': font_scale,
              'axes.titlesize': font_scale,
              'text.fontsize': font_scale,
              'legend.fontsize': font_scale * 4.0 / 4.0,
              'xtick.labelsize': font_scale,
              'ytick.labelsize': font_scale,
              'font.weight': 500,
              'axes.labelweight': 500,
              'text.usetex': False,
              #'figure.figsize': (8, 8 * y_scaling),
              #'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)

    fig, ax = plt.subplots(1, 1, figsize=(7, 7))

    ax.plot(x, gauss1, linewidth=3)

    ax.set_xlabel('Velocity (km/s)')
    ax.set_ylabel('Intensity')

    plt.savefig('gauss.png')

def plot_data():

    import numpy as np
    import matplotlib.pyplot as plt

    x = np.arange(-100, 100, 0.1)

    gauss1 = gauss(x, 1, 20, 10)
    gauss2 = gauss(x, 10, 5, 15)
    gauss3 = gauss(x, 15, 50, 4)
    gauss4 = gauss(x, 7, 60, 30)
    gauss5 = gauss(x, 50, 35, 4)

    noise = np.random.normal(0, 0.5, len(x))

    gauss_list = [gauss1, gauss2, gauss3, gauss4, gauss5, ]

    gauss_tot = np.zeros(len(x))
    for comp in gauss_list:
        gauss_tot += comp
    gauss_tot += noise

    scale = np.max(gauss_tot)

    gauss_tot /= scale

    for i, comp in enumerate(gauss_list):
        gauss_list[i] = comp / scale

    # Set up plot aesthetics
    plt.clf()
    plt.close()
    plt.rcdefaults()
    colormap = plt.cm.gist_ncar
    font_scale = 15
    params = {#'backend': .pdf',
              'axes.labelsize': font_scale,
              'axes.titlesize': font_scale,
              'text.fontsize': font_scale,
              'legend.fontsize': font_scale * 4.0 / 4.0,
              'xtick.labelsize': font_scale,
              'ytick.labelsize': font_scale,
              'font.weight': 500,
              'axes.labelweight': 500,
              'text.usetex': False,
              #'figure.figsize': (8, 8 * y_scaling),
              #'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)

    fig, ax = plt.subplots(1, 1, figsize=(7, 7))

    ax.plot(x, gauss_tot, linewidth=1, color='k')

    ax.set_xlabel('Velocity (km/s)')
    ax.set_ylabel('Intensity')

    fig.savefig('data.png')


    fig, ax = plt.subplots(1, 1, figsize=(7, 7))

    ax.plot(x, gauss_tot, linewidth=1, color='k')

    ax.set_xlabel('Velocity (km/s)')
    ax.set_ylabel('Intensity')
    for comp in gauss_list:
        ax.plot(x, comp, linewidth=1,)

    fig.savefig('data_comps.png')

def plot_derivs():

    import numpy as np
    import matplotlib.pyplot as plt

    x = np.arange(0, 100, 0.1)

    params = (x, 15, 50, 1)

    gauss1 = gauss(*params)
    gauss1d = gauss_1st_deriv(*params)
    gauss2d = gauss_2nd_deriv(*params)
    gauss3d = gauss_3rd_deriv(*params)
    gauss4d = gauss_4th_deriv(*params)

    deriv_list = [gauss1d, gauss2d, gauss3d, gauss4d,]

    for i, comp in enumerate(deriv_list):
        deriv_list[i] = comp / np.max(comp)

    # Set up plot aesthetics
    plt.clf()
    plt.close()
    plt.rcdefaults()
    colormap = plt.cm.gist_ncar
    font_scale = 15
    params = {#'backend': .pdf',
              'axes.labelsize': font_scale,
              'axes.titlesize': font_scale,
              'text.fontsize': font_scale,
              'legend.fontsize': font_scale * 4.0 / 4.0,
              'xtick.labelsize': font_scale,
              'ytick.labelsize': font_scale,
              'font.weight': 500,
              'axes.labelweight': 500,
              'text.usetex': False,
              #'figure.figsize': (8, 8 * y_scaling),
              #'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)

    fig, ax = plt.subplots(1, 1, figsize=(7, 7))
    ax.plot(x, gauss1, linewidth=3, color='k')
    ax.set_ylim(-2.5, 1.1)
    ax.set_xlabel('Velocity (km/s)')
    ax.set_ylabel('Intensity')

    plt.savefig('gauss_deriv0.png')

    for i, deriv in enumerate(deriv_list):

        ax.plot(x, deriv, linewidth=1)

        plt.savefig('gauss_deriv' + str(i+1) + '.png')

def main():

    plot_temps()
    plot_gauss()
    plot_data()
    plot_derivs()


if __name__ == '__main__':

    main()


