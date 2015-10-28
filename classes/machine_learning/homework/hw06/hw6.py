#!/usr/bin/python



def problem1():

    from scipy.interpolate import UnivariateSpline as fit_spline
    import numpy as np
    import matplotlib.pyplot as plt

    data = np.genfromtxt('hw6.dat', skip_header=1)

    x, y = data.T[0], data.T[1]

    k = 3

    lams = np.exp((-2, -1, 0, 1, 2))
    x_fits = []
    y_fits = []

    for i, lam in enumerate(lams):
        spline = fit_spline(x, y, s=lam, k=1)

        x_fit = np.linspace(x.min(), x.max(), 1000)
        y_fit = spline(x_fit)

        x_fits.append(x_fit)
        y_fits.append(y_fit)

    # Fit spline with GCV
    spline = fit_spline(x, y, k=k)
    lam_gcv = spline._data[6]

    x_fit_gcv = np.linspace(x.min(), x.max(), 1000)
    y_fit_gcv = spline(x_fit_gcv)

    # Plot!
    # Set up plot aesthetics
    plt.clf()
    plt.rcdefaults()
    colormap = plt.cm.gist_ncar
    color_cycle = [colormap(i) for i in np.linspace(0, 0.9, len(lams))]
    font_scale = 10
    params = {#'backend': .pdf',
              'axes.labelsize': font_scale,
              'axes.titlesize': font_scale,
              'text.fontsize': font_scale,
              'legend.fontsize': font_scale * 3 / 4.0,
              'xtick.labelsize': font_scale,
              'ytick.labelsize': font_scale,
              'font.weight': 500,
              'axes.labelweight': 500,
              'text.usetex': False,
              #'figure.figsize': (8, 8 * y_scaling),
              'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)
    fig = plt.figure(figsize=(6,4))

    ax = fig.add_subplot(111)

    ax.plot(x, y, linestyle='', marker='^', alpha=0.4)
    for i in xrange(len(x_fits)):
        x_fit = x_fits[i]
        y_fit = y_fits[i]
        lam = lams[i]
        ax.plot(x_fit, y_fit, label=r'$\lambda$ = {0:.1f}'.format(lam))
    ax.plot(x_fit_gcv, y_fit_gcv, color='r',
            label=r'GCV, $\lambda$ = {0:.0f}'.format(lam_gcv))
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$y$')
    ax.legend(loc='best')

    plt.savefig('problem1_fig.png', bbox_inches='tight')

def main():
    problem1()

if __name__ == '__main__':
    main()
