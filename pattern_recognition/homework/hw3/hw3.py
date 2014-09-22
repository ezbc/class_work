#!/usr/bin/python

import numpy as np

def problem_1():

    A = np.array([[1, 1],
                  [1, -1],
                  [1, 1]])

    b = [1, 1, 0]

    x = np.linalg.lstsq(A, b)

    print 'x = ', x[0]

def problem_2():

    import matplotlib.pyplot as plt

    m = 1000

    a = np.random.uniform(0, 1, size=m)
    A = np.matrix([a, np.ones((m))])
    A = A.T

    b = a**2 + np.random.normal(0, 1, size=m)
    B = np.matrix(b)
    B = B.T

    # Fit b
    b_fit = np.linalg.inv(A.T * A) * A.T * B

    # Plot!
    # Set up plot aesthetics
    plt.clf()
    plt.rcdefaults()
    colormap = plt.cm.gist_ncar
    #color_cycle = [colormap(i) for i in np.linspace(0, 0.9, len(flux_list))]
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
              #'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)
    fig = plt.figure(figsize=(3,2))

    ax = fig.add_subplot(111)

    ax.plot(a, b, linestyle='', marker='^', alpha=0.4)
    ax.plot(a, b_fit, color='r')
    ax.set_xlabel(r'$a_i$')
    ax.set_ylabel(r'$b_i$')

    plt.savefig('problem2c.png', bbox_inches='tight')

def main():
    problem_1()
    problem_2()

if __name__ == '__main__':
	main()
