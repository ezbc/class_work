#!/usr/bin/python



def calc_L(radius, temperature):

    ''' Calculates luminosity.
    '''
    import numpy as np

    sigma = 5.67e-8 # W m^-2 K^-4

    return 4*np.pi * sigma * radius**2 * temperature**4

def plot_HR(temperatures, luminosities, normalize_L=False,
        filename='.png'):

    # Import external modules
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib

    # Create figure
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)
    if normalize_L:
    	luminosities = luminosities / 3.8e26
    	print luminosities
        ax.scatter(temperatures, luminosities)
    elif not normalize_L:
        ax.scatter(temperatures, luminosities)

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlim(temperatures.max()+0.2*temperatures.max(),
            temperatures.min() - 0.2*temperatures.min())

    # Adjust asthetics
    ax.set_xlabel('Temperature (K)',
              size = 'small',
              family='serif')
    if normalize_L:
        ax.set_ylabel(r'Luminosity (M$_{\odot}$)',
              size = 'small',
              family='serif')
        ax.set_ylim(1e-6,1e7)
    elif not normalize_L:
        ax.set_ylabel('Luminosity (W)',
              size = 'small',
              family='serif')


    ax.grid(True)

    fig.show()
    plt.savefig(filename, bbox_inches='tight')

def main():

    import numpy as np

    radii = np.array([0.2,  0.6,  1,    2,     8,     0.7,  0.8,  5])
    radii *= 6.98e8
    temps = np.array([1100, 1200, 1700, 10000, 11000, 4000, 6000, 10e4])

    data = np.array([[0.2, 1100],
                     [0.6, 1200],
                     [1,   1700],
                     [2,   1e4],
                     [8,   11e4],
                     [0.7, 4000],
                     [0.8, 6000],
                     [5,   1e4],
                     [5,   10**4.1],
                     [2,   9000],
                     [5,   1e4],
                     [5,   1e4],
                     [1,   10**3.8],
                     [5,   1e4],
                     [1,   10**3.8],
                     [5,   1e4],
                     [5,   1e4],
                     [6,   1e4],
                     [1,   10**3.6],
                     [4,   1e4],
                     [9,   10**4.9],
                     [.98, 10**3.75],
                     [0.45,4500],
                     [7,   10**4.2],
                     [1,   1e4],
                     [3,   1e4],
                     [5,   1e4],
                     [10,  1e4],
                     [10**-.5, 10**3.4],
                     [1, 10**3.8],
                     [9, 10**4.3],
                     [0.3, 10**3.3],
                     [1, 10**3.7],
                     [10**0.6, 10**4.2],
                     [10**.5, 10**4],
                     [10**0.5, 10**3.4],
                     [6, 10**4.1],
                     [9, 10**4.3],
                     [15, 10**4.6],
                     [0.5, 10**3.4],
                     [5, 10**4.2],
                     [15.84, 39810],
                     [6, 10**4.2],
                     [2.5, 10**3.9],
                     [10**0.5, 10**3.4],
                     [5, 10**4],
                     [5, 10**4.1],
                     [3.16, 10**5],
                     [1, 10**3.75],
                     [1, 10**3.76],
                     [0.5, 10**3.4],
                     [0.5, 10**3.45],
                     [10**0.6, 10**4.1],
                     [10**0.6, 10**4],
                     [10**0.1, 10**3.8],
                     [10**0.5, 10**3.5],
                     [10**0.5, 10**3.3],
                     [0.9, 10**3.7],
                     [1, 10**3.75],
                     [10**0.5, 10**4],
                     [10**0.75, 10**4.7],
                     [10**0.75, 10**4.25],
                     [10**0.5, 10**4],
                     [10**0.4, 10**4],
                     [0.98,10**3.73],
                     [10**1.3, 10**4.5],
                     [10**0.5, 10**4],
                     [1.01, 10**3.74],
                     [10**0.4, 10**4],
                     [1.02, 10**3.74]
                     ])


    radii = data[:,0]* 6.98e8
    temps = data[:,1]

    luminosities = calc_L(radii, temps)

    plot_HR(temps, luminosities, filename='hr_plot_watts.png')
    plot_HR(temps, luminosities, normalize_L=True, filename='hr_plot_Lsun.png')

if __name__ == '__main__':
	main()

