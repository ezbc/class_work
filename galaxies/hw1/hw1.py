#!/usr/bin/python

def problem_1b(show=True):
    import matplotlib.pyplot as plt
    import numpy as np

    phi_star = 1.64e-2 # h^3 Mpc^-2
    h = 0.7
    m_star = -19.67 + 5*np.log10(h)
    alpha = -1.21
    m_array = np.arange(-24 + 5*np.log10(h), -10 + 5*np.log10(h),0.01)
    #print m_array

    prob_dens = schechter(m_array, alpha, m_star, phi_star)

    # Create figure
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)

    ax.plot((m_array - 5*np.log10(h)).T, (np.log10(prob_dens)).T )

    ax.set_ylim([-6,-1])
    ax.set_xlim([-13.1,-22.9])

    # Adjust asthetics
    ax.set_xlabel(r'M$_{\rm bj}$ - 5log(h)',
              size = 'small',
              family='serif')
    ax.set_ylabel(r'log($\phi$ / [h$^3$Mpc$^{-3}$mag$^{-1}$])',
              size = 'small',
              family='serif')
    ax.set_title('1b: Luminosity function of galaxies shown in 2.25')
    ax.grid(True)

    if True:
        plt.savefig('hw1.1b.png',bbox_inches='tight',dpi=500)
    if show:
        fig.show()

def schechter(m_array, alpha, m_star, phi_star):
    import numpy as np

    prob_dens = np.log(10) / 2.5 * phi_star * \
            (10**(-0.4*(m_array - m_star)))**(alpha+1) *\
            np.exp(-10**(-0.4*(m_array - m_star)))

    return prob_dens

def problem_1c(show=True, N=10**3):
    import matplotlib.pyplot as plt
    import numpy as np

    # define schechter function parameters
    phi_star = 1.64e-2 # h^3 Mpc^-2
    h = 0.7
    m_star = -19.67 + 5*np.log10(h)
    alpha = -1.21
    M_binsize = 0.1
    m_array = np.arange(-22.5 + 5*np.log10(h), -13.1 + 5*np.log10(h),
            M_binsize)

    # create distribution of M
    prob_dens = schechter(m_array, alpha, m_star, phi_star)

    # integrate
    prob = (prob_dens * M_binsize).sum()
    volume = N / prob * M_binsize

    # get number of galaxies with Magnitude M
    N_galaxies_dM = prob_dens * volume

    # galaxies = N_galaxies / (volume * (m_array[1] - m_array[0]))

    distances = np.random.random(N) * .98e9 + .020e6 # Mpc

    # Move each individual galaxy to given distance,
    # calculate apparent magnitude
    count = 0
    apparent_mags = np.zeros(distances.shape)
    for i in range(len(m_array)):
        for j in range(int(N_galaxies_dM[i])):
            apparent_mags[count + j] = m_array[i] + \
                    5*np.log10(distances[count + j]/10.)
        count += N_galaxies_dM[i]

    plot_problem_1c(apparent_mags, show=show)

def plot_problem_1c(apparent_mags, show=True):
    import matplotlib.pyplot as plt
    import numpy as np

    # Create figure
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)

    ax.hist(apparent_mags,bins=100, histtype='step', color='k')

    #ax.set_ylim([-6,-1])
    ax.set_xlim([27,10])

    # Adjust asthetics
    ax.set_xlabel(r'm$_{bj}$ [magnitudes]',
              size = 'small',
              family='serif')
    ax.set_ylabel(r'N',
              size = 'small',
              family='serif')
    ax.set_title('1c: Number density of galaxy apparent magnitudes.')
    ax.grid(True)

    if True:
        plt.savefig('hw1.1c.png',bbox_inches='tight')
    if show:
        fig.show()

def problem_1d(show=True, N=10**5):
    import matplotlib.pyplot as plt
    import numpy as np

    # define schechter function parameters
    phi_star = 1.64e-2 # h^3 Mpc^-2
    h = 0.7
    m_star = -19.67 + 5*np.log10(h)
    alpha = -1.21
    M_binsize = 0.1
    m_array = np.arange(-22.5 + 5*np.log10(h), -13.1 + 5*np.log10(h),
            M_binsize)

    # create distribution of M
    prob_dens = schechter(m_array, alpha, m_star, phi_star)

    # integrate
    prob = (prob_dens * M_binsize).sum()
    volume = N / prob * M_binsize

    # get number of galaxies with Magnitude M
    N_galaxies_dM = prob_dens * volume
    # galaxies = N_galaxies / (volume * (m_array[1] - m_array[0]))

    distances = np.random.random(N) * .98e9 + .020e6 # Mpc

    # Move each individual galaxy to given distance,
    # calculate apparent magnitude
    #apparent_mags = np.zeros(N)
    #for i, distance in enumerate(distances):
    #    for j in range(len(N_galaxies_dM)):
    #        apparent_mags[i] = m_array[j] + 5*np.log10(distance/10.)

    # Move each individual galaxy to given distance,
    # calculate apparent magnitude
    count = 0
    apparent_mags = np.zeros(distances.shape)
    for i in range(len(m_array)):
        for j in range(int(N_galaxies_dM[i])):
            apparent_mags[count + j] = m_array[i] + \
                    5*np.log10(distances[count + j]/10.)
        count += N_galaxies_dM[i]

    count = 0
    N_galaxies_cut = np.zeros(m_array.shape)
    for i in range(len(m_array)):
        for j in range(int(N_galaxies_dM[i])):
            if apparent_mags[count + j] < 20:
                N_galaxies_cut[i] += 1
        count += N_galaxies_dM[i]

    plot_problem_1d(m_array=m_array, N_galaxies_dM=N_galaxies_dM,
        N_galaxies_cut=N_galaxies_cut, apparent_mags=apparent_mags, show=show)

def plot_problem_1d(m_array=None, N_galaxies_dM=None,
        N_galaxies_cut=None, apparent_mags=None, show=True):
    import matplotlib.pyplot as plt
    import numpy as np

    # Create figure
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)

    ax.plot(m_array,N_galaxies_dM,
            color='g', label='Untrimmed',
            drawstyle='steps')

    ax.plot(m_array,N_galaxies_cut,
            color='r', label='Magnitude-limited',
            drawstyle='steps-pre')

    ax.set_xlim([-13.1,-22.9])
    ax.set_ylim([1,10**4])
    ax.set_yscale('log')

    # Adjust asthetics
    ax.set_xlabel(r'M [magnitude]',
              size = 'small',
              family='serif')
    ax.set_ylabel(r'N',
              size = 'small',
              family='serif')
    ax.set_title('1d: Histogram of Absolute Magnitudes')
    ax.grid(True)
    ax.legend(loc='lower left')

    if True:
        plt.savefig('hw1.1d.png',bbox_inches='tight',dpi=500)
    if show:
        fig.show()

def problems_1fgh(show=True, N=10**5, distance_bin=1):

    '''
    Parameters
    ----------
    distance_bin : float
        In units of Mpc

    '''

    import matplotlib.pyplot as plt
    import numpy as np

    # define schechter function parameters
    phi_star = 1.64e-2 # h^3 Mpc^-2
    h = 0.7
    m_star = -19.67 + 5*np.log10(h)
    alpha = -1.21
    M_binsize = 0.1
    m_array = np.arange(-22.5 + 5*np.log10(h), -13.1 + 5*np.log10(h),
            M_binsize)

    # create distribution of M
    prob_dens = schechter(m_array, alpha, m_star, phi_star)

    # integrate
    prob = (prob_dens * (m_array[1] - m_array[0])).sum()
    volume = N / prob * M_binsize

    # get number of galaxies with Magnitude M
    N_galaxies_dM_untrimmed = prob_dens * volume

    # Calulate random distances between 20 and 1000 Mpc
    distances = np.random.random(N) * .98e9 + .020e9 # pc
    indices = distances.argsort()

    # There must be at least 1 galaxy in the last slice
    N_1000 = 0
    while N_1000 == 0:
        # find number of galaxies at farthest distance
        N_1000 = distances[distances > 1e9 - distance_bin * 1e6].size
        distance_bin += 0.1

    # find volume of farthest distance
    volume_1000 = N_1000 / prob * M_binsize

    N_volume_bins = int((1000 - 20) / distance_bin)

    N_galaxies_dM_array = np.zeros(shape=(prob_dens.shape[0], N_volume_bins))

    volume_farther = volume_1000
    distance_farther = 1000 - distance_bin
    N_galaxies_dM_array[:,0] = N_1000

    # Define the number of galaxies in each volume bin
    for i in range(1,N_volume_bins):
        distance_closer = distance_farther - distance_bin
        volume_closer = volume_farther * \
                distance_closer**2 / distance_farther**2

        # get number of galaxies between Magnitude M and M+dM
        N_galaxies_dM_array[:,i] = prob_dens * volume_closer

        # perform again on next volume slice
        volume_farther = volume_closer
        distance_farther = distance_closer

    # Calculate total numbe rof galaxies
    N_galaxies = N_galaxies_dM_array.sum()

    print('Total number of galaxies in volume-limited sample = ' + \
            str(int(N_galaxies)))

    # Calculate the number of galaxies between absolute mag M and M + dM
    N_galaxies_dM = N_galaxies_dM_array.sum(axis=1)

    # Calulate random distances between 20 and 1000 Mpc now for the number of
    # galaxies in the volume limited sample
    distances = np.random.random(N_galaxies) \
            * .98e9 + .020e9 # pc

    # calculate apparent magnitude
    count = 0
    apparent_mags = np.zeros(distances.shape)
    for i in range(len(m_array)):
        for j in range(int(N_galaxies_dM[i])):
            apparent_mags[count + j] = m_array[i] + \
                    5*np.log10(distances[count + j]/10.)
        count += N_galaxies_dM[i]

    # Find galaxies with m < 20
    count = 0
    N_galaxies_cut = np.zeros(m_array.shape)
    distances_cut = []
    m_array_distances = []
    for i in range(len(m_array)):
        for j in range(int(N_galaxies_dM[i])):
            if apparent_mags[count + j] < 20:
                N_galaxies_cut[i] += 1
                distances_cut.append(distances[count + j])
                m_array_distances.append(m_array[i])
        count += N_galaxies_dM[i]

    # problem 1f
    ############
    plot_problem_1f(m_array=m_array, N_galaxies_cut=N_galaxies_cut,
            N_galaxies_dM=N_galaxies_dM,
            N_galaxies_dM_untrimmed=N_galaxies_dM_untrimmed, show=show)


    # problem 1g
    ############
    apparent_mags_binsize = 50
    apparent_mags_cut = apparent_mags[apparent_mags < 20]
    bins = len(apparent_mags_cut) / apparent_mags_binsize

    N_apparent_mags, apparent_mags_bin = \
            np.histogram(apparent_mags_cut, bins=bins)

    # calculate fraction near limiting magnitude
    limiting_indices = np.where(apparent_mags_bin[:-1] > 19.5)[0]
    sample_total_N = N_apparent_mags.sum()
    limiting_fraction = N_apparent_mags[limiting_indices[-1]] / \
            float(sample_total_N)
    print('Fraction of sample within 0.5 magnitudes of the limiting' + \
            ' magnitude = %.2f' % limiting_fraction)

    plot_problem_1g(apparent_mags_bin=apparent_mags_bin,
            N_apparent_mags=N_apparent_mags, show=show)


    # problem 1h
    ############
    plot_problem_1h(m_array_distances=np.asarray(m_array_distances),
            distances_cut=np.asarray(distances_cut), show=show)

def plot_problem_1f(m_array=None, N_galaxies_cut=None, N_galaxies_dM=None,
        N_galaxies_dM_untrimmed=None, show=True):

    import matplotlib.pyplot as plt
    import numpy as np
    # Create figure
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)

    ax.plot(m_array,N_galaxies_dM_untrimmed,
            color='b', label='Untrimmed',  linestyle='-',
            drawstyle='steps')
    ax.plot(m_array,N_galaxies_dM,
            color='r', label='Volume-limited',  linestyle='-',
            drawstyle='steps')
    ax.plot(m_array,N_galaxies_cut,
            color='g', label='Magnitude- and Volume-limited',  linestyle='-',
            drawstyle='steps')

    ax.set_xlim([-13.1,-22.9])
    ax.set_ylim([1,10**4])
    ax.set_yscale('log')

    # Adjust asthetics
    ax.set_xlabel(r'M [magnitude]',
              size = 'small',
              family='serif')
    ax.set_ylabel(r'N',
              size = 'small',
              family='serif')
    ax.set_title(r'1f: Luminosity function of magnitude limited sample')
    ax.grid(True)
    ax.legend(loc='lower left')

    if True:
        plt.savefig('hw1.1f.png',bbox_inches='tight',dpi=500)
    if show:
        fig.show()

def plot_problem_1g(apparent_mags_bin=None, N_apparent_mags=None, show=True):
    import matplotlib.pyplot as plt
    import numpy as np

    # Create figure
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)

    ax.plot(apparent_mags_bin[:-1], N_apparent_mags,
            color='g', label='Untrimmed', linestyle='-',
            drawstyle='steps-post')

    ax.set_xlim([17,20.5])
    #ax.set_yscale('log')

    # Adjust asthetics
    ax.set_xlabel(r'm$_{bj}$',
              size = 'small',
              family='serif')
    ax.set_ylabel(r'N',
              size = 'small',
              family='serif')
    ax.set_title(r'1g: Number of galaxies in volume-limited sample')
    ax.grid(True)

    if True:
        plt.savefig('hw1.1g.png',bbox_inches='tight',dpi=500)
    if show:
        fig.show()

def plot_problem_1h(m_array_distances=None, distances_cut=None, show=True):
    import matplotlib.pyplot as plt
    import numpy as np
    # Create figure
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)

    ax.scatter(m_array_distances, distances_cut/1e6,
            color='k', label='Untrimmed', alpha=0.1)

    #ax.set_xlim([-13.1,-22.9])
    #ax.set_yscale('log')

    # Adjust asthetics
    ax.set_xlabel(r'M [magnitudes]',
              size = 'small',
              family='serif')
    ax.set_ylabel(r'Distance [Mpc]',
              size = 'small',
              family='serif')
    ax.set_title(r'1h: Distance vs. M in volume-limited sample')
    ax.grid(True)

    if True:
        plt.savefig('hw1.1h.png',bbox_inches='tight',dpi=500)
    if show:
        fig.show()

def main():
    import numpy as np

    N = 10**5
    show = False

    #problem_1b(show=show)

    problem_1c(N=N, show=show)

    problem_1d(N=N, show=show)

    problems_1fgh(N=N, show=show)

if __name__ == '__main__':
    main()
