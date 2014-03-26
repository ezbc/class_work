#!/usr/bin/python

def calc_N_photons(mass):

    ''' Calculates number of photons emitted from a star with mass mass, in a
    year.

    Parameters
    ----------
    mass : float
        mass of star in M_sun

    Returns
    -------
    N_L : float
        Number of photons emitted from a star with mass mass, in a
        year.

    '''

    # Import external modules
    import numpy as np

    # Piecewise function
    if (mass > 18) and (mass < 41):
        log_N_L = 44.129 + 3.276 * np.log10(mass)
    elif mass >= 41:
        log_N_L = 47.170 + 1.417 * np.log10(mass)
    else:
        log_N_L = 1e-3000000000000000000000

    return 10**log_N_L

def calc_star_lifetime(mass):

    ''' Calculates lifetime of star in years.

    Parameters
    ----------
    mass : float
        Mass of star in M_sun.

    Returns
    -------
    t_ms : float
        Lifetime of star in years.

    '''

    # import external modules
    import numpy as np

    # piecewise function
    if (mass > 15) and (mass < 33):
        log_t_ms = 8.475 -1.160 * np.log10(mass)
    elif mass >= 33:
        log_t_ms = 7.472 - 0.501 * np.log10(mass)
    else:
        log_t_ms = 1e-30000000000000000000000000000

    sec2yr = 365 * 24 * 3600.

    return 10**log_t_ms

def calc_N_stars(mass=18, Gamma=-1.35):

    ''' Calculates number of stars of given mass in IMF.

    Parameters
    ----------
    mass : float
        Mass of star in M_sun.
    Gamma : float
        Gamma in IMF

    Returns
    -------
    N_stars : float
        Number of stars / normalization factor.

    '''

    N_stars = mass**(Gamma - 1)

    return N_stars

def calc_avg_N_photons(mass_limits=(18,100), Gamma=-1.35):

    ''' Calculates number of photons produced in a stars lifetime in a given
    mass range of the IMF.

    Parameters
    ----------
    mass_limits : tuple
        (lower limit, upper limit) in M_sun.
    Gamma : float
        Gamma in IMF

    Returns
    -------
    N_photons_avg : float
        Number of photons produced in a stars lifetime in a given
        mass range of the IMF.

    '''

    # import external modules
    import numpy as np
    from scipy.integrate import quad as integrate

    def integrand(mass=0, Gamma=Gamma):
        return calc_N_photons(mass) * calc_star_lifetime(mass)\
                * calc_N_stars(mass=mass, Gamma=Gamma)

    N_photons_total = integrate(integrand,
            mass_limits[0], mass_limits[1],
            args=(Gamma))[0]

    N_stars = integrate(calc_N_stars,
            mass_limits[0], mass_limits[1],
            args=(Gamma))[0]

    N_photons_avg = N_photons_total / N_stars * 365 * 24 * 3600

    return N_photons_avg

def calc_sfr_coeff(mass_limits_ionizing=(18,100), mass_limits_total=(0.1,100),
        Gamma=-1.35, N_photons_avg=2.52e63):

    ''' Calculates coefficient to convert Halpha luminosity into a SFR.

    Parameters
    ----------
    mass_limits_ionizing : tuple
        (lower limit, upper limit) in M_sun of stars creating Halpha.
    mass_limits_total : tuple
        (lower limit, upper limit) in M_sun of all stars.
    Gamma : float
        Gamma in IMF.
    N_photons_avg : float
        Average number of photons per star emitted in one year.

    Returns
    -------
    sfr_coeff : float
        SFR (M_sun) = sfr_coeff * Halpha luminosity

    '''

    # import external modules
    import numpy as np
    from scipy.integrate import quad as integrate


    N_stars_ionizing = integrate(calc_N_stars,
            mass_limits_ionizing[0], mass_limits_ionizing[1],
            args=(Gamma))[0]

    print('N %s' % N_stars_ionizing)

    # Number of ionizing photons per second by massive stars
    efficiency = 2/3.
    sec2yr = 365 * 24 * 3600.
    n = 7.37e11 * efficiency / N_photons_avg * sec2yr

    sfr_coeff = n \
            / (mass_limits_ionizing[1]**Gamma - \
                    mass_limits_ionizing[0]**Gamma) \
            * (Gamma / (Gamma + 1)) \
            * (mass_limits_total[1]**(Gamma+1) - \
                    mass_limits_total[0]**(Gamma+1))

    return sfr_coeff

def calc_energy_distribution():

    ''' Calculates energy distribution of Halpha luminosities

    Parameters
    ----------
    mass_limits_ionizing : tuple
        (lower limit, upper limit) in M_sun of stars creating Halpha.
    mass_limits_total : tuple
        (lower limit, upper limit) in M_sun of all stars.
    Gamma : float
        Gamma in IMF.
    N_photons_avg : float
        Average number of photons per star emitted in one year.

    Returns
    -------
    sfr_coeff : float
        SFR (M_sun) = sfr_coeff * Halpha luminosity

    '''

def calc_luminosities(mass_limits = (0.1, 100), binsize = 0.1,
        Gamma=-1.35):

    ''' Calculates Halpha luminosities of stars in a Salpeter IMF

    Parameters
    ----------
    mass_limits : tuple
        (lower limit, upper limit) in units of M_sun.
    binsize : float
        Binsize to create mass bins.

    Returns
    -------
    luminosities : array-like
        Halpha luminosities of each mass bin.
    masses : array-like
        Mass bins for luminosities.

    '''

    # Import external modules
    import numpy as np

    N_bins = (mass_limits[1] - mass_limits[0]) / binsize

    masses = np.linspace(mass_limits[0], mass_limits[1], N_bins)

    luminosities = np.zeros(masses.shape)

    efficiency = 2/3.

    for i in range(len(masses) - 1):
        luminosities[i] = calc_N_stars(mass = masses[i], Gamma=Gamma) \
                       * efficiency \
                       * calc_avg_N_photons(mass_limits = (masses[i],
                           masses[i] + binsize), Gamma = Gamma) \
                       * binsize \
                       / 7.37e11

    return luminosities, masses

def calc_energies(masses = None, luminosities = None, Gamma=-1.35):

    ''' Calculates energies given luminosities and masses.

    Parameters
    ----------
    luminosities : array-like
        Halpha luminosities of each mass bin.
    masses : array-like
        Mass bins for luminosities.
    Gamma : float
        Gamma in IMF.

    Returns
    -------
    energies : array-like
        Halpha energy output of stars in each mass bin.

    '''

    # Import external modules
    import numpy as np

    energies = np.zeros(masses.shape)

    # Create masses.
    for i in range(len(masses)):
        energies[i] = luminosities[i] * calc_star_lifetime(masses[i]) \
                      * calc_N_stars(mass=masses[i], Gamma=Gamma) \
                      * (masses[1] - masses[0])

    return energies

def cumulate_energies(energies = None):

    ''' Calculates cumulative energy distriubution.

    Parameters
    ----------
    energies : array-like
        Energies of each mass bin.

    Returns
    -------
    energies_cd : array-like
        Cumulative energy distribution in each mass bin.

    '''

    # Import external modules
    import numpy as np

    energies_cd = np.zeros(energies.shape)

    # Sum all elements after next element
    for i in range(len(energies)):
        energies_cd[i] = energies[i:-1].sum()

    energies_cd /= energies.sum()

    return energies_cd

def cumulate_luminosities(luminosities = None):

    ''' Calculates cumulative luminosity distriubution

    Parameters
    ----------
    luminosities : array-like
        Luminosities of each mass bin.

    Returns
    -------
    luminosities_cd : array-like
        Cumulative luminosity distribution in each mass bin.

    '''

    # Import external modules
    import numpy as np

    luminosities_cd = np.zeros(masses.shape)

    for i in range(len(masses)):
        luminosities_cd[i] = luminosities[i:-1].sum()

    luminosities_cd /= luminosities.sum()

    return luminosities_cd

def plot_energy_dist(masses=None, energy_cd=None):

    ''' Plots cumulative energy distribution of Halpha output from massive
    stars.

    Parameters
    ----------
    masses : array-like
        1-D array of masses of stars.
    energy_cd : array-like
        1-D array of cumalitve Halpha energy output for each stellar mass.

    Returns
    -------
    None

    '''

    # Import external modules
    import matplotlib.pyplot as plt
    import numpy as np

    # Create figure
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)

    # Plot data
    ax.plot(masses, energy_cd, color='k')

    # Set limits
    ax.set_ylim([0,1])
    ax.set_xlim([100,1])
    ax.set_xscale('log')

    # Adjust asthetics
    ax.set_xlabel(r'Mass (M$_\odot$)',
              size = 'large',
              family='serif')
    ax.set_ylabel(r'Cumulative Energy',
              size = 'large',
              family='serif')
    ax.set_title('4c: CDE as a Function of Mass')
    ax.grid(True)

    # Save and show figure
    if True:
        plt.savefig('hw2_4c.png',bbox_inches='tight')
    if True:
        fig.show()

def problem_4a():

    N_photons_avg = calc_avg_N_photons(mass_limits=(18,200), Gamma=-1.35)
    sfr_coeff = calc_sfr_coeff(mass_limits_ionizing=(18,200),
            mass_limits_total=(0.1,200), Gamma=-1.35,
            N_photons_avg=N_photons_avg)

    print('Number of ionizing photons per year: %s photons' % N_photons_avg)
    print('SFR: %s M_sun / yr * L_Halpha' % sfr_coeff)

def problem_4b():

    N_photons_avg = calc_avg_N_photons(mass_limits=(18,100), Gamma=-1.7)
    sfr_coeff = calc_sfr_coeff(mass_limits_ionizing=(18,100),
            mass_limits_total=(0.1,100), Gamma=-1.7,
            N_photons_avg=N_photons_avg)

    print('Number of ionizing photons per year: %s photons' % N_photons_avg)
    print('SFR: %s M_sun / yr * L_Halpha' % sfr_coeff)

def problem_4c():

    # Calculate Halpha luminosities of stars in a Salpeter IMF
    luminosities, masses = calc_luminosities(mass_limits = (0.1, 100),
            binsize = 0.01)

    # Caclulate energies from luminosities, given stellar lifetimes and a
    # Salpeter IMF
    energies = calc_energies(luminosities = luminosities,
            masses = masses)

    print('Luminosity max = %s erg s^-1' % luminosities.max())

    # Create cumulative distribution of energies
    energy_cd = cumulate_energies(energies = energies, masses = masses)

    # Finally plot
    plot_energy_dist(masses = masses, energy_cd = energy_cd)

def problem_4d():

    # Calculate Halpha luminosities of stars in a Salpeter IMF
    luminosities, masses = calc_luminosities(mass_limits = (0.1, 100),
            binsize = 0.01)

    # Create cumulative distribution of energies
    luminosities_cd = cumulate_luminosities(luminosities = luminosities,
            masses = masses)

    ionizing_mass = masses[luminosities_cd < 0.9][0]

    print('Mass above which radiates 90%% of total Halpha luminosity' \
            ' = %s M_sun' % ionizing_mass)

def main():

    print('\nProblem 4a:')
    problem_4a()

    print('\nProblem 4b:')
    problem_4b()

    print('\nProblem 4c:')
    problem_4c()

    print('\nProblem 4d:')
    problem_4d()

if __name__ == '__main__':
	main()



