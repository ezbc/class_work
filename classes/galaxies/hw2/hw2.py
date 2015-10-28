#!/usr/bin/python

def calc_g(omega):

    return 5/2. * omega * (1/70. + 209*omega/140. - omega**2 / 140 +
            omega**(4/7.))**-1

def calc_omega(omega_0, z):
    '''
    '''

    return omega_0 * (1 + z)**3 / (1 - omega_0 + (1 + z)**3 * omega_0)

def problem_1a():

    g_0 = calc_g(calc_omega(1, 0))
    g_2 = calc_g(calc_omega(1, 2))

    g = g_0 / g_2

    print('Factor of growth of perturbations where Omega_m = 1: %s' % g)

    g_0 = calc_g(calc_omega(.30, 0))
    g_2 = calc_g(calc_omega(.30, 2))

    g = g_0 / g_2

    print('Factor of growth of perturbations Omega_m = 0.3: %s' % g)

def plot_g():

    import matplotlib.pyplot as plt
    import numpy as np

    # Create figure
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)

    z_array = np.logspace(-4,4,1000)

    ax.plot(z_array,calc_omega(0.3,z_array))

    #ax.set_ylim([-6,-1])
    #ax.set_xlim([1e-4,1e3])
    ax.set_yscale('log')
    ax.set_xscale('log')

    # Adjust asthetics
    ax.set_xlabel(r'g(z)',
              size = 'small',
              family='serif')
    ax.set_ylabel(r'z',
              size = 'small',
              family='serif')
    ax.set_title('1b')
    ax.grid(True)
    ax.legend(loc='lower left')

    if True:
        plt.savefig('hw2_1b.png',bbox_inches='tight')
    if True:
        fig.show()

def problem_1b():

    plot_g()

    g_2 = calc_g(calc_omega(.30, 2))
    g_1000 = calc_g(calc_omega(.30, 1000))

    g = g_2 / g_1000

    print('Factor of growth of perturbations where Omega_m = 1: %s' % g)

def problem_2b():

    import numpy as np

    H_0 = 100 * 10**5 / (10**6 * 3e18) # 1/s
    G = 6.67e-8

    p_0 = 3 * H_0**2 / (8 * np.pi * G)

    denom = (p_0 * 178)

    print('H_0: %s' % H_0)
    print('Denominator: %s' % denom)

def problem_2c():

    import numpy as np

    H_0 = 100 * 10**5 / (10**6 * 3e18) # 1/s
    G = 6.67e-8

    p_0 = 3 * H_0**2 / (8 * np.pi * G)

    denom = (p_0 * 178)

    R = 3e23
    v = 220e5

    p_vir = 3 * v**2 / (4*np.pi * G * R**2)

    z = (p_vir / denom)**(1/3.)

    print('p_vir: %s' % p_vir)
    print('z: %s' % z)

def problem_3a():

    from scipy.integrate import quad as integrate

    def imf_N(mass=0, alpha=2.35):
        return mass**-alpha

    alpha = 2.35

    coeff = 25. / integrate(imf_N, 20, 100, args=(alpha))[0]

    def imf_M(mass=0, alpha=-2.35):
        return mass**-alpha * mass

    mass = coeff*integrate(imf_M, 0.1, 100, args=(alpha))[0]

    print('Mass of cluster: %s' % mass)

def problem_3b():

    import numpy as np
    from scipy.integrate import quad as integrate
    from sympy.solvers import solve

    def imf_N(mass=0, alpha=2.35):
        if mass <= 1:
            return (1/mass)*np.exp(-(np.log10(mass) -
                np.log10(0.22))**2 / (2*0.57**2))
        elif mass > 1:
            return mass**-alpha

    alpha = 2.35

    coeff_high = 25. / integrate(imf_N, 20, 100, args=(alpha))[0]

    mass = 1
    coeff_low = coeff_high * mass**-alpha / \
                ((1/mass)*np.exp(-(np.log10(mass) -
                np.log10(0.22))**2 / (2*0.57**2)))

    def imf_M(mass=0, alpha=-2.35):
        if mass <= 1:
            return coeff_low*(mass/mass)*np.exp(-(np.log10(mass) -
                np.log10(0.22))**2 / (2*0.57**2))
        elif mass > 1:
            return coeff_high * mass**-alpha * mass

    mass = integrate(imf_M, 0.1, 100, args = (alpha))[0]

    print('Mass of cluster: %s' % mass)

def problem_3c():

    import numpy as np
    from scipy.integrate import quad as integrate


    # Compute Salpeter mass
    def imf_N(mass=0, alpha=2.35):
        return mass**-alpha

    alpha = 2.35

    coeff = 30. / integrate(imf_N, 20, 100, args=(alpha))[0]

    def imf_M(mass=0, alpha=-2.35):
        return mass**-alpha * mass

    mass_a = coeff*integrate(imf_M, 0.1, 100, args=(alpha))[0]


    # Compute Chabrier mass
    def imf_N(mass=0, alpha=2.35):
        if mass <= 1:
            return (1/mass)*np.exp(-(np.log10(mass) -
                np.log10(0.22))**2 / (2*0.57**2))
        elif mass > 1:
            return mass**-alpha
    alpha = 2.35

    coeff_high = 25. / integrate(imf_N, 20, 100, args=(alpha))[0]

    mass = 1
    coeff_low = coeff_high * mass**-alpha / \
                ((1/mass)*np.exp(-(np.log10(mass) -
                np.log10(0.22))**2 / (2*0.57**2)))

    def imf_M(mass=0, alpha=-2.35):
        if mass <= 1:
            return coeff_low*(mass/mass)*np.exp(-(np.log10(mass) -
                np.log10(0.22))**2 / (2*0.57**2))
        elif mass > 1:
            return coeff_high * mass**-alpha * mass

    mass_b = integrate(imf_M, 0.1, 100, args = (alpha))[0]


    # now calculate velocity dipsersion
    min_vel_disp = mass_a - mass_b

    print('Min velocity resolution: %s' % min_vel_disp)

def main():

    print('\nProblem 1a:')
    problem_1a()

    print('\nProblem 1b:')
    problem_1b()

    print('\nProblem 2b:')
    problem_2b()

    print('\nProblem 2c:')
    problem_2c()

    print('\nProblem 3a:')
    problem_3a()

    print('\nProblem 3b:')
    problem_3b()

    print('\nProblem 3c:')
    problem_3c()

if __name__ == '__main__':
	main()

