#!/usr/bin/python

import numpy as np

class IVP_simulation():

    def __init__(self, scheme='FTCS', boundary_type='periodic', flow_speed=1,
            x_range=(0,1), t_range=(0,3), delta_x=0.01, delta_t=0.005,
            initial_condition=0):

        self.scheme = scheme
        self.boundary_type = boundary_type
        self.flow_speed = flow_speed
        self.x_range = x_range
        self.t_range = t_range
        self.delta_x = delta_x
        self.delta_t = delta_t
        self.grid = self.__create_grid(x_range, t_range, delta_x, delta_t)
        self.x_grid = self.__calc_x_grid(x_range, delta_x)
        self.t_grid = self.__calc_time_grid(t_range, delta_t)
        self.inital_condition = \
                self.__set_initial_conditions(initial_condition)
        self.t = 0

    def __create_grid(self, x_range, t_range, delta_x, delta_t):
        x_size = (x_range[1] - x_range[0]) / delta_x
        t_size = (t_range[1] - t_range[0]) / delta_t

        return np.empty((x_size, t_size))

    def __set_initial_conditions(self, initial_condition):
        # initial condition at time = 0
        self.__set_y(initial_condition(self.x_grid), 0)

    def __get_y(self, time):
        return self.grid[:,self.__get_time_index(time)]

    def __set_y(self, x_values, time):
        self.grid[:,self.__get_time_index(time)] = x_values

    def __get_time_index(self, time):
        return np.abs(self.t_grid - time).argmin()

    def __calc_x_grid(self, x_range, delta_x):
        x_size = (self.x_range[1] - self.x_range[0]) / self.delta_x
        x_grid = np.linspace(self.x_range[0], self.x_range[1], x_size)
        return x_grid

    def __calc_time_grid(self, t_range, delta_t):
        time_size = (self.t_range[1] - self.t_range[0]) / self.delta_t
        time_grid = np.linspace(self.t_range[0], self.t_range[1], time_size)
        return time_grid

    def run_simulation(self, scheme=None, boundary_type=None, flow_speed=None,
            x_range=None, t_range=None, delta_x=None, delta_t=None,
            initial_condition=None):

        '''
        Parameters
        ----------
        scheme : str
            Scheme for time steps. Options are:
                'FTCS' : Forward Time Center Space
                'FTBS' : Forward Time Backward Space
                'FTFS' : Forward Time Forward Space
        boundary_type : str
            Type of boundary condition solution. Options are:
                'periodic' : periodic
        flow_speed : float
            The flow speed.
        x_range : tuple
            Spatial range.
        t_range : tuple
            Time range.
        delta_x : float
            Spatial resolution.
        delta_t : float
            Time resolution.

        '''

        # Check for set parameters
        change_grid = False
        if scheme is not None:
            self.scheme = scheme
        if boundary_type is not None:
            self.boundary_type = boundary_type
        if flow_speed is not None:
            self.flow_speed = flow_speed
        if x_range is not None:
            self.x_range = x_range
            change_grid = True
        if t_range is not None:
            self.t_range = t_range
            change_grid = True
        if delta_x is not None:
            self.delta_x = delta_x
            change_grid = True
        if delta_t is not None:
            self.delta_t = delta_t
            change_grid = True
        if change_grid:
            self.grid = __create_grid(self.x_range, self.t_range,
                self.delta_x, self.delta_t)
            self.x_grid = self.__calc_x_grid(x_range, delta_x)
            self.t_grid = self.__calc_time_grid(t_range, delta_t)
        if initial_condition is not None:
            self.initial_condition = initial_condition

        # Complete simulation
        if self.scheme.lower() == 'ftcs':
            for i in range(len(self.t_grid)):
                self.step_FTCS()
        elif self.scheme.lower() == 'ftbs':
            for i in range(len(self.t_grid)):
                self.step_FTBS()
        elif self.scheme.lower() == 'ftfs':
            for i in range(len(self.t_grid)):
                self.step_FTFS()
        elif self.scheme.lower() == 'btcs':
            for i in range(len(self.t_grid)):
                self.step_BTCS()

    def step_FTCS(self):

        # get y values at all x at one time
        y = self.__get_y(self.t)

        # initialize next time step
        y_tp1 = np.empty((y.shape))

        for k in range(len(self.x_grid)):
            if k == len(self.x_grid) - 1:
                y_kp1 = y[0]
                x_kp1 = self.x_grid[k] + self.delta_x
                y_km1 = y[k - 1]
                x_km1 = self.x_grid[k - 1]
            elif k == 0:
                y_km1 = y[len(self.x_grid) - 1]
                x_km1 = self.x_grid[k] - self.delta_x
                y_kp1 = y[k + 1]
                x_kp1 = self.x_grid[k + 1]
            else:
                y_kp1 = y[k + 1]
                y_km1 = y[k - 1]
                x_kp1 = self.x_grid[k + 1]
                x_km1 = self.x_grid[k - 1]

            y_tp1[k] = y[k] - self.flow_speed * (self.delta_t) * \
                    (y_kp1 - y_km1) / (x_kp1 - x_km1)

        # Set next time step y values, and increase time
        self.__set_y(y_tp1, self.t + self.delta_t)
        self.t += self.delta_t

    def step_FTBS(self):
        # get y values at all x at one time
        y = self.__get_y(self.t)

        # initialize next time step
        y_tp1 = np.empty((y.shape))

        for k in range(len(self.x_grid)):
            if k == len(self.x_grid) - 1:
                y_kp1 = y[0]
                y_k = y[k]
                y_km1 = y[k - 1]
            elif k == 0:
                y_km1 = y[-1]
                y_k = y[k]
                y_kp1 = y[k + 1]
            else:
                y_kp1 = y[k + 1]
                y_km1 = y[k - 1]
                y_k = y[k]

            y_tp1[k] = y[k] - self.flow_speed * (self.delta_t) * \
                    (y_k - y_km1) / self.delta_x

        # Set next time step y values, and increase time
        self.__set_y(y_tp1, self.t + self.delta_t)
        self.t += self.delta_t

    def step_FTFS(self):
        # get y values at all x at one time
        y = self.__get_y(self.t)

        # initialize next time step
        y_tp1 = np.empty((y.shape))

        for k in range(len(self.x_grid)):
            if k == len(self.x_grid) - 1:
                y_kp1 = y[0]
                x_kp1 = self.x_grid[k] + self.delta_x
                y_k = y[k]
                x_k = self.x_grid[k]
                y_km1 = y[k - 1]
                x_km1 = self.x_grid[k - 1]
            elif k == 0:
                y_km1 = y[len(self.x_grid) - 1]
                x_km1 = self.x_grid[k] - self.delta_x
                y_k = y[k]
                x_k = self.x_grid[k]
                y_kp1 = y[k + 1]
                x_kp1 = self.x_grid[k + 1]
            else:
                y_kp1 = y[k + 1]
                y_km1 = y[k - 1]
                y_k = y[k]
                x_k = self.x_grid[k]
                x_kp1 = self.x_grid[k + 1]
                x_km1 = self.x_grid[k - 1]

            y_tp1[k] = y[k] - self.flow_speed * (self.delta_t) * \
                    (y_kp1 - y_k) / (x_kp1 - x_k)

        # Set next time step y values, and increase time
        self.__set_y(y_tp1, self.t + self.delta_t)
        self.t += self.delta_t

    def step_BTCS(self):

        ''' Backward-time Center-Space scheme. Solves using linear algebra.
        '''
        # get y values at all x at one time
        y = self.__get_y(self.t)

        # initialize next time step
        y_tp1 = np.empty((y.shape))

        K = len(self.x_grid)
        M = np.zeros((K,K))
        alpha = self.flow_speed * self.delta_t / self.delta_x

        for i in range(K):
            for j in range(K):
                if i == j:
                    M[i,j] = 1
                    if i==0:
                        M[i+1,j] = -alpha
                        M[-1,j] = alpha
                    elif i==K-1:
                        M[0,j] = -alpha
                        M[i-1,j] = alpha
                    else:
                        M[i+1,j] = -alpha
                        M[i-1,j] = alpha

        #
        y_tp1 = np.linalg.solve(M, y)

        # Set next time step y values, and increase time
        self.__set_y(y_tp1, self.t + self.delta_t)
        self.t += self.delta_t

        return None

    def plot_slice(self, times, limits=None, savedir='./', filename=None,
            show=True, title='', additional_sims=None, additional_labels=None):

        ''' Plots
        '''

        # Import external modules
        import numpy as np
        import math
        import pyfits as pf
        import matplotlib.pyplot as plt
        import matplotlib
        from mpl_toolkits.axes_grid1 import ImageGrid

        # Set up plot aesthetics
        plt.clf()
        plt.rcdefaults()
        fontScale = 10
        params = {#'backend': 'png',
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
                      nrows_ncols=(2,len(times)/2),
                      ngrids = len(times),
                      direction='row',
                      axes_pad=1,
                      aspect=False,
                      share_all=True,
                      label_mode='All')

        x_analytic = np.linspace(self.x_range[0],self.x_range[1],1e5)
        def calc_y(x, time):
            return np.sin(2*np.pi * (x - self.flow_speed * time))

        # save current grid if additional ones being plotted
        grid_save = self.grid

        markers = ['s','o','*']

        for i, time in enumerate(times):
            ax = grid[i]

            if additional_sims is None:
                ax.plot(self.x_grid, self.__get_y(time),
                        color='k',
                        markersize=3,
                        marker='o',
                        linestyle='None',
                        label='Numerical')
            elif additional_sims is not None:
                ax.plot(self.x_grid, self.__get_y(time),
                        #color='k',
                        markersize=3,
                        marker='^',
                        linestyle='None',
                        label='Numerical: %s' % additional_labels[0])
                for j, sim in enumerate(additional_sims):
                    ax.plot(sim.x_grid, sim.__get_y(time),
                            #color='k',
                            markersize=3,
                            marker=markers[j],
                            linestyle='None',
                            label='Numerical: %s' % additional_labels[j+1])

            # plot analytic solution
            y_analytic = calc_y(x_analytic, time)
            ax.plot(x_analytic, y_analytic,
                    color='k',
                    linestyle='-',
                    label='Analytic')

            if limits is not None:
                ax.set_xlim(limits[0],limits[1])
                ax.set_ylim(limits[2],limits[3])

            # Adjust asthetics
            ax.set_xlabel('x',)
            ax.set_ylabel(r'y(x, %s)' % time)
            ax.grid(True)
            ax.legend(loc='lower left')
            ax.set_title('t = %s' % time)

        # reset grid
        self.grid = grid_save

        if filename is not None:
            plt.savefig(savedir + filename,bbox_inches='tight')
        if show:
            fig.show()

def plot_amplitude(alpha_array = (-1,1), G_values = (1), amp_functions = None,
        limits=None, savedir='./', filename=None, show=True, title=''):

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
        params = {#'backend': 'png',
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
                      nrows_ncols=(2,1),
                      ngrids = 2,
                      direction='row',
                      axes_pad=1,
                      aspect=False,
                      share_all=True,
                      label_mode='All')

        colors = ['k','b','g','r','c']
        linestyles = ['-','--','-.','-','-']
        letters = ['a','b']

        for i, amp_function in enumerate(amp_functions):
            ax = grid[i]

            for j, G in enumerate(G_values):
                ax.plot(alpha_array, amp_function(alpha_array, G),
                        color = colors[j],
                        label = 'G = %s' % G,
                        linestyle = '-')

            if limits is not None:
                ax.set_xlim(limits[0],limits[1])
                ax.set_ylim(limits[2],limits[3])

            # Adjust asthetics
            ax.set_xlabel(r'$\alpha$',)
            ax.set_ylabel(r'Amplitude')
            ax.annotate('(%s)' % letters[i],
                    xy = (0.9, 0.1),
                    xycoords='axes fraction',
                    textcoords='axes fraction')
            ax.grid(True)
            ax.legend(loc='lower left')

        if filename is not None:
            plt.savefig(savedir + filename,bbox_inches='tight')
        if show:
            fig.show()

def problem_1():
    def initial_condition(x):
        return np.sin(2 * np.pi * x)

    ftcs_sim_ftbs = IVP_simulation(scheme = 'FTCS',
            boundary_type = 'periodic',
            flow_speed = 1,
            x_range = (0,1),
            t_range = (0,3),
            delta_x = 0.01,
            delta_t = 0.005,
            initial_condition = initial_condition)

    ftcs_sim_ftbs.run_simulation()

    savedir = '/home/elijah/class_2014_spring/fluids/midterm/'
    savedir = '/usr/users/ezbc/Desktop/fluids/midterm/'
    times = [0, 1, 2, 3,]
    ftcs_sim_ftbs.plot_slice(times, savedir=savedir,
                filename='q1_ftcs.png',
                title = 'FTCS simulation slices',
                show=False)

def problem_2():
    def initial_condition(x):
        return np.sin(2 * np.pi * x)

    ftbs_sim = IVP_simulation(scheme = 'FTBS',
            boundary_type = 'periodic',
            flow_speed = 1,
            x_range = (0,1),
            t_range = (0,3),
            delta_x = 0.01,
            delta_t = 0.005,
            initial_condition = initial_condition)

    ftfs_sim = IVP_simulation(scheme = 'FTFS',
            boundary_type = 'periodic',
            flow_speed = 1,
            x_range = (0,1),
            t_range = (0,3),
            delta_x = 0.01,
            delta_t = 0.005,
            initial_condition = initial_condition)

    ftbs_sim.run_simulation()
    ftfs_sim.run_simulation()

    savedir = '/home/elijah/class_2014_spring/fluids/midterm/'
    savedir = '/usr/users/ezbc/Desktop/fluids/midterm/'

    times = [0, 1, 2, 3,]
    ftbs_sim.plot_slice(times, savedir=savedir,
                filename='q2_ftbs.png',
                title = 'FTBS simulation slices',
                show=False)
    ftfs_sim.plot_slice(times, savedir=savedir,
                filename='q2_ftfs.png',
                title = 'FTFS simulation slices',
                show=False)

def problem_3():
    def initial_condition(x):
        return np.sin(2 * np.pi * x)

    ftbs_sim = IVP_simulation(scheme = 'FTBS',
            boundary_type = 'periodic',
            flow_speed = -1,
            x_range = (0,1),
            t_range = (0,3),
            delta_x = 0.01,
            delta_t = 0.005,
            initial_condition = initial_condition)

    ftfs_sim = IVP_simulation(scheme = 'FTFS',
            boundary_type = 'periodic',
            flow_speed = -1,
            x_range = (0,1),
            t_range = (0,3),
            delta_x = 0.01,
            delta_t = 0.005,
            initial_condition = initial_condition)

    ftbs_sim.run_simulation()
    ftfs_sim.run_simulation()

    savedir = '/home/elijah/class_2014_spring/fluids/midterm/'
    savedir = '/usr/users/ezbc/Desktop/fluids/midterm/'
    times = [0, 1, 2, 3,]
    ftbs_sim.plot_slice(times, savedir=savedir,
                filename='q3_ftbs.png',
                title = 'FTBS simulation slices',
                show=False)
    ftfs_sim.plot_slice(times, savedir=savedir,
                filename='q3_ftfs.png',
                title = 'FTFS simulation slices',
                show=False)

def problem_5():

    def amp1(alpha, G):
        return (1 + 2*G*alpha*(alpha + 1))**0.5

    def amp2(alpha, G):
        return (1 + 2*G*alpha*(alpha - 1))**0.5

    amp_tuple = (amp1, amp2)
    alpha_array = np.linspace(-3, 3, 1e4)
    G_values = (0.1, 0.3, 0.5, 0.7, 0.9)

    savedir = '/usr/users/ezbc/Desktop/fluids/midterm/'

    plot_amplitude(alpha_array = alpha_array, G_values = G_values,
            amp_functions = amp_tuple,
            limits = [-2, 2, 0.7, 1.2],
            savedir = savedir,
            filename = 'q5.png',
            title = 'FTFS and FTBS Amplitudes',
            show=False)

def problem_6():
    def initial_condition(x):
        return np.sin(2 * np.pi * x)

    ftbs_sim = IVP_simulation(scheme = 'FTBS',
            boundary_type = 'periodic',
            flow_speed = 1,
            x_range = (0,1),
            t_range = (0,3),
            delta_x = 0.01,
            delta_t = 0.02,
            initial_condition = initial_condition)

    ftbs_sim.run_simulation()

    savedir = '/home/elijah/class_2014_spring/fluids/midterm/'
    savedir = '/usr/users/ezbc/Desktop/fluids/midterm/'
    times = [0, 1, 2, 3,]
    ftbs_sim.plot_slice(times, savedir=savedir,
                filename='q6_ftbs.png',
                title = 'FTBS simulation slices',
                show=False)

def problem_7():
    def initial_condition(x):
        return np.sin(2 * np.pi * x)

    alphas = (0.25, 0.5, 0.75)
    delta_x = 0.01
    flow_speed = 1
    x_range = (0,1)
    t_range = (0,3)

    sim_list = []
    additional_labels = []

    savedir = '/home/elijah/class_2014_spring/fluids/midterm/'
    #savedir = '/usr/users/ezbc/Desktop/fluids/midterm/'
    times = [0, 1, 2, 3,]

    for i, alpha in enumerate(alphas):
        delta_t = np.abs(alpha * delta_x)

        ftbs_sim = IVP_simulation(scheme = 'FTBS',
                boundary_type = 'periodic',
                flow_speed = flow_speed,
                x_range = x_range,
                t_range = t_range,
                delta_x = delta_x,
                delta_t = delta_t,
                initial_condition = initial_condition)

        ftbs_sim.run_simulation()

        sim_list.append(ftbs_sim)

        additional_labels.append(r'$\alpha$ = %s' % alpha)

    sim_list[0].plot_slice(times, savedir=savedir,
                filename='q7.png',
                title = 'FTBS simulation slices',
                show=False,
                additional_sims = sim_list[1:],
                additional_labels = additional_labels)

def problem_8():
    def initial_condition(x):
        return np.sin(2 * np.pi * x)

    alphas = (0.1, 0.5, 1, 2)
    delta_x = 0.01
    flow_speed = 1
    x_range = (0,1)
    t_range = (0,3)

    sim_list = []
    additional_labels = []

    savedir = '/home/elijah/class_2014_spring/fluids/midterm/'
    #savedir = '/usr/users/ezbc/Desktop/fluids/midterm/'
    times = [0, 1, 2, 3,]

    for i, alpha in enumerate(alphas):
        delta_t = np.abs(alpha * delta_x)

        btcs_sim = IVP_simulation(scheme = 'BTCS',
                boundary_type = 'periodic',
                flow_speed = flow_speed,
                x_range = x_range,
                t_range = t_range,
                delta_x = delta_x,
                delta_t = delta_t,
                initial_condition = initial_condition)

        btcs_sim.run_simulation()

        sim_list.append(btcs_sim)

        additional_labels.append(r'$\alpha$ = %s' % alpha)

    sim_list[0].plot_slice(times, savedir=savedir,
                filename='q8.png',
                title = 'BTCS simulation slices',
                show=False,
                additional_sims = sim_list[1:],
                additional_labels = additional_labels)

def main():
    #problem_1()
    #problem_2()
    #problem_3()
    #problem_5()
    #problem_6()
    #problem_7()
    problem_8()

if __name__ == '__main__':
    main()




