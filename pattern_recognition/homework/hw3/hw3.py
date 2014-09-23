#!/usr/bin/python

import numpy as np

def chunks(data, chunk_size):

    l = data
    n = chunk_size

    if n < 1:
        n = 1
    #return [l[i:i + n] for i in range(0, len(l), n)]
    for i in xrange(0, len(l), n):
        yield l[i:i+n]

def chunk(input, size):
    return map(None, *([iter(input)] * size))

def split_seq(iterable, size):

    import itertools

    it = iter(iterable)
    item = list(itertools.islice(it, size))
    while item:
        yield item
        item = list(itertools.islice(it, size))
'''
Problems
'''

def problem_1():

    A = np.array([[1, 1],
                  [1, -1],
                  [1, 1]])

    b = [1, 1, 0]

    x = np.linalg.lstsq(A, b)

    print 'x = ', x[0]

def problem_2c():

    import matplotlib.pyplot as plt

    m = 1000

    a = np.random.uniform(0, 1, size=m)
    A = np.matrix([np.ones((m)), a]).T

    b = a**2 + np.random.normal(0, 1, size=m)
    B = np.matrix(b).T

    # Fit x
    x = np.linalg.inv(A.T * A) * A.T * B

    # create b_fit
    a_fit = np.arange(0,1,0.01)
    b_fit = x[0, 0] + x[1, 0] * a_fit**2

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
    ax.plot(a_fit, b_fit, color='r')
    ax.set_xlabel(r'$a_i$')
    ax.set_ylabel(r'$b_i$')

    plt.savefig('problem2c_fig.png', bbox_inches='tight')

def problem_2d():

    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import ImageGrid

    noises = np.arange(0.1, 1.9, 0.3)
    a_list = []
    b_list = []
    a_fit_list = []
    b_fit_list = []

    for i, noise in enumerate(noises):
        m = 1000

        a = np.random.uniform(0, 1, size=m)
        A = np.matrix([np.ones((m)), a]).T

        b = a**2 + np.random.normal(0, noise, size=m)
        B = np.matrix(b).T

        # Fit x
        x = np.linalg.inv(A.T * A) * A.T * B

        # create b_fit
        a_fit = np.arange(0,1,0.01)
        b_fit = x[0, 0] + x[1, 0] * a_fit**2

        a_list.append(a)
        b_list.append(b)
        a_fit_list.append(a_fit)
        b_fit_list.append(b_fit)

    # Plot!
    # Set up plot aesthetics
    plt.clf()
    plt.rcdefaults()
    colormap = plt.cm.gist_ncar
    #color_cycle = [colormap(i) for i in np.linspace(0, 0.9, len(flux_list))]
    font_scale = 8
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
    fig = plt.figure(figsize=(7,2))

    imagegrid = ImageGrid(fig, (1,1,1),
                 nrows_ncols=(1, len(a_list)),
                 ngrids=len(a_list),
                 axes_pad=0,
                 aspect=False,
                 label_mode='L',
                 share_all=True)


    for i in xrange(len(a_list)):

        a = a_list[i]
        b = b_list[i]
        a_fit = a_fit_list[i]
        b_fit = b_fit_list[i]

        ax = imagegrid[i]
        #ax = fig.add_subplot(100 + len(a_list)*10 + i)
        ax.plot(a, b, linestyle='', marker='^', markersize=3, alpha=0.4)
        ax.plot(a_fit, b_fit, color='r')
        ax.set_xlabel(r'$a_i$')
        ax.set_ylabel(r'$b_i$')
        ax.set_xlim([0.01, 0.99])

        ax.annotate(r'$\sigma$ = {0:.2f}'.format(noises[i]),
                xytext=(0.5, 0.85),
                xy=(0.1, 0.9),
                textcoords='axes fraction',
                xycoords='axes fraction',
                color='k',
                fontsize=7,
                bbox=dict(boxstyle='round',
                          facecolor='w',
                          alpha=0.3),
                horizontalalignment='right',
                verticalalignment='bottom',

                )

    plt.savefig('problem2d_fig.png', bbox_inches='tight')

def problem_4():

    import scipy.io
    data = scipy.io.loadmat('face_emotion_data.mat')

    A = np.matrix(data['A'])
    b = np.matrix(data['b'])

    x = np.linalg.lstsq(A, b)[0]

    # 4a answer
    print '4a'
    print x
    print ''

    # 4e
    A_chunked = np.asarray(chunk(np.asarray(A), 8))
    b_chunked = np.asarray(chunk(np.asarray(b), 8))

    false_positives = 0.0
    false_negatives = 0.0

    for i in xrange(0, 8):
        indices_1 = np.arange(0, 8, 1).tolist()
        del indices_1[i]

        indices_1 = np.asarray(indices_1)

        # Exclude 8 elements
        A_sub = A_chunked[:, indices_1, :]
        b_sub = b_chunked[:, indices_1]

        # Reshape to original format
        A_sub = np.reshape(A_sub,
                           (A_sub.shape[0]*A_sub.shape[1], A_sub.shape[2]))
        b_sub = np.reshape(b_sub,
                           (b_sub.shape[0]*b_sub.shape[1], b_sub.shape[2]))

        # Convert arrays to matrices
        A_sub = np.matrix(A_sub)
        b_sub = np.matrix(b_sub)

        # Solve for weights
        x = np.linalg.lstsq(A_sub, b_sub)[0]

        # Check error of weights on holdout set
        A_holdout = A_chunked[:, i, :]
        b_holdout = b_chunked[:, i]
        A_holdout = np.matrix(A_holdout)
        b_holdout = np.matrix(b_holdout)

        estimates = A_holdout * x

        # Check for positive errors
        false_positives += len(b_holdout[estimates > 0] < 0) / 16.0
        false_negatives += len(b_holdout[estimates < 0] > 0) / 16.0

    false_positives /= 8.0
    false_negatives /= 8.0

    print '4f'
    print 'Errors with 9 features:'
    print 'Error on positives, error on negatives', \
          false_positives, false_negatives

    # Perform with only 3 features
    A_chunked = np.asarray(chunk(np.asarray(A[:, (0, 3, 8)]), 8))
    b_chunked = np.asarray(chunk(np.asarray(b), 8))

    false_positives = 0.0
    false_negatives = 0.0

    for i in xrange(0, 8):
        indices_1 = np.arange(0, 8, 1).tolist()
        del indices_1[i]

        indices_1 = np.asarray(indices_1)

        # Exclude 8 elements
        A_sub = A_chunked[:, indices_1, :]
        b_sub = b_chunked[:, indices_1]

        # Reshape to original format
        A_sub = np.reshape(A_sub,
                           (A_sub.shape[0]*A_sub.shape[1], A_sub.shape[2]))
        b_sub = np.reshape(b_sub,
                           (b_sub.shape[0]*b_sub.shape[1], b_sub.shape[2]))

        # Convert arrays to matrices
        A_sub = np.matrix(A_sub)
        b_sub = np.matrix(b_sub)

        # Solve for weights
        x = np.linalg.lstsq(A_sub, b_sub)[0]

        # Check error of weights on holdout set
        A_holdout = A_chunked[:, i, :]
        b_holdout = b_chunked[:, i]
        A_holdout = np.matrix(A_holdout)
        b_holdout = np.matrix(b_holdout)

        estimates = A_holdout * x

        # Check for positive errors
        false_positives += len(b_holdout[estimates > 0] < 0) / 16.0
        false_negatives += len(b_holdout[estimates < 0] > 0) / 16.0

    false_positives /= 8.0
    false_negatives /= 8.0

    print 'Errors with 3 features:'
    print 'Error on positives, error on negatives', \
          false_positives, false_negatives

def main():
    #problem_1()
    #problem_2c()
    #problem_2d()
    problem_4()

if __name__ == '__main__':
	main()
