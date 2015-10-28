# Elijah test

import pickle
import sys
import matplotlib.pyplot as plt


import os
os.system('export PYTHONPATH=$PYTHONPATH:/home/ezbc/classes/machine_learning/project/gausspy_module/gausspy/')
import gausspy

import gausspy.gp

#results = pickle.load(open('data/spline_fits.pickle'))

path = 'data/HT2003_data_test100.pickle'
data = pickle.load(open(path))

g = gausspy.gp.GaussianDecomposer()
g.load_training_data(path)

print(g.training_data)

g.p['SNR_thresh'] = 5.
g.set('phase', 'two')
g.set('alpha1', -0.88+2.0)
g.set('alpha2', 0.73+2.0)

g.train(alpha1_initial = -0.88+2.0, alpha2_initial=0.73+2.0, plot=False, verbose=False, mode='c', eps=0.)
#g.batch_decomposition()




