#!/usr/bin/python

#2c)

import numpy as np
from math import pi

# constants
e = 1.602e-10 # erg
h = 6.626e-27 # erg s
c = 2.9979e10 # cm/s
me = 9.109e-28 # g
mp = 1.6726e-24 # g

# values
gp = 5.6
ge = 2.0
nu = 1420e6 # Hz

def einstein(gp=None,ge=None,nu=None):
    return (e*gp/(2*mp*c) + e*ge/(2*me*c))**2 * 8*pi**2*nu**3*h/c**3

Aif = einstein(gp=gp,ge=ge,nu=nu) # s
Aif_years = Aif * 365*24*3600

#2d)

# contsants
k = 1.38065e-16 # erg/K
mh = mp + me

# values
mtot = 1e6 * 1.99e33 # g
ntot = mtot / mh
g1 = 1
g2 = 3
T = 10e4 # K

n1 = g1*ntot * np.exp(h*nu/(k*T)) / (g2*(1+g1/g2*np.exp(h*nu/(k*T))))

n2 = ntot - n1

L = h*nu*Aif * n2








