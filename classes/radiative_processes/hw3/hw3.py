#!/usr/bin/python

# 2c

import numpy as np
from math import pi,exp
from scipy.integrate import quad as integrate

# SI units
e = 1.602e-19
me = 9.31e-31
c = 2.98e8
h = 6.652e-34
k = 1.38e-23
T = 5e3 * e / k
energy_low = 0.5e3 * e
energy_high = 3e3 * e
nu_low = energy_low / h
nu_high = energy_high / h
r = 1.52e22

constant = 4*pi*r**2 * 4*pi*8*e**6/(h*3*me*c**3) * (2/(pi*k*me))**0.5 * T**-0.5

def integrand(nu,T):
    return exp(-h*nu/(k*T)/nu)

result = integrate(integrand, nu_low, nu_high, args=(T))[0]

