# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 12:19:34 2017

@author: olmo
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.cosmology import Planck15
from astropy.cosmology import z_at_value
import astropy.units as u
import random
from scipy import interpolate

font = {'size': 12}

matplotlib.rc('font', **font)


templist = np.loadtxt('temperature_cooling_functions.dat')
print (templist)
lambda_zplus0pt5 = np.loadtxt('lambda_cooling_zplus0pt5.dat')
lambda_zsolar = np.loadtxt('lambda_cooling_zsolar.dat')
lambda_zminus0pt5 = np.loadtxt('lambda_cooling_zminus0pt5.dat')
lambda_zminus1 = np.loadtxt('lambda_cooling_zminus1.dat')
lambda_zminus1pt5 = np.loadtxt('lambda_cooling_zminus1pt5.dat')
lambda_zminus2 = np.loadtxt('lambda_cooling_zminus2.dat')
lambda_zminus3 = np.loadtxt('lambda_cooling_zminus3.dat')
lambda_znull = np.loadtxt('lambda_cooling_znull.dat')

zlist = np.array([0.5, 0., -0.5, -1., -1.5, -2., -3.])
new_interp_list = np.arange(-3., 0.5, 0.1)

X = zlist
Y = []

for i in range(len(templist)-1):

	X = np.concatenate((X, zlist))

for i in range(len(templist)):
	for j in range(len(zlist)):
		Y.append(templist[i])

Y = np.array((Y))

lambdaf = np.column_stack((lambda_zplus0pt5, lambda_zsolar, lambda_zminus0pt5, lambda_zminus1, lambda_zminus1pt5, lambda_zminus2, lambda_zminus3))


data = templist

for j in range(len(new_interp_list)):

	intlist = []

	for i in range(len(templist)):

		f_interp_xfix = interpolate.interp1d(zlist, lambdaf[i], kind = 'cubic')

		intlist.append(f_interp_xfix(new_interp_list[j]))	
		
	intlist = np.array(intlist)
	data = np.column_stack((data, intlist))


#f_interp = interpolate.interp2d(10.**X, 10.**Y, 10.**lambdaf)

#lambda_zminus1pt3 = np.log10(f_interp(10.**-1.3, 10.**templist))

plt.plot(templist, data[:,3])
plt.plot(templist, data[:,13])
plt.plot(templist, data[:,33])
plt.plot(templist, lambda_zplus0pt5, '--')
plt.plot(templist, lambda_zsolar, '--')
plt.plot(templist, lambda_zminus0pt5, '--')
plt.plot(templist, lambda_zminus1, '--')
plt.plot(templist, lambda_zminus1pt5, '--')
plt.plot(templist, lambda_zminus2, '--')
plt.plot(templist, lambda_zminus3, '--')


np.savetxt('lambda_cooling_tot.dat', data, fmt='%.4e')

plt.ylim(bottom = -23.6, top = -20.1)
plt.show()



