"""
============================================
 spectra broadening using Gaussian function:

g(x)=1/(sigma*sqrt(2pi))*exp(-0.5*(x-mu)**2 / simga**2)

============================================

by zito, 2020/11/4
Email: shenzt@vip.henu.edu.cn

"""

import numpy as np
import matplotlib.pyplot as plt
import math

sigma = 0.2 # eV
grids = list(range(0, 201))
x = [-float(grid)/100-500 for grid in grids]

with open('O_1s.dat') as file_object:
	rows = file_object.readlines()

mu = []
for row in rows:
	mu.append(float(row))

y = []
#gau_tests = []
for x_grid in x:
	gau = 0.0
	for row in rows:
		gau = gau + 1 / (sigma * np.sqrt(2*math.pi)) *\
np.exp(-0.5*(x_grid-float(row))**2 / sigma**2)
	y.append(gau)
#	gau_test = 1 / (sigma * np.sqrt(2*math.pi)) *\
#np.exp(-0.5*(x_grid+501.0)**2 / sigma**2)
#	gau_tests.append(gau_test)

plt.plot(x, y, color='red', linewidth=2.0, label="O1s")
#plt.plot(x,gau_tests)
plt.eventplot(mu, lineoffsets=5, linelengths=10)
plt.xlabel('E (eV)', fontsize=24)
plt.ylabel('Intensity (a.u.)', fontsize=24)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.xlim(-502, -500)
plt.ylim(0, 120)

plt.legend(fontsize=20)
plt.show()
