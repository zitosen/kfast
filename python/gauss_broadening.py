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

sigma1 = 0.35 # eV
sigma2 = 0.35
sigma3 = 0.35
xmin = -506.0
xmax = -496.0
grids = list(range(0, 1001))
# total spectra
x1 = [-float(grid)/100-496 for grid in grids]

with open('o1s_all.dat') as file_object:
	rows = file_object.readlines()

mu = []
for row in rows:
	mu.append(float(row))

y1 = []
for x1_grid in x1:
	gau = 0.0
	for row in rows:
		gau = gau + 1 / (sigma1 * np.sqrt(2*math.pi)) *\
np.exp(-0.5*(x1_grid-float(row))**2 / sigma1**2)
	y1.append(gau)
# base spectra
x2 = [-float(grid)/100-496 for grid in grids]

with open('o1s_base.dat') as file_object:
        rows = file_object.readlines()

y2 = []
for x2_grid in x2:
        gau = 0.0
        for row in rows:
                gau = gau + 1 / (sigma2 * np.sqrt(2*math.pi)) *\
np.exp(-0.5*(x2_grid-float(row))**2 / sigma2**2)
        y2.append(gau)
# defect spectra
x3 = [-float(grid)/100-496 for grid in grids]

with open('o1s_defect.dat') as file_object:
        rows = file_object.readlines()

y3 = []
for x3_grid in x3:
        gau = 0.0
        for row in rows:
                gau = gau + 1 / (sigma3 * np.sqrt(2*math.pi)) *\
np.exp(-0.5*(x3_grid-float(row))**2 / sigma3**2)
        y3.append(gau)
#
plt.plot(x1, y1, color='black', linewidth=4.0, label="Total")
plt.plot(x2, y2, color='red', linewidth=4.0, label="Pristine")
plt.plot(x3, y3, color='blue', linewidth=4.0, label="Defect")
#
plt.eventplot(mu, lineoffsets=5, linelengths=10)
plt.xlabel('E (eV)', fontsize=30)
plt.ylabel('Intensity (a.u.)', fontsize=30)
plt.tick_params(axis='both', which='major', labelsize=26)
plt.xlim(xmin, xmax)
plt.ylim(0, 120)

plt.legend(fontsize=26)
plt.show()
