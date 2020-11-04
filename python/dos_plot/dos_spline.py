"""
=================
Plot DOS
=================

by zito, 2020/10/26
Email: shenzt@vip.henu.edu.cn

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator
from scipy.interpolate import make_interp_spline


with open('PDOS_SUM.dat', encoding='GBK') as file_object:
    rows = file_object.readlines()

x = []
y = []
for row in rows:
    row_split = row[:].split()
    if row_split[0] != '#Energy':
        x.append(float(row_split[0]))
        y.append(float(row_split[10]))

x = np.array(x)
y = np.array(y)
xp = np.linspace(x.min(), x.max(), 1000)
y_smooth = make_interp_spline(x, y)(xp)

plt.plot(xp, y_smooth, c='red')
plt.tick_params(axis='both', which='major', labelsize=16)
x_major_locator = MultipleLocator(1)
y_major_locator = MultipleLocator(10)
plt.xlabel('E (eV)', fontsize=20)
plt.ylabel('DOS', fontsize=20)
ax = plt.gca()
ax.plot([0, 0], [0, 50], c='black', dashes=[4, 2])

ax.xaxis.set_major_locator(x_major_locator)
ax.yaxis.set_major_locator(y_major_locator)
plt.xlim(-2, 2)
plt.ylim(0, 50)

plt.show()