"""
=================
Plot DOS
=================

by zito, 2020/10/26
Email: shenzt@vip.henu.edu.cn

"""

import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator
from scipy.interpolate import make_interp_spline


with open('PDOS_SUM.dat', encoding='GBK') as file_object:
    rows = file_object.readlines()

x = []
dos_tot = []
dos_px = []
dos_py = []
dos_pz = []
for row in rows:
    row_split = row[:].split()
    if row_split[0] != '#Energy':
        x.append(float(row_split[0]))
        dos_tot.append(float(row_split[10]))
        dos_py.append(float(row_split[2]))
        dos_pz.append(float(row_split[3]))
        dos_px.append(float(row_split[4]))

plt.plot(x, dos_tot, c='red', linewidth=2.5)
plt.plot(x, dos_px, c='blue', linewidth=2.5)
plt.plot(x, dos_py, c='green', linewidth=2.5)
plt.plot(x, dos_pz, c='yellow', linewidth=2.5)
plt.tick_params(axis='both', which='major', labelsize=24)
x_major_locator = MultipleLocator(1)
y_major_locator = MultipleLocator(10)
plt.xlabel('E (eV)', fontsize=30)
plt.ylabel('DOS', fontsize=30)
ax = plt.gca()
ax.plot([0, 0], [0, 50], c='black', dashes=[4, 2], linewidth=2)

ax.xaxis.set_major_locator(x_major_locator)
ax.yaxis.set_major_locator(y_major_locator)
plt.xlim(-2, 2)
plt.ylim(0, 50)

plt.show()