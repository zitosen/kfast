"""
=================
Cell Performance
=================

Create multiple scatter plots with different
star symbols to show cell performance

by zito, 2020/07/17
Email: shenzt@vip.henu.edu.cn

modified by zito, 2020/07/24

"""

import matplotlib.pyplot as plt

filename = input("Give me your filename: \n")
cell_performance = input("Which cell performance do you want to see?\n Voc: 1 | Jsc: 2 | PCE: 3 | FF: 4 \n")
cell_performance = int(cell_performance)

with open(filename, encoding='GBK') as file_object:
    lines = file_object.readlines()

y1 = []
y2 = []
y3 = []
y4 = []
y5 = []
y6 = []
y1_floats = []
y2_floats = []
y3_floats = []
y4_floats = []
y5_floats = []
y6_floats = []
x1 = []
x2 = []
x3 = []
x4 = []
x5 = []
x6 = []
for line in lines:
    if line[2] == '1' and line[1] == '-':
        line_split = line[:].split()
        y1.append(line_split[cell_performance])
        if line[0] == '1':
            x1.append(1)
        elif line[0] == '2':
            x1.append(2)
        elif line[0] == '3':
            x1.append(3)
        elif line[0] == '4':
            x1.append(4)
        elif line[0] == '5':
            x1.append(5)
        elif line[0] == '6':
            x1.append(6)
        elif line[0] == '7':
            x1.append(7)
        elif line[0] == '8':
            x1.append(8)
        elif line[0] == '9':
            x1.append(9)
    elif line[2] == '2' and line[1] == '-':
        line_split = line[:].split()
        y2.append(line_split[cell_performance])
        if line[0] == '1':
            x2.append(1)
        elif line[0] == '2':
            x2.append(2)
        elif line[0] == '3':
            x2.append(3)
        elif line[0] == '4':
            x2.append(4)
        elif line[0] == '5':
            x2.append(5)
        elif line[0] == '6':
            x2.append(6)
        elif line[0] == '7':
            x2.append(7)
        elif line[0] == '8':
            x2.append(8)
        elif line[0] == '9':
            x2.append(9)
    elif line[2] == '3' and line[1] == '-':
        line_split = line[:].split()
        y3.append(line_split[cell_performance])
        if line[0] == '1':
            x3.append(1)
        elif line[0] == '2':
            x3.append(2)
        elif line[0] == '3':
            x3.append(3)
        elif line[0] == '4':
            x3.append(4)
        elif line[0] == '5':
            x3.append(5)
        elif line[0] == '6':
            x3.append(6)
        elif line[0] == '7':
            x3.append(7)
        elif line[0] == '8':
            x3.append(8)
        elif line[0] == '9':
            x3.append(9)
    elif line[2] == '4' and line[1] == '-':
        line_split = line[:].split()
        y4.append(line_split[cell_performance])
        if line[0] == '1':
            x4.append(1)
        elif line[0] == '2':
            x4.append(2)
        elif line[0] == '3':
            x4.append(3)
        elif line[0] == '4':
            x4.append(4)
        elif line[0] == '5':
            x4.append(5)
        elif line[0] == '6':
            x4.append(6)
        elif line[0] == '7':
            x4.append(7)
        elif line[0] == '8':
            x4.append(8)
        elif line[0] == '9':
            x4.append(9)
    elif line[2] == '5' and line[1] == '-':
        line_split = line[:].split()
        y5.append(line_split[cell_performance])
        if line[0] == '1':
            x5.append(1)
        elif line[0] == '2':
            x5.append(2)
        elif line[0] == '3':
            x5.append(3)
        elif line[0] == '4':
            x5.append(4)
        elif line[0] == '5':
            x5.append(5)
        elif line[0] == '6':
            x5.append(6)
        elif line[0] == '7':
            x5.append(7)
        elif line[0] == '8':
            x5.append(8)
        elif line[0] == '9':
            x5.append(9)
    elif line[2] == '6' and line[1] == '-':
        line_split = line[:].split()
        y6.append(line_split[cell_performance])
        if line[0] == '1':
            x6.append(1)
        elif line[0] == '2':
            x6.append(2)
        elif line[0] == '3':
            x6.append(3)
        elif line[0] == '4':
            x6.append(4)
        elif line[0] == '5':
            x6.append(5)
        elif line[0] == '6':
            x6.append(6)
        elif line[0] == '7':
            x6.append(7)
        elif line[0] == '8':
            x6.append(8)
        elif line[0] == '9':
            x6.append(9)

for line in lines:
    if line[3] == '1' and line[2] == '-':
        line_split = line[:].split()
        y1.append(line_split[cell_performance])
        if line[1] == '0':
            x1.append(10)
        elif line[1] == '1':
            x1.append(11)
        elif line[1] == '2':
            x1.append(12)
    elif line[3] == '2' and line[2] == '-':
        line_split = line[:].split()
        y2.append(line_split[cell_performance])
        if line[1] == '0':
            x2.append(10)
        elif line[1] == '1':
            x2.append(11)
        elif line[1] == '2':
            x2.append(12)
    elif line[3] == '3' and line[2] == '-':
        line_split = line[:].split()
        y3.append(line_split[cell_performance])
        if line[1] == '0':
            x3.append(10)
        elif line[1] == '1':
            x3.append(11)
        elif line[1] == '2':
            x3.append(12)
    elif line[3] == '4' and line[2] == '-':
        line_split = line[:].split()
        y4.append(line_split[cell_performance])
        if line[1] == '0':
            x4.append(10)
        elif line[1] == '1':
            x4.append(11)
        elif line[1] == '2':
            x4.append(12)
    elif line[3] == '5' and line[2] == '-':
        line_split = line[:].split()
        y5.append(line_split[cell_performance])
        if line[1] == '0':
            x5.append(10)
        elif line[1] == '1':
            x5.append(11)
        elif line[1] == '2':
            x5.append(12)
    elif line[3] == '6' and line[2] == '-':
        line_split = line[:].split()
        y6.append(line_split[cell_performance])
        if line[1] == '0':
            x6.append(10)
        elif line[1] == '1':
            x6.append(11)
        elif line[1] == '2':
            x6.append(12)

for y1_string in y1:
    y1_floats.append(float(y1_string))
for y2_string in y2:
    y2_floats.append(float(y2_string))
for y3_string in y3:
    y3_floats.append(float(y3_string))
for y4_string in y4:
    y4_floats.append(float(y4_string))
for y5_string in y5:
    y5_floats.append(float(y5_string))
for y6_string in y6:
    y6_floats.append(float(y6_string))

z1 = []
z2 = []
z3 = []
z4 = []
z5 = []
z6 = []
for y1_float in y1_floats:
    z1.append(y1_float * 100)
for y2_float in y2_floats:
    z2.append(y2_float * 100)
for y3_float in y3_floats:
    z3.append(y3_float * 100)
for y4_float in y4_floats:
    z4.append(y4_float * 100)
for y5_float in y5_floats:
    z5.append(y5_float * 100)
for y6_float in y6_floats:
    z6.append(y6_float * 1002)

fig, axs = plt.subplots(3, 2)
if cell_performance == 1:
    fig.suptitle(filename + ': Voc', fontsize=16)
    axs[0, 0].set(ylim=(0.8, 1.2))
    axs[0, 0].scatter(x1, y1_floats, s=100, c=z1, marker=">")
    axs[0, 0].grid(True)
    axs[0, 1].set(ylim=(0.8, 1.2))
    axs[0, 1].scatter(x4, y4_floats, s=100, c=z4, marker=">")
    axs[0, 1].grid(True)
    axs[1, 0].set(ylim=(0.8, 1.2))
    axs[1, 0].scatter(x2, y2_floats, s=100, c=z2, marker=(5, 2))
    axs[1, 0].grid(True)
    axs[1, 1].set(ylim=(0.8, 1.2))
    axs[1, 1].scatter(x5, y5_floats, s=100, c=z5, marker=(5, 2))
    axs[1, 1].grid(True)
    axs[2, 0].set(ylim=(0.8, 1.2))
    axs[2, 0].scatter(x3, y3_floats, s=100, c=z3, marker=(5, 1))
    axs[2, 0].grid(True)
    axs[2, 1].set(ylim=(0.8, 1.2))
    axs[2, 1].scatter(x6, y6_floats, s=100, c=z6, marker=(5, 1))
    axs[2, 1].grid(True)
elif cell_performance == 2:
    fig.suptitle(filename + ': Jsc', fontsize=16)
    axs[0, 0].set(ylim=(20, 30))
    axs[0, 0].scatter(x1, y1_floats, s=100, c=z1, marker=">")
    axs[0, 0].grid(True)
    axs[0, 1].set(ylim=(20, 30))
    axs[0, 1].scatter(x4, y4_floats, s=100, c=z4, marker=">")
    axs[0, 1].grid(True)
    axs[1, 0].set(ylim=(20, 30))
    axs[1, 0].scatter(x2, y2_floats, s=100, c=z2, marker=(5, 2))
    axs[1, 0].grid(True)
    axs[1, 1].set(ylim=(20, 30))
    axs[1, 1].scatter(x5, y5_floats, s=100, c=z5, marker=(5, 2))
    axs[1, 1].grid(True)
    axs[2, 0].set(ylim=(20, 30))
    axs[2, 0].scatter(x3, y3_floats, s=100, c=z3, marker=(5, 1))
    axs[2, 0].grid(True)
    axs[2, 1].set(ylim=(20, 30))
    axs[2, 1].scatter(x6, y6_floats, s=100, c=z6, marker=(5, 1))
    axs[2, 1].grid(True)
elif cell_performance == 3:
    fig.suptitle(filename + ': PCE', fontsize=16)
    axs[0, 0].set(ylim=(0.1, 0.2))
    axs[0, 0].scatter(x1, y1_floats, s=100, c=z1, marker=">")
    axs[0, 0].grid(True)
    axs[0, 1].set(ylim=(0.1, 0.2))
    axs[0, 1].scatter(x4, y4_floats, s=100, c=z4, marker=">")
    axs[0, 1].grid(True)
    axs[1, 0].set(ylim=(0.1, 0.2))
    axs[1, 0].scatter(x2, y2_floats, s=100, c=z2, marker=(5, 2))
    axs[1, 0].grid(True)
    axs[1, 1].set(ylim=(0.1, 0.2))
    axs[1, 1].scatter(x5, y5_floats, s=100, c=z5, marker=(5, 2))
    axs[1, 1].grid(True)
    axs[2, 0].set(ylim=(0.1, 0.2))
    axs[2, 0].scatter(x3, y3_floats, s=100, c=z3, marker=(5, 1))
    axs[2, 0].grid(True)
    axs[2, 1].set(ylim=(0.1, 0.2))
    axs[2, 1].scatter(x6, y6_floats, s=100, c=z6, marker=(5, 1))
    axs[2, 1].grid(True)
if cell_performance == 4:
    fig.suptitle(filename + ': FF', fontsize=16)
    axs[0, 0].set(ylim=(0.5, 0.8))
    axs[0, 0].scatter(x1, y1_floats, s=100, c=z1, marker=">")
    axs[0, 0].grid(True)
    axs[0, 1].set(ylim=(0.5, 0.8))
    axs[0, 1].scatter(x4, y4_floats, s=100, c=z4, marker=">")
    axs[0, 1].grid(True)
    axs[1, 0].set(ylim=(0.5, 0.8))
    axs[1, 0].scatter(x2, y2_floats, s=100, c=z2, marker=(5, 2))
    axs[1, 0].grid(True)
    axs[1, 1].set(ylim=(0.5, 0.8))
    axs[1, 1].scatter(x5, y5_floats, s=100, c=z5, marker=(5, 2))
    axs[1, 1].grid(True)
    axs[2, 0].set(ylim=(0.5, 0.8))
    axs[2, 0].scatter(x3, y3_floats, s=100, c=z3, marker=(5, 1))
    axs[2, 0].grid(True)
    axs[2, 1].set(ylim=(0.5, 0.8))
    axs[2, 1].scatter(x6, y6_floats, s=100, c=z6, marker=(5, 1))
    axs[2, 1].grid(True)

plt.show()
