from pyfd.aux import set_plt_params, plt
import numpy as np
#from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

c = []
error = []
with open('parameter_test.txt') as f:
    counter = 0;
    for line in f:
        c.append([])
        x = repr(line).split()
        error.append(float(x[5][:8]))
        for i in range(1, 5):
            c[counter].append(float(x[i].rstrip(']')))
        counter+=1

ci = np.array(c).T

c_concat = []

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

#for x in ci:
    #c1 = x[0]**2 + x[1]**2
    #c2 = x[2]**2 + x[3]**2 
    #c_concat.append([np.sqrt(c1), np.sqrt(c2)])

x = (ci[0] / np.mean(ci[0])) ** 2 + (ci[1] / np.mean(ci[1]))**2
y = (ci[2] / np.mean(ci[2])) ** 2 + (ci[3] / np.mean(ci[3]))**2

X = []
Y = []
Z = []
colorArray = []
alphaArray = []
it = 0
mult = 9.5
for e in error:
    if e < 0.10:
        X.append(x[it])
        Y.append(y[it])
        Z.append(float(error[it]))
        colorArray.append([1. - mult * error[it], 0., 0.])
    it += 1

X = np.log(X)
Y = np.log(Y)


ax.set_xlabel('breakup parameters')
ax.set_ylabel('coalescence parameters')
ax.set_zlabel('error')

ax.bar3d(X, Y, np.zeros_like(X), 0.5, 0.75, Z, alpha=0.75,
    color=colorArray, edgecolor="none", zorder=10)

#ax.scatter(X, Y, Z, color=colorArray)

N = 20

ct = 0
for tempX in [X, Y]:
    ymin = min(tempX)
    ymax = max(tempX)

    newX = np.linspace(ymin, ymax, N)
    newY = np.zeros_like(newX)

    for yi in tempX:
        it = int((yi - ymin) / (ymax - ymin) * N)
        it = min(it, N - 1)
        newY[it] += 1.0 / len(Y) / 5
    #fig = plt.figure()

    if ct == 1:
        ax.bar(newX, newY, zdir='x', zs=-20)
    else:
        ax.bar(newX, newY, zdir='y', zs=-30)
    ct += 1

plt.show()
