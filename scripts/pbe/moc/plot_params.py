import numpy as np
from mpl_toolkits.mplot3d import Axes3D
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

x = (ci[0] / np.mean(ci[0])) + (ci[1] / np.mean(ci[1]))
y = (ci[2] / np.mean(ci[2])) + (ci[3] / np.mean(ci[3]))

ax.bar3d(x, y, np.zeros_like(x), 1, 1, np.array(error), cmap='coolwarm', alpha=0.1)
ax.set_xscale('log')
ax.set_yscale('log')

plt.show()
