from numpy import arange, sum, exp, linspace, sqrt, pi, zeros, genfromtxt
from itertools import cycle
from moc import MOCSolution
import matplotlib.pyplot as plt
from fluid import fluid
import os.path
import numpy as np

pbe_solutions = dict()

F = fluid('coulaloglou', 1, 3)
F.C1 = F.C1 * 1.0
F.C2 = F.C2 * 1.0
F.C3 = F.C3 * 1.0  # e07
F.C4 = F.C4 * 1.0
print F.V

# 20 times residence time should be enough to get convergence
t = arange(0.0, 10, 1e-01)
# data structure: each row contains pair of values Re, d32
# cases are aranged with increasing Re
d0 = 0.5e-03
# convert diameter to volume
v0 = pi / 6.0 * d0 ** 3
# standard deviation
s0 = v0 / 4.0
# this is a breakup dominated case so we don't need large diameters
vmax = 3 * v0
g = 40

dv = vmax / g
v = dv + dv * arange(g)
Ninit = F.alpha / v0 * F.V\
    * 1.0 / s0 / sqrt(2.0 * pi)\
    * exp(- (v - v0) ** 2 / 2 / s0 ** 2)

pbe_solutions[0] = MOCSolution(
    Ninit, t, dv,
    #Q=F.Q,
    #gamma=F.gamma,
    #beta=F.beta,
    pdf='number'
)

m1_init = sum(Ninit[:] * v[:])
normInit = sum(Ninit)
meanV_init = m1_init / normInit
print "initial m1: ", m1_init
print "initial d: ", (6.0 * meanV_init / pi) ** (1.0 / 3.0)

for n in sorted(pbe_solutions):
    fig = plt.figure()
    ax = fig.gca()
    markers = cycle(['o', 's', 'v', '*', '.', ','])

    d = (6.0 * pbe_solutions[n].xi / pi) ** (1.0 / 3.0) * 1e03
    Nplot = pbe_solutions[n].N[-1] / sum(pbe_solutions[n].N[-1])
    ax.plot(
        d, Nplot, "+",
        marker=next(markers),
        label="MOC with N={0}".format(n))

    meanV = sum(pbe_solutions[n].xi * pbe_solutions[n].N[-1])
    m1 = sum(pbe_solutions[n].xi * pbe_solutions[n].N[-1])
    norm = sum(pbe_solutions[n].N[-1])
    meanV = m1 / norm

    print "--------------"
    print 'case ', n
    print "m1 / initial m1: ", m1 / m1_init
    meanD = (6.0 * meanV / pi) ** (1.0 / 3.0)
    print "d=", meanD

    ax.grid()
    ax.set_xlabel('diameter[mm]')
    ax.set_ylabel('frequency')

ax.plot(
    d, Ninit / sum(Ninit), "+",
    marker=next(markers),
    label="Ninit".format(n))
ax.legend(loc='best')

coalescence = zeros(v.shape[0])
#breakup = F.gamma(pbe_solutions[0].xi) * pbe_solutions[0].N[-1]
breakup = F.gamma(pbe_solutions[0].xi)
N = pbe_solutions[0].N[-1]
for i in arange(v.shape[0]):
    coalescence[i] =\
        sum(F.Q(pbe_solutions[0].xi[i], pbe_solutions[0].xi) * N)

fig = plt.figure()
ax = fig.gca()

ax.plot(
    d, coalescence, "+",
    marker=next(markers),
    label="coalescence frequency".format(n))
ax.plot(
    d, breakup, "+",
    marker=next(markers),
    label="breakup frequency".format(n))
ax.legend(loc="best")
#ax.set_xlim(0, 0.6)
ax.set_yscale('log')

plt.show()