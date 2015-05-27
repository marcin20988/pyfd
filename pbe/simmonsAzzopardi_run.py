from numpy import arange, sum, exp, linspace, sqrt, pi, zeros, genfromtxt
from itertools import cycle
from moc import MOCSolution
import matplotlib.pyplot as plt
from fluid import fluid
import os.path
import numpy as np

pbe_solutions = dict()

F = fluid('simmonsAzzopardi')
F.C1 = F.C1 * 1.0
F.C2 = F.C2 * 1.0
F.C3 = F.C3 * 1.0
F.C4 = F.C4 * 1.0
print F.V

t = arange(0.0, 1.0, 1e-02)
# data structure: each row contains pair of values Re, d32
# cases are aranged with increasing Re
d0 = 500e-06
# convert diameter to volume
v0 = pi / 6.0 * d0 ** 3
# standard deviation
s0 = v0 / 4.0
# experimentally measured diameters are between 200-1000 microns
vmax = 8 * v0
g = 180

dv = vmax / g
v = dv + dv * arange(g)
Ninit = F.alpha / v0 * F.V / (max(v) - min(v))\
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


#breakup = F.gamma(v)
#for i in arange(v.shape[0]):
    #coalescence[i] =\
        #sum(F.Q(v[i], v) * pbe_solutions[0].delta_xi * pbe_solutions[0].xi)

coalescence = zeros(v.shape[0])
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
ax.set_yscale('log')

plt.show()
