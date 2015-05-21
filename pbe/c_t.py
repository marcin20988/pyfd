from numpy import arange, sum, exp, linspace, sqrt, pi, zeros
from itertools import cycle
from moc import MOCSolution
import matplotlib.pyplot as plt
from fluid import fluid

F = fluid('galinat', 4)
F.C1 = F.C1 * 3.0
F.C2 = F.C2 * 35.0

t = arange(0.0, 20.0 * F.theta, 1e-03)   # residence time is 10 minutes
d0 = 2.89e-03   # let's try 1mm diameter; given in cm
v0 = pi / 6.0 * d0 ** 3     # given in cm^3
s0 = v0 / 12.0
vmax = 1.2 * v0
grids = [80]


pbe_solutions = dict()
for g in grids:
    dv = vmax / g
    v = dv + dv * arange(g)
    Ninit = F.alpha / v0 / s0 / sqrt(2.0 * pi)\
        * exp(- (v - v0) ** 2 / 2 / s0 ** 2)
    pbe_solutions[g] = MOCSolution(
        Ninit, t, dv,
        #Q=Q,
        gamma=F.gamma,
        beta=F.beta,
        theta=F.theta,
        pdf='density'
    )

fig = plt.figure()
ax = fig.gca()
markers = cycle(['o', 's', 'v', '*', '.', ','])

m0_init = sum(Ninit[:] * v[:])
normInit = sum(Ninit)
meanV_init = m0_init / normInit
print "initial m0: ", m0_init
print "initial d: ", (6.0 * meanV_init / pi) ** (1.0 / 3.0)

for n in sorted(pbe_solutions):
    d = (6.0 * pbe_solutions[n].xi / pi) ** (1.0 / 3.0)
    ax.plot(
        d, pbe_solutions[n].N[-1] / (vmax / n), "+",
        marker=next(markers),
        label="MOC with N={0}".format(n))

    meanV = sum(pbe_solutions[n].xi * pbe_solutions[n].N[-1])
    m0 = sum(pbe_solutions[n].xi * pbe_solutions[n].N[-1])
    norm = sum(pbe_solutions[n].N[-1])
    meanV = m0 / norm

    print "m0 / initial m0: ", m0 / m0_init
    print "d=", (6.0 * meanV / pi) ** (1.0 / 3.0)

ax.plot(
    d, Ninit / (vmax / n), "+",
    marker=next(markers),
    label="init N".format(n))
ax.legend(loc='upper right', shadow=True)
ax.set_xlabel('Volume')
ax.set_ylabel('N/N0')

omega = zeros(v.shape[0])
#for i in arange(v.shape[0]):
    #omega[i] = sum(Q(v[i], v) * pbe_solutions[n].N[-1]) * pbe_solutions[n].delta_xi

fig = plt.figure()
ax = fig.gca()
ax.plot(
    v, F.gamma(v), "+",
    marker=next(markers),
    label="breakup".format(n))
ax.legend(loc='upper right', shadow=True)
plt.show()
