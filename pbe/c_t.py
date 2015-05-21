from numpy import arange, sum, exp, linspace, sqrt, pi, zeros
from itertools import cycle
from moc import MOCSolution
import matplotlib.pyplot as plt
from fluid import fluid

F = fluid('galinat')

print F.thetas
print F.epsilons

# goal:
# to reproduce figure 2 from Coulaloglou and Tavlarides paper

t = arange(0.0, 1e-02, 1e-03)   # residence time is 10 minutes
d0 = 0.025e-02   # let's try 1mm diameter; given in cm
v0 = pi / 6.0 * d0 ** 3     # given in cm^3
s0 = v0 / 3.0
vmax = 5.0 * v0
grids = [40]

# impeller speed and diameter
D = 0.1	 # 10.0     # 10 cm; given in cm
Nstar = 310.0 / 60.0  # 310 rpm; given in 1/s

# volume fraction
Phi = 0.1
rhod = 820  # 0.82  # TODO: check this; given in g / cm^3
sigma = 42.0e-03  # TODO: check this; given in dynes / cm [= 0.001 N / m]
muc = 1e-03  # TODO: check this; given in g / cm / s [= 1 poise = 0.1 Pa s]
rhoc = 1000  # TODO: check this; in g/ cm ^ 3
theta = 10.0 * 60


#------------------------------------------------------------------
# function EXACTLY as given in the Coulaloglou and Tavlarides paper:

# droplet daughter distribution:
def beta(xi2, xi1):
    return 2.0 * 2.4 / xi1 * exp(- 4.5 * (2.0 * xi2 - xi1) ** 2 / xi1 ** 2)


# breakup rate:
def gamma(xi):
    #C1 = 0.336
    #C2 = 0.106
    C1 = 0.4
    C2 = 0.08
    C = C1 * xi ** (-2.0 / 9.0) * D ** (2.0 / 3.0) * Nstar / (1.0 + Phi)
    exp_argument = - C2 * sigma * (1.0 + Phi) ** 2\
        / (rhod * xi ** (5.0 / 9.0) * D ** (4.0 / 3.0) * Nstar ** 2)
    return C * exp(exp_argument)


# coalescence rate:
def Q(xi1, xi2):
    C3 = 2.8e-06
    C4 = 1.83e09

    #C3 = 2.32e-06
    #C4 = 1.2e09

    #C3 = 2.32
    #C4 = 1.2e09

    dRatio = xi1 ** (1.0 / 3.0) * xi2 ** (1.0 / 3.0)\
        / (xi1 ** (1.0 / 3.0) + xi2 ** (1.0 / 3.0))
    dRatio = dRatio ** 4

    exp_argument = - C4 * muc * rhoc * D ** 2 / sigma ** 2
    exp_argument *= Nstar ** 3 / (1.0 + Phi) ** 3 * dRatio

    C = C3 * (xi1 ** (2.0 / 3.0) + xi2 ** (2.0 / 3.0)) * D ** (2.0 / 3.0)\
        * (xi1 ** (2.0 / 9.0) + xi2 ** (2.0 / 9.0)) ** 0.5\
        * Nstar / (1.0 + Phi)
    return exp(exp_argument) * C
    #return 1e-09
    #return xi1 ** (2.0 / 3.0)


pbe_solutions = dict()
for g in grids:
    dv = vmax / g
    v = dv + dv * arange(g)
    #Ninit = (N0 / v0) * (v / v0) * exp(-v / v0) * dv
    #Ninit = Phi / v0 / s0 / sqrt(2.0 * pi) * exp(- (v - v0) ** 2 / 2 / s0)
    Ninit = Phi / v0 / s0 / sqrt(2.0 * pi) * exp(- (v - v0) ** 2 / 2 / s0 ** 2)
    #Ninit /= sum(Ninit)
    pbe_solutions[g] = MOCSolution(
        Ninit, t, dv,
        Q=Q,
        gamma=gamma,
        beta=beta,
        theta=theta,
	pdf='density'
    )

fig = plt.figure()
ax = fig.gca()
markers = cycle(['o', 's', 'v', '*', '.', ','])
#v = linspace(0, vmax, 100000)

m0_init = sum(Ninit[:] * v[:])
normInit = sum(Ninit)
meanV_init = m0_init / normInit
print "initial m0: ", m0_init
print "initial d: ", (6.0 * meanV_init / pi) ** (1.0 / 3.0)

for n in sorted(pbe_solutions):
    ax.plot(
        pbe_solutions[n].xi, pbe_solutions[n].N[-1] / (vmax / n), "+",
        #pbe_solutions[n].xi, pbe_solutions[n].N[-1], "+",
        #pbe_solutions[n].xi, pbe_solutions[n].N[-1], "+",
        marker=next(markers),
        label="MOC with N={0}".format(n))

    meanV = sum(pbe_solutions[n].xi * pbe_solutions[n].N[-1])
    m0 = sum(pbe_solutions[n].xi * pbe_solutions[n].N[-1])
    norm = sum(pbe_solutions[n].N[-1])
    meanV = m0 / norm

    print "m0 / initial m0: ", m0 / m0_init
    print "d=", (6.0 * meanV / pi) ** (1.0 / 3.0)

ax.plot(
    v, Ninit / (vmax / n), "+",
    marker=next(markers),
    label="init N".format(n))
ax.legend(loc='upper right', shadow=True)
ax.set_xlabel('Volume')
ax.set_ylabel('N/N0')
#plt.show()

omega = zeros(v.shape[0])
for i in arange(v.shape[0]):
    #for j in arange(v.shape[0]):
    omega[i] = sum(Q(v[i], v) * pbe_solutions[n].N[-1]) * pbe_solutions[n].delta_xi
#omega /= norm

fig = plt.figure()
ax = fig.gca()
#ax.set_xlim(0, 6.5e-05)
#ax.set_ylim(0, 0.8)
ax.plot(
    v, gamma(v), "+",
    marker=next(markers),
    label="breakup".format(n))
ax.plot(
    v, omega, "+",
    marker=next(markers),
    label="coalescence".format(n))
ax.legend(loc='upper right', shadow=True)
plt.show()
