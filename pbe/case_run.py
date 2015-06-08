from numpy import arange, sum, exp, linspace, sqrt, pi, zeros, genfromtxt
from itertools import cycle
from moc import MOCSolution
import matplotlib.pyplot as plt
from fluid import fluid
import os.path
import numpy as np
from scipy.optimize import minimize, fmin
from sys import argv, exit


def case_error(C):
    pbe_solutions = dict()
    F = fluid(argv[1], (int)(argv[2]), (int)(argv[3]))

    F.C1 = C[0]
    F.C2 = C[1]
    F.C3 = C[2]
    F.C4 = C[3]

    t = F.timeRange
    t = arange(0, 1.0e-02, 1e-02)
    Er = 0.0
    for i in np.arange(2):
      v0 = F.v0[i]
      s0 = F.s0[i]
      vmax = F.vMax[i]
      g = F.numberOfClasses
      print "initial d and v: ", i, F.d0[i], v0

      dv = vmax / g
      v = dv + dv * arange(g)
      Ninit = F.alpha / v0 * F.V\
	  * 1.0 / s0 / sqrt(2.0 * pi)\
	  * exp(- (v - v0) ** 2 / 2 / s0 ** 2)

      pbe_solutions[0] = MOCSolution(
	  Ninit, t, dv,
	  Q=F.Q,
	  gamma=F.gamma,
	  beta=F.beta,
	  pdf='number'
      )

      N = pbe_solutions[0].N[-1]
      d = (6.0 * pbe_solutions[0].xi / pi) ** (1.0 / 3.0) * 1e03

      m1 = sum(N[:] * v[:])
      norm = sum(N)
      meanV = m1 / norm
      dMean = (6.0 * meanV / pi) ** (1.0 / 3.0)

      fig = plt.figure()
      ax = fig.gca()
      ax.plot(
	  d, N, "+",
	  label="result")
      ax.plot(
	  d, Ninit, "+",
	  label="Ninit")

      coalescence = zeros(v.shape[0])
      breakup = F.gamma(pbe_solutions[0].xi)
      N = pbe_solutions[0].N[-1]

      for j in arange(v.shape[0]):
	  coalescence[j] =\
	      sum(F.Q(pbe_solutions[0].xi[j], pbe_solutions[0].xi) * N)

      fig2 = plt.figure()
      ax = fig2.gca()
      ax.plot(
	  d, coalescence, "+",
	  label="coalescence frequency")
      ax.plot(
	  d, breakup, "+",
	  label="breakup frequency")
      ax.legend(loc="best")
      ax.set_yscale('log')

      print '--------------'
      print F.C1, F.C2
      print dMean, F.expectedD
      Er += sqrt((dMean - F.expectedD) ** 2) / F.expectedD
    return Er

# just to get initial values for a case
#C0 = np.array([0.32, 0.45, 1.15e-09, 3.2e15])
#C0 = np.array([0.0407, 0.3527, 6.16e-09, 1.94e16])
#C0 = np.array([0.05, 0.005, 2.0e-07, 1.0e+015])

#C0 = np.array([0.01, 0.005, 1.0e-13, 1.0e+014])
C0 = np.array([0.31, 0.45, 1.0e-08, 1.0e+014])
#F = fluid('simmonsAzzopardi', 1, 3)
#C0 = np.array([5.47530498e-01, 1.60420466e-01, 6.91294878e-09,
              #9.46910002e+17])
#C0
C = case_error(C0)
plt.show()
#for i in arange(1):
    #F = fluid('galinat', i, 0)
    #print F.Re, F.St, F.Ca

#print C
