from numpy import arange, sum, exp, linspace, sqrt, pi, zeros, genfromtxt
from itertools import cycle
from moc import MOCSolution
import matplotlib.pyplot as plt
from fluid import fluid
import os.path
import numpy as np
from scipy.optimize import minimize, fmin, fmin_powell, differential_evolution
from sys import argv


def case_error(C):
    pbe_solutions = dict()
    F = fluid(argv[1], (int)(argv[2]), (int)(argv[3]))

    F.C1 = C[0]
    F.C2 = C[1]
    F.C3 = C[2]
    F.C4 = C[3]
    print '--------------'
    print 'constants: ', F.C1, F.C2, F.C3, F.C4

    Er = 0.0
    for i in np.arange(F.d0.shape[0]):
      t = F.timeRange
      v0 = F.v0[i]
      s0 = F.s0[i]
      vmax = F.vMax[i]
      g = F.numberOfClasses
      print F.d0
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
      m1 = sum(N[:] * v[:])
      m1Init = sum(Ninit[:] * v[:])
      norm = sum(N)
      meanV = m1 / norm
      dMean = (6.0 * meanV / pi) ** (1.0 / 3.0)
      if m1 / m1Init < 0.5 or m1 / m1Init > 1.5:
	print 'mass not conserved'
	return 1e12

      print 'mass conservations: ', m1 / m1Init
      print 'calculated, experimental d: ', dMean, F.expectedD
      Er += sqrt((dMean - F.expectedD) ** 2) / F.expectedD
    return Er

#C = fmin(case_error, C0, full_output=True, maxiter=2)
C0 = np.array([1e-02, 0.005, 1.0e-12, 1.0e+014])
C = case_error(C0)
#bounds = [(0, 1), (0, 1), (1e-08, 1e-17), (1e08, 1e18)]
#C = differential_evolution(case_error, bounds, maxiter=50, popsize=10)

print C
