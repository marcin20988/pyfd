from numpy import genfromtxt, abs, array, pi, exp
import numpy as np
from scipy.optimize import minimize, differential_evolution
from pyfd.pbe.moc import KarabelasSolution, KarabelasSolutionHighViscosity
import time
import pickle


class angeli_experiment:
    def __init__(self, U, d32, theta=1200.):
        self.theta = theta
        self.U = U
        self.d32 = d32

def error_function(C, experiment):
    v0s = array([0.5, 1.5]) * pi / 6 * experiment.d32**3

    mp = array(C)
    mp[:]=exp(0.1 * mp[:])
    mp[3] *= 1e12
    #mp = np.append(mp, 7.2e12)
    mp[2] *= 1e-02
    mp[1] *= 0.1
    mp[0] *= 1.

    pbe_solutions = [
        KarabelasSolution(
            M=20, v0=v0, U=experiment.U, theta=experiment.theta,
            model_parameters=mp)
        for v0 in v0s]
    error = sum(
        [abs(s.d32 - experiment.d32) / experiment.d32
            for s in pbe_solutions])
    print('constants={0}, error={1}'.format(C, error))
    print('d=[{0}, {1}], expd={2}'.format(pbe_solutions[0].d32, pbe_solutions[1].d32, experiment.d32))
    return error

experiments = []
Us = [2.98, 2.57, 2.22, 1.84]
ds = [407.0, 455.0, 794.0, 1066.0]
for i in range(4):
    experiments.append(angeli_experiment(Us[i], ds[i] * 1e-06))

# original C&T parameters where multiplying (Nstar**3 * D**2)
# and epsilon = 0.407 * NStar**3 * D**2
# so
s=0.407
c0 = [0.4 * s**(-1./3.), 0.08 / s**(-2./3.), 2.8 * s**(-1./3.), 1.83 * s]  # CT original constants

c0 = [1., 1., 1., 1.]  # CT original constants
results = []
for e in experiments:
    res = dict()

    Copt = minimize(
        lambda c: error_function(c, e), c0,
        method='L-BFGS-B',
        options={'disp': False, 'ftol': 0.001, 'maxiter': 30})
    #Copt = differential_evolution(
        #lambda c: error_function(c, e), 
        #bounds=[(0, 1e06)] * 4,
        #maxiter=1000,
        #polish=True
        #)
    res['best_fit'] = Copt.x
    res['setup'] = {'U': e.U, 'd32': e.d32, 'theta': e.theta}
    results.append(res)

with open('data/karabelas_optimization_results.pickle', 'wb') as f:
    pickle.dump(results, f)
