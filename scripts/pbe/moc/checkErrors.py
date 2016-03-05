import pickle
from pyfd.pbe.moc import CTSolution, AngeliSolution
from numpy import exp, array, pi

sets = ['angeli', 'ct']
#sets = ['ct']
multipliers = {'ct': [0.1, 0.1, 1., 1e12],
        'angeli': [1., 1., 1., 1e10]}

results = dict()

with open('data/ct_optimization_results.pickle', 'rb') as f:
    results['ct'] = pickle.load(f)

with open('data/angeli_optimization_results.pickle', 'rb') as f:
    results['angeli'] = pickle.load(f)


C1 = []
for key in sets:
    for data in results[key]:

        c = array(data['best_fit'])
        c[:] = exp(c[:])
        for i in range(4):
            c[i] *= multipliers[key][i]
        print('================= ', key)
        print(c)

        s = data['setup']
        v0 = pi / 6 * s['d32']**3

        if key == 'ct':
            nstar_to_eps = 0.407
            pbe_solutions = \
                CTSolution(
                    M=40, v0=v0, Nstar=s['Nstar'], phi=s['phi'],
                    model_parameters=c)
            c[0] = c[0] * nstar_to_eps ** (1./3.)
            c[1] = c[1] / nstar_to_eps ** (2./3.)
            c[2] = c[2] * nstar_to_eps ** (1./3.)
            c[3] = c[3] / nstar_to_eps
        elif key == 'angeli':
            pbe_solutions = AngeliSolution(
                M=40, v0=v0, U=s['U'], phi=s['phi'], theta=s['theta'],
                model_parameters=c)

        print (pbe_solutions.d32, s['d32'], (pbe_solutions.d32 - s['d32']) / s['d32'])
