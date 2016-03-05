import pickle
from pyfd.pbe.moc import CTSolution, AngeliSolution, KarabelasSolution, KarabelasSolutionHighViscosity
from numpy import exp, array, pi
import matplotlib.pyplot as plt

sets = ['angeli', 'ct', 'karabelas', 'karabelas_high']
#sets = ['karabelas']
multipliers = {'ct': [0.1, 0.1, 1., 1e12],
        'angeli': [1., 1., 1., 1e10],
        'karabelas': [0.1, 0.1, 1., 1e12],
        'karabelas_high': [0.1, 0.001, 1., 1e12]}

results = dict()

with open('data/ct_optimization_results.pickle', 'rb') as f:
    results['ct'] = pickle.load(f)

with open('data/angeli_optimization_results.pickle', 'rb') as f:
    results['angeli'] = pickle.load(f)

with open('data/karabelas_optimization_results.pickle', 'rb') as f:
    results['karabelas'] = pickle.load(f)

with open('data/karabelas_high_optimization_results.pickle', 'rb') as f:
    results['karabelas_high'] = pickle.load(f)

res = {'c1': [], 'c2': [], 'c3': [], 'c4': [], 'Re': [], 'St': [], 'Ca': [], 'We': []}
for key in sets:
    for data in results[key]:
        x = dict()

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

            Re = pbe_solutions.Nstar * pbe_solutions.D ** 2 / pbe_solutions.contProperties['mu']\
                * pbe_solutions.contProperties['rho']

        elif key == 'angeli':
            pbe_solutions = AngeliSolution(
                M=40, v0=v0, U=s['U'], phi=s['phi'], theta=s['theta'],
                model_parameters=c)
            Re = pbe_solutions.Re

        elif key =='karabelas':
            pbe_solutions = \
            KarabelasSolution(
                M=40, v0=v0, U=s['U'], theta=s['theta'],
                model_parameters=c)
            Re = pbe_solutions.Re

        elif key =='karabelas_high':
            pbe_solutions = \
            KarabelasSolutionHighViscosity(
                M=40, v0=v0, U=s['U'], theta=s['theta'],
                model_parameters=c)
            Re = pbe_solutions.Re

        cont = pbe_solutions.contProperties
        disp = pbe_solutions.dispProperties

        St = 2. * cont['rho'] / (2. * disp['rho'] + cont['rho'])\
            * 2. / 9. * pbe_solutions.d32 ** 2. / pbe_solutions.D ** 2 * Re
        c[2] = c[2] * pbe_solutions.Vt

        res['c1'].append(c[0])
        res['c2'].append(c[1])
        res['c3'].append(c[2])
        res['c4'].append(c[3])
        res['Re'].append(Re)
        res['St'].append(St)

        print (pbe_solutions.d32, s['d32'], (pbe_solutions.d32 - s['d32']) / s['d32'])
