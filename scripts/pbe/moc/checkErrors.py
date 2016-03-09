import pickle
from pyfd.pbe.moc import CTSolution, AngeliSolution, KarabelasSolution, KarabelasSolutionHighViscosity, SASolution
from numpy import exp, array, pi, sqrt

sets = ['sa', 'angeli', 'ct', 'karabelas']
#sets = ['angeli']
multipliers = {'ct': [0.1, 0.1, 1., 1e12],
        'angeli': [0.01, 0.01, 1000., 1e12],
        'karabelas': [1., 0.1, 0.01, 1e12],
        'karabelas_high': [0.1, 0.1, 100., 1e12],
        'sa': [0.1, 0.1, 1., 1e12]}

results = dict()

with open('data/ct_optimization_results.pickle', 'rb') as f:
    results['ct'] = pickle.load(f)

with open('data/angeli_optimization_results.pickle', 'rb') as f:
    results['angeli'] = pickle.load(f)

with open('data/karabelas_optimization_results.pickle', 'rb') as f:
    results['karabelas'] = pickle.load(f)

with open('data/karabelas_high_optimization_results.pickle', 'rb') as f:
    results['karabelas_high'] = pickle.load(f)

with open('data/sa_optimization_results.pickle', 'rb') as f:
    results['sa'] = pickle.load(f)

res = {'c1': [], 'c2': [], 'c3': [], 'c4': [], 'Re': [], 'St': [], 'Ca': [], 'We': [], 'key': []}
for key in sets:
    for data in results[key]:
        x = dict()

        c = array(data['best_fit'])
        if key == 'karabelas':
            c[:] = exp(0.1 * c[:])
        elif key == 'angeli':
            c[:] = exp(0.5 * c[:])
        else:
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

            Ca = pbe_solutions.contProperties['mu'] * pbe_solutions.Nstar * pbe_solutions.D\
                / pbe_solutions.sigma * sqrt(pbe_solutions.contProperties['rho'] / pbe_solutions.dispProperties['rho'])

        elif key == 'angeli':
            pbe_solutions = AngeliSolution(
                M=40, v0=v0, U=s['U'], phi=s['phi'], theta=s['theta'],
                model_parameters=c)
            Re = pbe_solutions.Re
            Ca = pbe_solutions.contProperties['mu'] * s['U'] / pbe_solutions.sigma

        elif key =='karabelas':
            #c[3] = 7.2e12
            pbe_solutions = \
            KarabelasSolution(
                M=40, v0=v0, U=s['U'], theta=s['theta'],
                model_parameters=c)
            Re = pbe_solutions.Re
            Ca = pbe_solutions.contProperties['mu'] * s['U'] / pbe_solutions.sigma


        elif key =='karabelas_high':
            pbe_solutions = \
            KarabelasSolutionHighViscosity(
                M=40, v0=v0, U=s['U'], theta=s['theta'],
                model_parameters=c)
            Re = pbe_solutions.Re
            Ca = pbe_solutions.contProperties['mu'] * s['U'] / pbe_solutions.sigma

        elif key == 'sa':
            pbe_solutions = \
            SASolution(
                M=20, v0=v0, U=s['U'], phi=s['phi'], theta=s['theta'],
                model_parameters=c)
            Ca = pbe_solutions.contProperties['mu'] * s['U'] / pbe_solutions.sigma
            Re = pbe_solutions.Re

        cont = pbe_solutions.contProperties
        disp = pbe_solutions.dispProperties

        print('tank volume = {0}'.format(pbe_solutions.Vt))

        St = 2. * cont['rho'] / (2. * disp['rho'] + cont['rho'])\
            * 2. / 9. * pbe_solutions.d32 ** 2. / pbe_solutions.D ** 2 * Re

        We = 2.0 * pbe_solutions.epsilon ** (2. / 3.) * cont['rho']\
            * pbe_solutions.d32 ** (5. / 3.) / pbe_solutions.sigma

        c[2] = c[2] * pbe_solutions.Vt

        res['c1'].append(c[0])
        res['c2'].append(c[1])
        res['c3'].append(c[2])
        res['c4'].append(c[3])
        res['Re'].append(Re)
        res['St'].append(St)
        res['We'].append(We)
        res['Ca'].append(Ca)
        res['key'].append(key)

        print (pbe_solutions.d32, s['d32'], (pbe_solutions.d32 - s['d32']) / s['d32'])

with open('data/nondims.pickle', 'wb') as f:
    pickle.dump(res, f)
