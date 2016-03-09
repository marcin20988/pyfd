import pickle
from pyfd.pbe.moc import CTSolution, AngeliSolution, KarabelasSolution, KarabelasSolutionHighViscosity, SASolution
from numpy import exp, array, pi, sqrt

sets = ['sa', 'angeli', 'ct', 'karabelas']

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

        # raw average
        c = [0.2592, 0.1073, 1.6894, 4.9410e13]

        # ct only
        #c = [0.09652, 0.1282, 0.022, 7.198e13]

        # karabelas only
        #c = [0.1117, 0.1456, 0.067, 1.27e13]

        s = data['setup']
        v0 = pi / 6 * s['d32']**3

        if key == 'ct':
            c[2] /= 0.012
            pbe_solutions = \
                CTSolution(
                    M=40, v0=v0, Nstar=s['Nstar'], phi=s['phi'],
                    model_parameters=c)


        elif key == 'angeli':
            c[2] /= 0.0042977
            pbe_solutions = AngeliSolution(
                M=40, v0=v0, U=s['U'], phi=s['phi'], theta=s['theta'],
                model_parameters=c)

        elif key =='karabelas':
            c[2] /= 0.06384
            pbe_solutions = \
            KarabelasSolution(
                M=40, v0=v0, U=s['U'], theta=s['theta'],
                model_parameters=c)

        elif key =='karabelas_high':
            c[2] /= 0.06384
            pbe_solutions = \
            KarabelasSolutionHighViscosity(
                M=40, v0=v0, U=s['U'], theta=s['theta'],
                model_parameters=c)

        elif key == 'sa':
            c[2] /= 0.0140276
            pbe_solutions = \
            SASolution(
                M=40, v0=v0, U=s['U'], phi=s['phi'], theta=s['theta'],
                model_parameters=c)

        print key, (pbe_solutions.d32, s['d32'], (pbe_solutions.d32 - s['d32']) / s['d32'])
