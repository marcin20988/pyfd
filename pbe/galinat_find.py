from numpy import arange, sum, exp, linspace, sqrt, pi, zeros, genfromtxt
from itertools import cycle
from moc import MOCSolution
import matplotlib.pyplot as plt
from fluid import fluid
import os.path
import numpy as np
import sys

pbe_solutions = dict()
cases = arange(5)

minError = 1.0e12
C1Best = - 1.0
C2Best = - 1.0

C1s = np.array([0.005, 0.5, 0.025])
C2s = np.array([0.005, 0.1, 0.005])

for C1 in C1s:
    for C2 in C2s:

        sufix = '_C1=' + repr(C1) + '_C2=' + repr(C2)
        print sufix
#sys.exit('exit')

        for c in cases:
            F = fluid('galinat', c)
            # expression for epsilon proposed by Galinat et al seems to be
            # with correct order of magnitude, however I find it about
            # twice as big as it should be (comparing to 3D CFD)
            #F.epsilon /= 2.0
            F.C1 = C1
            F.C2 = C2
            #print F.V

            # 20 times residence time should be enough to get convergence
            t = arange(0.0, 5.0 * F.theta, 1e-01)
            # data structure: each row contains pair of values Re, d32
            # cases are aranged with increasing Re
            d0s = genfromtxt('validationData/orifice/d32_upstream.txt',
                             delimiter=',')
            d0 = d0s[c][1] * 1e-03
            # convert diameter to volume
            v0 = pi / 6.0 * d0 ** 3
            # standard deviation
            s0 = v0 / 8.0
            # this is a breakup dominated case so we don't need large diameters
            vmax = 1.4 * v0
            g = 40

            dv = vmax / g
            v = dv + dv * arange(g)
            Ninit = F.alpha / v0 * F.V\
                * 1.0 / s0 / sqrt(2.0 * pi)\
                * exp(- (v - v0) ** 2 / 2 / s0 ** 2)
            pbe_solutions[c] = MOCSolution(
                Ninit, t, dv,
                gamma=F.gamma,
                beta=F.beta,
                theta=F.theta,
                pdf='number'
            )

        m1_init = sum(Ninit[:] * v[:])
        normInit = sum(Ninit)
        meanV_init = m1_init / normInit
        print "initial m1: ", m1_init
        print "initial d: ", (6.0 * meanV_init / pi) ** (1.0 / 3.0)

        dRe = zeros(cases.shape)
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
            dRe[n] = meanD

            fname = 'validationData/orifice/' + repr(n) + '/downstream.txt'
            if os.path.isfile(fname):
                data = genfromtxt(fname, delimiter=',')
                # experimental and numerical data are discretised on different
                # number of point so to present them on a single plot they need
                # to be multiplied by appropriate ratio of node points
                y = data.T[1] * (float)(data.shape[0]) / g
                x = data.T[0]
                w = 0.1

                ax.bar(
                    x, y, width=w, color=(0.8, 0.8, 0.8, 0.8),
                    label="experiments".format(n))
            ax.grid()
            ax.set_xlabel('diameter[mm]')
            ax.set_ylabel('frequency')

            plt.savefig('validationData/sweep/orifice/hist_' + sufix
                        + '___' + repr(n) + '.pdf')

        exp_dRe = genfromtxt('validationData/orifice/d32_downstream.txt',
                             delimiter=',')

        fig = plt.figure()
        ax = fig.gca()

        y = exp_dRe.T[1]
        ax.plot(F.Re[0:-1], y, '+', marker='+', label="Galinat et al")
        ax.plot(F.Re[0:dRe.shape[0]], dRe * 1e03, '+', marker='o',
                label='PBE')

        error = 0.0
        for iError in arange(dRe.shape[0]):
            error += (y[iError] - dRe[iError] * 1e03) ** 2
        if error < minError:
            C1Best = C1
            C2Best = C2
            minError = error
        ax.grid()
        plt.savefig('validationData/sweep/orifice/dRe' + sufix + '.pdf')

        fig = plt.figure()
        ax = fig.gca()
        ax.plot(
            d, F.gamma(v), "+",
            marker=next(markers),
            label="breakup frequency".format(n))
        ax.legend(loc="best")
        ax.set_yscale('log')
        plt.savefig('validationData/sweep/orifice/breakup' + sufix + '.pdf')

        #plt.show()

print "minimum error: ", minError
print "best match for C1: ", C1Best
print "best match for C2: ", C2Best
