import numpy as np
from numpy import arange, sum, exp, linspace, sqrt, pi, zeros
import sys


class fluid:
    """
        class to store properties needed for breakup and coalescence models

    """

    def __init__(self, name, caseNr=0, caseNr2=0):
        if name is "galinat":
            self.rhoc = 996.0
            self.rhod = 683.7
            self.mud = 4.5e-04
            self.muc = 8.2e-04
            self.sigma = 4.7e-02
            self.dpMax = np.array([173.0, 363.0, 566.0, 706.0, 871.0, 1120.0])
            self.U = np.array([0.118, 0.177, 0.236, 0.275, 0.314, 0.354])
            self.Re = np.array([4250, 6400, 8600, 10000, 11500, 12900])
            # orifice ratio
            self.beta_or = 0.5
            # volume fraction
            self.alpha = 0.02
            # pipe diameter, and length
            self.D = 0.03
            self.L = 1.0
            # residence time; this is a though one...
            # we'll take it equal to the time needed for the
            # mean flow to travel length equal to three time
            # the orifice thickness
            # orifice thickness is 5mm
            self.thetas = 2.0 * 0.005 / self.U[:]
            self.epsilons = 1.0 / self.rhoc * self.dpMax[:] * self.U[:]\
                / 2.0 / self.D * (1.0 / self.beta_or ** 2 - 1.0)
            self.epsilon = self.epsilons[caseNr]
            self.theta = self.thetas[caseNr]
            self.V = pi * (self.D / 2.0) ** 2 * self.L
        elif name is "simmonsAzzopardi":
            self.rhoc = 797.0
            self.muc = 1.8e-03
            self.rhod = 1166.0
            self.mud = 1.6e-03
            self.sigma = 1.0e-02
            self.U = 2.71
            # volume fraction
            self.alpha = 0.117
            # pipe diameter, length and volume
            self.D = 0.063
            self.L = 4.5
            self.V = pi * (self.D / 2.0) ** 2 * self.L
            self.epsilon = 0.082
            self.theta = None
        elif name is "coulaloglou":
            self.rhoc = 1000.0
            self.muc = 1.0e-03
            self.rhod = 972.0
            self.mud = 1.3e-03
            self.sigma = 42.82e-03
            # volume fraction
            self.alphas = np.array([0.05, 0.1, 0.15])
            self.alpha = self.alphas[caseNr]
            # impeller speed: (in 1 / min)
            self.Nstars = np.array([190.0, 220.0, 250.0, 280.0, 310.0])
            # convert to 1 / second:
            self.Nstar = self.Nstars[caseNr2] / 60.0
            # impeller diameter: (10cm)
            self.Dstar = 0.1
            # tank volume (12l)
            self.V = 12.0e-03
            self.epsilon = 0.407 * self.Nstar ** 3 * self.Dstar ** 2
            # residence time is 10 minutes
            self.theta = 10.0 * 60.0
        else:
            sys.exit("Valid cases are: 'galinat', 'simmonsAzzopardi', 'coulaloglou'")

        # default values from Coulaloglou and Tavlarides
        self.C1 = 0.00487
        self.C2 = 0.08
        self.C3 = 2.17e-16
        self.C4 = 2.28e13

    def gamma(self, xi):
        C = self.C1 * xi ** (-2.0 / 9.0) * self.epsilon ** (1.0 / 3.0)\
            / (1.0 + self.alpha)
        exp_argument = - self.C2 * self.sigma * (1.0 + self.alpha) ** 2 \
            / (self.rhod * xi ** (5.0 / 9.0) * self.epsilon ** (2.0 / 3.0))
        return C * exp(exp_argument)

    # droplet daughter distribution:
    def beta(self, xi2, xi1):
        return 2.0 * 2.4 / xi1 * exp(- 4.5 * (2.0 * xi2 - xi1) ** 2 / xi1 ** 2)

    # coalescence rate:
    def Q(self, xi1, xi2):
        dRatio = xi1 ** (1.0 / 3.0) * xi2 ** (1.0 / 3.0)\
            / (xi1 ** (1.0 / 3.0) + xi2 ** (1.0 / 3.0))
        dRatio = dRatio ** 4

        exp_argument = - self.C4 * self.muc * self.rhoc * self.epsilon\
            / (1.0 + self.alpha) ** 3 * dRatio / self.sigma ** 2

        C = self.C3 * (xi1 ** (1.0 / 3.0) + xi2 ** (1.0 / 3.0)) ** 2\
            * (xi1 ** (2.0 / 9.0) + xi2 ** (2.0 / 9.0)) ** 0.5\
            * self.epsilon ** (1.0 / 3.0) / (1.0 + self.alpha) / self.V
        return exp(exp_argument) * C
