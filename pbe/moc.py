from numpy import arange, zeros
from scipy.integrate import odeint
import sys

"""
Method of classes
"""


class MOCSolution:
    """
    Based on Brooks and Hidy uniform discretisation

    """
    def RHS(
        self, N, t
    ):
        dNdt = zeros(self.number_of_classes)

        if self.gamma is not None and self.beta is not None:
            # Death breakup term
            dNdt -= N * self.gamma(self.xi)
            # Birth breakup term
            for i in arange(self.number_of_classes):
                for j in arange(i + 1, self.number_of_classes):
                    dNdt[i] = \
                        self.beta(self.xi[i], self.xi[j]) \
                        * self.gamma(self.xi[j]) \
                        * N[j] * self.delta_xi

        if self.Q is not None:
            for i in arange(self.number_of_classes):
                # Birth coalescence term
                for j in arange(1, i):
                    dNdti = 0.5 * N[i - j] * N[j] \
                        * self.Q(self.xi[j], self.xi[i - j])
                    if self.pdf == "number":
                        dNdt[i] += dNdti
                    elif self.pdf == "density":
                        dNdt[i] += dNdti * self.delta_xi
                    else:
                        sys.exit(
                            "Available pdf types are 'density' and 'number'")
                # Death coalescence term
                for j in arange(self.number_of_classes):
                    dNdti = N[i] * N[j] * self.Q(self.xi[i], self.xi[j])
                    if self.pdf == "number":
                        dNdt[i] -= dNdti
                    elif self.pdf == "density":
                        dNdt[i] -= dNdti * self.delta_xi
                    else:
                        sys.exit(
                            "Available pdf types are 'density' and 'number'")

        if self.theta is not None:
            dNdt -= (N - self.N0) / self.theta
        return dNdt

    def __init__(self, N0, t, xi0, beta=None, gamma=None, Q=None, theta=None,
                 pdf="number"):
        self.number_of_classes = N0.shape[0]
        # Kernels setup
        self.beta = beta  # Daughter particle distribution
        self.gamma = gamma  # Breakup frequency
        self.Q = Q  #
        # inflow and outflow replaced with relaxation to equilibrium
        # process with relaxation time equal to residence time theta
        self.theta = theta
        self.N0 = N0
        # choose pdf formulation: number or number density
        self.pdf = pdf
        # Uniform grid
        self.xi = xi0 + xi0 * arange(self.number_of_classes)
        self.delta_xi = xi0
        # Solve procedure
        self.N = odeint(lambda NN, t: self.RHS(NN, t), N0, t)
