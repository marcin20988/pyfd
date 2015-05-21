import numpy as np


class fluid:
    """
        class to store properties needed for breakup and coalescence models

    """

    def __init__(self, name, caseNr=0):
        if name is "galinat":
            self.rhoc = 996.0
            self.rhod = 683.7
            self.mud = 4.5e-04
            self.muc = 8.2e-04
            self.sigma = 4.7e-02
            self.dpMax = np.array([173.0, 363.0, 566.0, 706.0, 871.0, 1120.0])
            self.U = np.array([0.118, 0.177, 0.236, 0.275, 0.314, 0.354])
            # orifice ratio
            self.beta = 0.5
            # volume fraction
            self.alpha = 0
            # pipe diameter
            self.D = 0.03
            # residence time; this is a though one...
            # we'll take it equal to the time needed for the
            # mean flow to travel length equal to three time
            # the orifice thickness
            # orifice thickness is 5mm
            self.thetas = 3.0 * 0.005 / self.U[:]
            self.epsilons = 1.0 / self.rhoc * self.dpMax[:] * self.U[:]\
                / 2.0 / self.D * (1.0 / self.beta ** 2 - 1.0)
            self.epsilon = self.epsilons[caseNr]
            self.theta = self.thetas[caseNr]
        else:
            sys.exit("Valid cases are: 'galinat'")
