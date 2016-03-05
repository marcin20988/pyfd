from case_class import CaseSolution
from numpy import zeros, pi


class KarabelasSolution(CaseSolution):
    def __init__(
            self,
            M=20,
            U=1.1,  # [rps] impeller revolutions
            v0=5e-10,  # [cm^3]
            model_parameters=None,
            theta=600.):

        # pipe diamter and length
        self.D = 5.04e-02
        self.L = 32.0

        contProperties = dict()
        dispProperties = dict()
        domainProperties = dict()

        # oil
        contProperties['mu'] = 1.8e-3  # [P = kg * m^-1 s^-1]
        contProperties['rho'] = 798.  # [kg/cm3]

        # calculate turbulent properties
        Re = U * self.D / contProperties['mu'] * contProperties['rho']
        self.Re = Re
        I = 0.16 * Re ** (-1. / 8.)
        u_rms = U * I
        k = 3. / 2. * u_rms ** 2
        L_t = 0.038 * self.D
        contProperties['epsilon'] = 0.09 * k ** (3./2.) / L_t
        contProperties['Re'] = Re
        # water solution
        dispProperties['sigma'] = 33.1e-03  # [P = kg * m^-1 s^-1]
        dispProperties['rho'] = 1000.  # [kg/m3]
        dispProperties['phi'] = 0.2

        dispProperties['vMax'] = v0 * 3.

        # Feed distribution
        dispProperties['v0'] = v0
        dispProperties['sigma0'] = v0 / 10

        # Feed
        domainProperties['theta'] = theta
        domainProperties['V'] = pi * self.L * (self.D / 2) ** 2
        domainProperties['M'] = M

        CaseSolution.__init__(
            self, dispProperties, contProperties, domainProperties,
            model_parameters=model_parameters)

class KarabelasSolutionHighViscosity(CaseSolution):
    def __init__(
            self,
            M=20,
            U=1.1,  # [rps] impeller revolutions
            v0=5e-10,  # [cm^3]
            model_parameters=None,
            theta=600.):

        # pipe diamter and length
        self.D = 5.04e-02
        self.L = 32.0

        contProperties = dict()
        dispProperties = dict()
        domainProperties = dict()

        # oil
        contProperties['mu'] = 16.0e-3  # [P = kg * m^-1 s^-1]
        contProperties['rho'] = 890.  # [kg/cm3]

        # calculate turbulent properties
        Re = U * self.D / contProperties['mu'] * contProperties['rho']
        I = 0.16 * Re ** (-1. / 8.)
        u_rms = U * I
        k = 3. / 2. * u_rms ** 2
        L_t = 0.038 * self.D
        contProperties['epsilon'] = 0.09 * k ** (3./2.) / L_t
        contProperties['Re'] = Re
        # water solution
        dispProperties['sigma'] = 34.0e-03
        dispProperties['rho'] = 1000.  # [kg/m3]
        dispProperties['phi'] = 0.2

        dispProperties['vMax'] = v0 * 3.

        # Feed distribution
        dispProperties['v0'] = v0
        dispProperties['sigma0'] = v0 / 10

        # Feed
        domainProperties['theta'] = theta
        domainProperties['V'] = pi * self.L * (self.D / 2) ** 2
        domainProperties['M'] = M

        CaseSolution.__init__(
            self, dispProperties, contProperties, domainProperties,
            model_parameters=model_parameters)
