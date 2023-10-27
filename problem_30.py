'''
Plot Debye screening length vs concentration:

lam_D = 1/k_D = eps*kB*T/(e^2 (z_+^2*rho_+ + z_-^2*rho_-))

Plot Activity Coefficient vs concentration:

gamma = exp(-(z*e)^2 * k_D / (8pi*eps*kB*T))
'''

import numpy as np
import matplotlib.pyplot as plt


class problem_30:

    def __init__(self):

        # constants in SI units
        self.T = 25.0+273.15
        self.kB = 1.380649e-23
        self.kBT = self.kB*self.T
        self.e = 1.60218e-19
        self.eps = 78.54*8.854187e-12
        self.N_A = 6.0221408e23
        self.Liter2M3 = 0.001

        # concentration in mol/L
        self.C = np.linspace(0.0, 0.1, 200)

    def k_D(self):

        # Since Na+ and Cl-, so z_+ = z_- = 1, and rho_+ + rho_- = rho0

        rho = self.C*self.N_A/self.Liter2M3
        k_D = np.sqrt(2.0*rho*self.e**2/(self.eps*self.kBT))
        return k_D

    def gamma(self):

        k_D = self.k_D()
        gamma = np.exp(-k_D*self.e**2/(8.0*np.pi*self.eps*self.kBT))        
        return gamma

    def plot(self):

        fig = plt.figure(figsize=(16, 9))
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)

        lam_D = 1.0/self.k_D()
        gamma = self.gamma()

        ax1.plot(self.C, lam_D, label=r"$\lambda_D$")
        ax1.set_ylabel(r"$\lambda_D [m]$")
        ax1.set_xlabel(r"$Concentration [mol/L]$")
        ax2.plot(self.C, gamma, label=r"$\gamma_{\pm}$")
        ax2.set_ylabel(r"$\gamma_{\pm}$")
        ax2.set_xlabel(r"$Concentration [mol/L]$")
        ax1.legend()
        ax2.legend()
        plt.show()
        return fig


prob30 = problem_30()
prob30.plot()
