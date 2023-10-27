'''
Make an Arrhenius-style plot of p0 on a log scale vs 1/T
'''

import matplotlib.pyplot as plt
import numpy as np

def plot_pressure_vs_iTemp(eps, Temp_bounds):
    # Arrhenius-style plot for pressure vs temperature assuming coverage is 1/2

    # Define constants

    # boltzmann constant
    kB = 1.38e-23

    # planck constant
    h = 6.626e-34

    # mass of argon
    m = 39.948 * 1.66054e-27

    # Binding energy
    eps = eps

    # inverse Temperature with kB
    iTemp = 1.0 / (np.logspace(np.log10(Temp_bounds[0]), np.log10(Temp_bounds[1]), 1000))

    beta = iTemp / kB
    # thermal de Broglie wavelength
    lam = np.sqrt(h*h*beta / (2.0*np.pi*m))

    # pressure [Torr]
    p = np.exp(-eps*beta) / (beta*lam**3) * 0.00750062

    fig, ax = plt.subplots()
    ax.plot(iTemp, p)
    ax.set_xlabel(r"$T^{-1} [K^{-1}]$")
    ax.set_ylabel(r"$p_0$ [Torr]")
    ax.set_yscale("log")
    ax.set_xticks(np.linspace(iTemp[0], iTemp[-1], 8))
    ax.grid()
    plt.show()

# 99 meV
eps = 99.0* 1.60218e-22
temperature_bound = [30, 80]
plot_pressure_vs_iTemp(eps, temperature_bound)
