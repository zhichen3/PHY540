'''
Part a)
Make Poisson distribution for mu = 100.
Compare it with Gaussian distribution for
mean = 100, and sigma = 10
'''

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import factorial

class problem_10:
    def __init__(self, mu, mean, sigma):

        # mu: average number of events in interval T for Poisson
        self.mu = mu

        # mean: mean for Gaussian
        self.mean = mean

        # sigma: standard deviation for Gaussian
        self.sigma = sigma

    def _draw_poisson(self, size):
        # Draw samples from Poisson distribution

        generator = np.random.default_rng()
        return generator.poisson(lam=self.mu, size=size)

    def _gaussian(self, n):
        # Compute the Gaussian distribution given an array of n

        return np.exp(-0.5*((n - self.mean)/self.sigma)**2) / (self.sigma*np.sqrt(2.0*np.pi))

    def _plot(self):
        # Plot Poisson Distribution and Gaussian

        fig, ax1 = plt.subplots()

        poisson_samples = self._draw_poisson(1000)
        ax1.hist(poisson_samples, bins=20, density=True, label="Poisson")

        x_gaussian = np.linspace(25.0, 175.0, 1000)
        y_gaussian = self._gaussian(x_gaussian)

        ax1.legend(loc="upper right")
        ax1.plot(x_gaussian, y_gaussian, color='red', label="Gaussian")
        # ax2=ax1.twinx()
        # ax2.plot(x_gaussian, y_gaussian, color='red', label="Gaussian")
        # ax2.set_ylim([0.0, y_gaussian.max()])
        # ax2.legend(loc="upper left")

        plt.show()
        return fig

    def part_a(self):
        # Do Part a of the problem
        
        self._plot()


prob_10 = problem_10(100, 100, 10)
prob_10.part_a()

        
