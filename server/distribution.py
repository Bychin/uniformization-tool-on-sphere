import math
import numpy as np
from scipy.special import erf

import autograd.numpy as anp
from autograd.scipy.special import erf as aerf

# TODO add description and 3-D checks 
class Distribution():
    def __init__(self, mean, cov):
        self.mean = mean
        self.lam = self.__getFullInvCovMatrix(cov) # lam is inverse matrix of cov
        self.lam_det = np.linalg.det(self.lam)
        self.mean_norm = math.sqrt(self.__lamInnerProduct(self.mean, self.mean))
        self.max = (None, -1)  # point of PDF's maximum and it's value

    def __getFullInvCovMatrix(self, cov):
        cov_matrix = np.array([[cov[0], cov[1], cov[2]], [cov[1], cov[3], cov[4]], [cov[2], cov[4], cov[5]]])
        return np.linalg.inv(cov_matrix)

    def __lamInnerProduct(self, x, y):
        return np.sum([[self.lam[i][j] * x[i] * y[j] for j in range(3)] for i in range(3)])
    
    # u should be 3-D vector
    def calc(self, u):
        u_norm = math.sqrt(self.__lamInnerProduct(u, u))
        z = self.__lamInnerProduct(self.mean, u) / u_norm

        c = math.exp(-0.5 * self.mean_norm * self.mean_norm) * math.sqrt(self.lam_det) / (4 * math.pi * (u_norm ** 3))
        value = c * (z * math.sqrt(2 / math.pi) + math.exp(0.5 * z * z) * (1 + z * z) * (1 + erf(z / math.sqrt(2))))

        return value

    def calcForGradient(self, u):
        u_norm = anp.sqrt(anp.sum(anp.array([anp.array([self.lam[i][j] * u[i] * u[j] for j in range(3)]) for i in range(3)])))
        z = anp.sum(anp.array([anp.array([self.lam[i][j] * self.mean[i] * u[j] for j in range(3)]) for i in range(3)])) / u_norm

        c = anp.exp(-0.5 * self.mean_norm * self.mean_norm) * anp.sqrt(self.lam_det) / (4 * anp.pi * (u_norm ** 3))
        value = c * (z * anp.sqrt(2 / anp.pi) + anp.exp(0.5 * z * z) * (1 + z * z) * (1 + aerf(z / anp.sqrt(2))))

        return value
