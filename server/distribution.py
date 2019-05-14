import math
import numpy as np
from scipy.special import erf


# TODO add description and 3-D checks 
class Distribution():
    def __init__(self, mean, cov):
        self.mean = mean
        self.lam = self.__getFullInvCovMatrix(cov) # lam is inverse matrix of cov
        self.lam_det = np.linalg.det(self.lam)

    def __getFullInvCovMatrix(self, cov):
        cov_matrix = np.array([[cov[0], cov[1], cov[2]], [cov[1], cov[3], cov[4]], [cov[2], cov[4], cov[5]]])
        return np.linalg.inv(cov_matrix)

    def __lamInnerProduct(self, x, y):
        return np.sum([[self.lam[i][j] * x[i] * y[j] for j in range(3)] for i in range(3)])
    
    # u should be 3-D vector
    def calc(self, u):
        u_norm = math.sqrt(self.__lamInnerProduct(u, u))
        m_norm = math.sqrt(self.__lamInnerProduct(self.mean, self.mean))
        z = self.__lamInnerProduct(self.mean, u) / u_norm

        c = math.exp(-0.5 * m_norm * m_norm) * math.sqrt(self.lam_det) / (4 * math.pi * (u_norm ** 3))
        return c * (z * math.sqrt(2 / math.pi) + math.exp(0.5 * z * z) * (1 + z * z) * (1 + erf(z / math.sqrt(2))))
