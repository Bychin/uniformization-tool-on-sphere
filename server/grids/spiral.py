import math
import numpy as np


class SpiralGrid:
    def __init__(self, N, function):
        self.__generateGrid(N)
        self.__calculateFunc(function)

        self.k = 4 * math.pi / N  # Area of an elementary part of a sphere for every grid value

    def __generateGrid(self, N):
        self.points = []
        theta = 0

        for k in range(1, N+1):
            h = -1 + (2 * k - 2) / (N - 1)
            phi = math.acos(h)

            s_phi = math.sin(phi)
            c_phi = math.cos(phi)

            if k != 1 and k != N:
                theta = (theta + 3.8 / math.sqrt(N * (1 - h * h))) % (2 * math.pi)
            else:
                theta = 0

            s_theta = math.sin(theta)
            c_theta = math.cos(theta)

            x = s_phi * c_theta
            y = s_phi * s_theta
            z = c_phi

            self.points.append((x, y, z))

    def __calculateFunc(self, function):
        self.function = {point : function(point) for point in self.points}  # TODO values
        self.data = np.fromiter(self.function.values(), dtype=float)

    def calculateIntegralInsideIsoline(self, c):
        valuesInsideIsoline = np.array(list(filter(lambda x: x > c, self.data))) #self.data[self.data > c]
        return np.mean(valuesInsideIsoline) * self.k * len(valuesInsideIsoline)
