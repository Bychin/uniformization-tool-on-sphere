import math
import numpy as np
from multiprocessing import cpu_count, Pool

import cfg


SPIRAL_GRID_POINTS = cfg.get_param("spiral_grid_points")


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
        self.function = function

        agents = 4 if cpu_count() > 4 else cpu_count()
        chunksize = math.ceil(len(self.points) / agents)
        with Pool(processes=agents) as pool:
            values = pool.map(function, self.points, chunksize)
        self.values = {self.points[i] : values[i] for i in range(len(self.points))}

        self.data = np.fromiter(self.values.values(), dtype=float)

    def calculateIntegralInsideIsoline(self, c):
        valuesInsideIsoline = np.array(list(filter(lambda x: x > c, self.data)))

        if len(valuesInsideIsoline) == 0:  # if c is pdf's maximum
            return 0
        return np.mean(valuesInsideIsoline) * self.k * len(valuesInsideIsoline)
