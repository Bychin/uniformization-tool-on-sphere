import math
import numpy as np
import time

from tornado import web

from distribution import Distribution
from grids.spiral import SpiralGrid, SPIRAL_GRID_POINTS
from grids.classic import ClassicGrid, CLASSIC_GRID_DIV


class IsolineAPI(web.RequestHandler):
    def initialize(self, database):
        self.database = database

    def __getArgs(self):
        self.mean = self.get_argument("mean", None, True)
        self.cov = self.get_argument("cov", None, True)
        self.ratios = self.get_argument("ratio", None, True)

    def __validate(self):
        if self.mean is None or self.cov is None: # or self.ratios is None:
            raise ValueError("not enough arguments")

        self.mean = tuple([float(coord) for coord in self.mean.split(',')])
        #length = np.sqrt(np.sum([x * x for x in self.mean]))
        #self.mean = tuple(x / length for x in self.mean)
        self.cov = tuple(float(cov) for cov in self.cov.split(','))
        # TODO len check

        self.ratios = [float(ratio) for ratio in self.ratios.split(',')]
        # TODO <1 etc checks

    def __getGrids(self):
        key = (self.mean, self.cov)

        if not (key in self.database):
            distribution = Distribution(self.mean, self.cov)

            self.database[key] = (
                SpiralGrid(SPIRAL_GRID_POINTS, distribution.calc),
                ClassicGrid(CLASSIC_GRID_DIV, distribution.calc)
            )

        self.spiral_grid, self.classic_grid = self.database[key]

    def __getConstForRatio(self, ratio):
        f1 = 0
        f2 = np.max(self.spiral_grid.data)
        c = 0

        integral = 10  # TODO
        eps = 0.001
        i = 0
        while math.fabs(integral - ratio) > eps:
            i += 1
            if i > 50:  # TODO
                break

            c = (f1 + f2) / 2
            integral = self.spiral_grid.calculateIntegralInsideIsoline(c)

            if integral > ratio:
                f1 = c
            else:
                f2 = c

        return c

    def __getIsolines(self):
        self.isolines = {}

        for ratio in self.ratios:
            c = self.__getConstForRatio(ratio)
            points = self.classic_grid.getIsolineCoords(c)
            self.isolines[ratio] = points

    def get(self):
        start = time.time()
        self.set_status(200)
        self.set_header("Access-Control-Allow-Origin", "*") # TODO

        try:
            self.__getArgs()
            self.__validate()

            self.__getGrids()
            print("IsolineAPI: Grids generation and function evaluation time:", time.time() - start)

            start = time.time()
            self.__getIsolines()
            print("IsolineAPI: Isolines evaluation time:", time.time() - start)

        except Exception as e:
            self.finish({"code": 500, "body": str(e)})
            return

        self.finish({"code": 200, "isolines": self.isolines})
