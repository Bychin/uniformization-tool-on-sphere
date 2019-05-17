import time

import numpy as np
from tornado import escape, web

from distribution import Distribution
from grids.spiral import SpiralGrid, SPIRAL_GRID_POINTS
from grids.classic import ClassicGrid, CLASSIC_GRID_DIV


class StatsAPI(web.RequestHandler):
    def initialize(self, database):
        self.database = database

    def __getArgs(self, body):
        data = escape.json_decode(body)

        self.points = data.get("points")
        self.mean = data.get("mean")
        self.cov = data.get("cov")

    def __validate(self):
        if self.mean is None or self.cov is None or self.points is None:
            raise ValueError("not enough arguments")

        self.points = np.reshape([float(point) for point in self.points.split(',')], (-1, 3))
        self.mean = tuple([float(coord) for coord in self.mean.split(',')])
        #length = np.sqrt(np.sum([x * x for x in self.mean]))
        #self.mean = tuple(x / length for x in self.mean)
        self.cov = tuple(float(cov) for cov in self.cov.split(','))
        # TODO len check

    def __getGrids(self):
        key = (self.mean, self.cov)

        if not (key in self.database):
            distribution = Distribution(self.mean, self.cov)

            self.database[key] = (
                SpiralGrid(SPIRAL_GRID_POINTS, distribution.calc),
                ClassicGrid(CLASSIC_GRID_DIV, distribution.calc)
            )

        self.spiral_grid, self.classic_grid = self.database[key]

    def __calculateTStats(self, value):
        return self.spiral_grid.calculateIntegralInsideIsoline(value)

    def __calculateSStats(self, point):
        pass

    def post(self):
        start = time.time()
        self.set_status(200)
        self.set_header("Access-Control-Allow-Origin", "*") # TODO

        try:
            self.__getArgs(self.request.body)
            self.__validate()

            self.__getGrids()
            print("StatsAPI: Grids generation and function evaluation time:", time.time() - start)

            start = time.time()

            t = []
            s = []

            self.pointsValues = {tuple(point) : self.spiral_grid.function(point) for point in self.points}
            for point, value in self.pointsValues.items():
                t.append(self.__calculateTStats(value))
                s.append(self.__calculateSStats(point))

            print("StatsAPI: Stats calculation time:", time.time() - start)

        except Exception as e:
            self.finish({"code": 500, "body": str(e)})
            return

        self.finish({"code": 200, "t": t, "s": s})

    # Something bad here, but it's fast workaround right now (CORS policy)
    def options(self):
        self.set_status(200)
        self.set_header("Access-Control-Allow-Origin", "*")
        self.set_header("Access-Control-Allow-Headers", "*")
        if not self.request.body:
            return
