from tornado.ioloop import IOLoop
from tornado import web
from tornado.httpserver import HTTPServer
import time
import numpy as np
import math
from multiprocessing import Process, Queue

from distribution import Distribution
from grids.spiral import SpiralGrid
from grids.classic import ClassicGrid


SPIRAL_GRID_POINTS = 100000
CLASSIC_GRID_DIV = 360


class IsolineAPI(web.RequestHandler):
    def initialize(self, database):
        self.database = database
        self.response_body = {} # TODO remove?

    def _getArgs(self):
        self.mean = self.get_argument("mean", None, True)
        self.cov = self.get_argument("cov", None, True)
        self.ratios = self.get_argument("ratio", None, True)

    def _validate(self):
        if self.mean is None or self.cov is None: # or self.ratios is None:
            raise ValueError("not enough arguments")

        self.mean = tuple([float(coord) for coord in self.mean.split(',')])
        #length = np.sqrt(np.sum([x * x for x in self.mean]))
        #self.mean = tuple(x / length for x in self.mean)
        self.cov = tuple(float(cov) for cov in self.cov.split(','))
        # TODO len check

        self.ratios = [float(ratio) for ratio in self.ratios.split(',')]
        # TODO <1 etc checks

    def _getGrids(self):
        key = (self.mean, self.cov)

        if not (key in self.database):
            distribution = Distribution(self.mean, self.cov)

            self.database[key] = (
                SpiralGrid(SPIRAL_GRID_POINTS, distribution.calc),
                ClassicGrid(CLASSIC_GRID_DIV, distribution.calc)
            )

        self.spiral_grid, self.classic_grid = self.database[key]

    def _getConstForRatio(self, ratio):
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

    def _getIsolines(self):
        self.isolines = {}

        for ratio in self.ratios:
            c = self._getConstForRatio(ratio)
            points = self.classic_grid.getIsolineCoords(c)
            self.isolines[ratio] = points

    def get(self):
        start = time.time()
        self.set_status(200)
        self.set_header("Access-Control-Allow-Origin", "*") # TODO

        try:
            self._getArgs()
            self._validate()

            self._getGrids()
            print("Grids generation and function evaluation time:", time.time() - start)

            start = time.time()
            self._getIsolines()
            print("Isolines evaluation time:", time.time() - start)

        except Exception as e:
            self.finish({"code": 500, "body": str(e)})
            return

        self.finish({"code": 200, "isolines": self.isolines})


def make_app(database):
    return web.Application([
        (r"/api/isoline", IsolineAPI, dict(database=database)),
    ])


if __name__ == "__main__":
    database = {}

    app = make_app(database)
    server = HTTPServer(app)
    server.listen(8080)
    IOLoop.current().start()
