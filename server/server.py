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
CLASSIC_GRID_DIV = 180


class IsolineAPI(web.RequestHandler):
    def initialize(self, database):
        self.database = database
        self.response_body = {} # TODO remove?

    def _getArgs(self):
        self.mean = self.get_argument("mean", None, True)
        self.cov = self.get_argument("cov", None, True)
        self.ratio = self.get_argument("ratio", None, True)

    def _validate(self):
        if self.mean is None or self.cov is None or self.ratio is None:
            raise ValueError("not enough arguments")

        self.mean = tuple([float(coord) for coord in self.mean.split(',')])
        length = np.sqrt(np.sum([x * x for x in self.mean]))
        self.mean = tuple(x / length for x in self.mean)
        self.cov = tuple(float(cov) for cov in self.cov.split(','))
        # TODO len check

        self.ratio = float(self.ratio)

    def _getGrids(self):
        key = (self.mean, self.cov)

        if not (key in self.database):
            distribution = Distribution(self.mean, self.cov)

            spiral_grid = SpiralGrid(SPIRAL_GRID_POINTS, distribution.calc)
            classic_grid = ClassicGrid(CLASSIC_GRID_DIV, distribution.calc)

            self.database[key] = (spiral_grid, classic_grid)

        return self.database[key]

    def get(self):
        start = time.time()
        self.set_status(200)
        self.set_header("Access-Control-Allow-Origin", "*") # TODO

        #try:
        self._getArgs()
        self._validate()
        spiral_grid, classic_grid = self._getGrids()

        #except Exception as e:
        #    self.response_body["code"] = 400
        #    self.response_body["body"] = str(e)
        #    self.finish(self.response_body)
        #    return

        # [0|1/2|1] in cycle

        print(time.time() - start)

        ### Integral calculations ###
        # k = 4 * math.pi / SPIRAL_GRID_POINTS
        # data = [spiral_grid.function[point] for point in spiral_grid.points]
        # print(np.mean(data) * k * len(data))

        self.finish({"code": 200, "trapeziums": classic_grid.trapeziums})


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
