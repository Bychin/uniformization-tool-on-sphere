from tornado.ioloop import IOLoop
from tornado import web
from tornado.httpserver import HTTPServer


SPIRAL_GRID_POINTS = 100000


class MainHandler(web.RequestHandler):
    def get(self):
        self.write("Hello, world")


class IsolineAPI(web.RequestHandler):
    def _validate(self):
        if self.mean is None or self.cov is None or self.ratio is None:
            raise ValueError("not enough arguments")

        self.mean = [int(coord) for coord in self.mean.split(',')]
        self.cov = [int(cov) for cov in self.cov.split(',')]
        self.ratio = int(self.ratio)

    def get(self):
        try:
            self.mean = self.get_argument("mean", None, True)
            self.cov = self.get_argument("cov", None, True)
            self.ratio = self.get_argument("ratio", None, True)

            self._validate()

        except Exception as e:
            self.finish({"code": 400, "body": str(e)})
            return

        # calc ggp on uni grid

        # [0|1/2|1] in cycle

        # store calcs on uni (also make std) to map[uid]
        # return c and uid

        self.finish({"code": 200, "body": "test"})


def make_app():
    return web.Application([
        (r"/api/isoline", IsolineAPI),
        (r"/", MainHandler),
    ])


if __name__ == "__main__":
    # TODO init grids to as global vars

    app = make_app()
    server = HTTPServer(app)
    server.listen(8080)
    IOLoop.current().start()
