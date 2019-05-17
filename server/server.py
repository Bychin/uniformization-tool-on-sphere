from tornado.ioloop import IOLoop
from tornado import web
from tornado.httpserver import HTTPServer

from api.isoline import IsolineAPI
from api.stats import StatsAPI


def make_app(database):
    return web.Application([
        (r"/api/isoline", IsolineAPI, dict(database=database)),
        (r"/api/stats", StatsAPI, dict(database=database)),
    ])


if __name__ == "__main__":
    database = {}

    app = make_app(database)
    server = HTTPServer(app)
    server.listen(8080)
    IOLoop.current().start()
