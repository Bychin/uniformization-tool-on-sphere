import time
import math

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

    # TODO classmethod?
    def __distanceBetweenPoints(self, a, b):
        return math.sqrt(np.sum([(a[i] - b[i]) ** 2 for i in range(3)]))

    def __insertPointIntoIsoline(self, point, isoline_points):
        point_index = None

        distances = np.array([self.__distanceBetweenPoints(iso_point, point) for iso_point in isoline_points])
        min_distances_by_indices = distances.argsort()[:2]
        indices_diff = int(math.fabs(min_distances_by_indices[0] - min_distances_by_indices[1]))
        if indices_diff == len(isoline_points) - 1:
            point_index = len(isoline_points)
        elif indices_diff == 1:
            point_index = np.max(min_distances_by_indices)
        else:
            print("WARNING: 2 min distances are not nearby!", isoline_points, distances, min_distances_by_indices, point)
            raise False  # TODO
        isoline_points.insert(point_index, point)  # TODO check it by drawing 3 points explicitly

        return point_index, isoline_points

    def __calculateSStats(self, point, value):
        isoline_points = self.classic_grid.getIsolineCoords(value)
        point_index, isoline_points = self.__insertPointIntoIsoline(point, isoline_points)

        theta, phi = self.classic_grid.getAnglesOfPoint(point)

        mean_len = np.linalg.norm(self.mean)
        mean_normed = [x / mean_len for x in self.mean]
        theta_mean, phi_mean = self.classic_grid.getAnglesOfPoint(mean_normed)  # TODO change mean with actual pdf's maximum
        if mean_normed[0] == 0 and mean_normed[1] == 0:  # if mean vector points to one of the poles
            theta_mean = 0

        # TODO add comments to this hard logic
        isoline_points_angles = [self.classic_grid.getAnglesOfPoint(iso_point) for iso_point in isoline_points]  # TODO rename
        indices_inside_zero_phi_mean = [i for i in range(len(isoline_points_angles)) if 0 <= isoline_points_angles[i][1] <= phi_mean]
        indices_inside_zero_pi_minus_phi_mean = [i for i in range(len(isoline_points_angles)) if 0 <= isoline_points_angles[i][1] <= (math.pi - phi_mean)]

        diffs_to_theta_mean = [math.fabs(isoline_points_angles[i][0] - theta_mean) for i in indices_inside_zero_phi_mean]
        diffs_to_theta_plus_pi_mean = [math.fabs(isoline_points_angles[i][0] - ((theta_mean + math.pi) % (2 * math.pi))) for i in indices_inside_zero_pi_minus_phi_mean]

        min_diffs_to_theta_mean = np.array(diffs_to_theta_mean).argsort()
        min_diffs_to_theta_plus_pi_mean = np.array(diffs_to_theta_plus_pi_mean).argsort()

        intersection_points = []
        if len(min_diffs_to_theta_mean) > 0:
            intersection_points.append(isoline_points_angles[indices_inside_zero_phi_mean[min_diffs_to_theta_mean[0]]])
            intersection_points.append(isoline_points_angles[indices_inside_zero_phi_mean[min_diffs_to_theta_mean[1]]])  # TODO what if only 1 point?
        if len(min_diffs_to_theta_plus_pi_mean) > 0:
            if len(intersection_points) > 0:
                if diffs_to_theta_mean[min_diffs_to_theta_mean[0]] > diffs_to_theta_plus_pi_mean[min_diffs_to_theta_plus_pi_mean[0]]:
                    intersection_points[0] = (isoline_points_angles[indices_inside_zero_pi_minus_phi_mean[min_diffs_to_theta_plus_pi_mean[0]]])
                    intersection_points[1] = (isoline_points_angles[indices_inside_zero_pi_minus_phi_mean[min_diffs_to_theta_plus_pi_mean[1]]])
            else:
                intersection_points.append(isoline_points_angles[indices_inside_zero_pi_minus_phi_mean[min_diffs_to_theta_plus_pi_mean[0]]])
                intersection_points.append(isoline_points_angles[indices_inside_zero_pi_minus_phi_mean[min_diffs_to_theta_plus_pi_mean[1]]])

        intersection_point = np.mean(intersection_points, axis=0)
        intersection_point_coords = self.classic_grid.getCoordsOfPoint(intersection_point)

        intersection_point_index, isoline_points = self.__insertPointIntoIsoline(intersection_point_coords, isoline_points)

        return {"isoline": isoline_points, "points": [intersection_point_coords]}


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
                s.append(self.__calculateSStats(point, value))

            print("StatsAPI: Stats calculation time:", time.time() - start)

        except Exception as e:
            self.finish({"code": 500, "body": str(e)})
            raise e  # TODO remove (debug)
            return

        self.finish({"code": 200, "t": t, "s": s})

    # Something bad here, but it's fast workaround right now (CORS policy)
    def options(self):
        self.set_status(200)
        self.set_header("Access-Control-Allow-Origin", "*")
        self.set_header("Access-Control-Allow-Headers", "*")
        if not self.request.body:
            return
