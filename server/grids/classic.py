import math
import numpy as np


class ClassicGrid:
    def __init__(self, div, function):
        self.div = div
        self.__generateAndCalculateGrid(function)

    def __generateAndCalculateGrid(self, function):
        self.points = {}  # (i, j) pair to point's (x, y, z)
        self.function = {}

        angle = math.pi / self.div
        for i in range(self.div):  # i goes around the sphere
            theta = 2 * i * angle
            si = math.sin(theta)
            ci = math.cos(theta)

            for j in range(self.div + 1):  # j goes from N to S of the unit sphere
                phi = j * angle
                sj = math.sin(phi)
                cj = math.cos(phi)

                point = (ci * sj, si * sj, cj)  # X, Y, Z
                pointIndex = (i, j)

                self.points[pointIndex] = point
                self.function[point] = function(point)

        # trapezium is stored as 4 vertices from A to D anticlockwise:
        #   A -<- D          A (== D)
        #  /       \   or   / \
        # B --->--- C      B - C

        self.trapeziums = []

        for i in range(self.div):
            column = []

            for j in range(self.div):
                right_i = i + 1 if i < self.div - 1 else 0

                A = (i, j)
                B = (i, j + 1)
                C = (right_i, j + 1)
                D = (right_i, j)

                column.append((A, B, C, D))

            self.trapeziums.append(column)
