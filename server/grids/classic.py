import math
import numpy as np


class ClassicGrid:
    def __init__(self, div, function):
        self.div = div
        self.function = function
        self.__generateAndCalculateGrid(function)

    def __generateAndCalculateGrid(self, function):
        self.points = {}  # (i, j) pair to point's (x, y, z)
        self.values = {}

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
                self.values[point] = function(point)

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

    # getIsolineCoords returns points' coordinates for isoline with value 'c'.
    # It uses marching squares algorithm (https://en.wikipedia.org/wiki/Marching_squares)
    def getIsolineCoords(self, c):
        points = []

        initial_trapezium_indices = None
        end_side = None
        end_point = None

        for i in range(len(self.trapeziums)):
            for j in range(len(self.trapeziums[i])):
                trapezium = self.trapeziums[i][j]
                index = self.__getTrapeziumIndex(trapezium, c)

                if index not in [0, 15, 5, 10]:
                    end_side, end_point = self.__processSegment(trapezium, index, -1, c)
                    points.append(end_point)
                    initial_trapezium_indices = (i, j)
                    break

        #if end_side == 2:
        #    print("end side 2")
        #    print(end_point, initial_trapezium_indices, c, len(self.trapeziums), len(self.trapeziums[initial_trapezium_indices[1]]))

        side, trapezium_indices = self.__getNextTrapeziumIndices(end_side, initial_trapezium_indices)

        #print("here")

        while trapezium_indices != initial_trapezium_indices:
            #print(trapezium_indices)
            trapezium = self.trapeziums[trapezium_indices[0]][trapezium_indices[1]]
            index = self.__getTrapeziumIndex(trapezium, c)
            side, point = self.__processSegment(trapezium, index, side, c)
            points.append(point)
            side, trapezium_indices = self.__getNextTrapeziumIndices(side, trapezium_indices)

        return points


    def __getNextTrapeziumIndices(self, prev_side, prev_trapezium_indices):
        #print("__getNextTrapeziumIndices")
        if prev_side == 0:  # we must get the one above
            return 2, self.__getUpperTrapeziumIndices(prev_trapezium_indices)

        elif prev_side == 1:  # get the right one
            return 3, self.__getRightTrapeziumIndices(prev_trapezium_indices)

        elif prev_side == 2:  # get the lower one
            return 0, self.__getLowerTrapeziumIndices(prev_trapezium_indices)

        elif prev_side == 3:  # get the left one
            return 1, self.__getLeftTrapeziumIndices(prev_trapezium_indices)

    def __getUpperTrapeziumIndices(self, prev_trapezium_indices):
        i, j = prev_trapezium_indices
        assert j != 0
        return i, j - 1

    def __getRightTrapeziumIndices(self, prev_trapezium_indices):
        i, j = prev_trapezium_indices
        if i == len(self.trapeziums) - 1:
            i = -1
        return i + 1, j

    def __getLowerTrapeziumIndices(self, prev_trapezium_indices):
        i, j = prev_trapezium_indices
        assert j != len(self.trapeziums[i]) - 1
        return i, j + 1

    def __getLeftTrapeziumIndices(self, prev_trapezium_indices):
        i, j = prev_trapezium_indices
        if i == 0:
            i = len(self.trapeziums)
        return i - 1, j

    # Linear interpolation to find intersection point with value c between vertices with values a and b
    # where x1 = f^(-1)(a); x2 = f^(-1)(b); sign((a - c) * (b - c)) should be -1.
    # Note that x1 and x2 - 3-D coordinates.
    def __findIntersection(self, x1, x2, a, b, c):
        k = (c - a) / (b - a)
        if k > 1.:
            print(k, c, a, b)
            assert k < 1. # TODO
        
        x_int = (1 - k) * x1[0] + k * x2[0]
        y_int = (1 - k) * x1[1] + k * x2[1]
        z_int = (1 - k) * x1[2] + k * x2[2]

        return (x_int, y_int, z_int)

    # start_side is 0 - AD, 1 - DC, 2 - CB, 3 - BA
    # end_side, end_point are returned
    def __processSegment(self, trapezium, index, start_side, const):

        # 0 -- 1    1 -- 0
        # |    | or |    |
        # 1 -- 1    0 -- 0
        if index == 7 or index == 8:
            if start_side == 0:
                end_side = 3
                # calc A and B
                point_1 = self.points[trapezium[0]]
                point_2 = self.points[trapezium[1]]

                end_point = self.__findIntersection(point_1, point_2, self.values[point_1], self.values[point_2], const)
                return end_side, end_point
            else:
                end_side = 0
                # calc A and D
                point_1 = self.points[trapezium[0]]
                point_2 = self.points[trapezium[3]]

                end_point = self.__findIntersection(point_1, point_2, self.values[point_1], self.values[point_2], const)
                return end_side, end_point

        # 0 -- 1    1 -- 0
        # |    | or |    |
        # 0 -- 0    1 -- 1
        elif index == 4 or index == 11:
            if start_side == 0:
                end_side = 1
                # calc C and D
                point_1 = self.points[trapezium[2]]
                point_2 = self.points[trapezium[3]]

                end_point = self.__findIntersection(point_1, point_2, self.values[point_1], self.values[point_2], const)
                return end_side, end_point
            else:
                end_side = 0
                # calc A and D
                point_1 = self.points[trapezium[0]]
                point_2 = self.points[trapezium[3]]

                end_point = self.__findIntersection(point_1, point_2, self.values[point_1], self.values[point_2], const)
                return end_side, end_point

        # 0 -- 0    1 -- 1
        # |    | or |    |
        # 0 -- 1    1 -- 0
        elif index == 2 or index == 13:
            if start_side == 1:
                end_side = 2
                # calc B and C
                point_1 = self.points[trapezium[1]]
                point_2 = self.points[trapezium[2]]

                end_point = self.__findIntersection(point_1, point_2, self.values[point_1], self.values[point_2], const)
                return end_side, end_point
            else:
                end_side = 1
                # calc C and D
                point_1 = self.points[trapezium[2]]
                point_2 = self.points[trapezium[3]]

                end_point = self.__findIntersection(point_1, point_2, self.values[point_1], self.values[point_2], const)
                return end_side, end_point

        # 0 -- 0    1 -- 1
        # |    | or |    |
        # 1 -- 0    0 -- 1
        elif index == 1 or index == 14:
            if start_side == 2:
                end_side = 3
                # calc A and B
                point_1 = self.points[trapezium[0]]
                point_2 = self.points[trapezium[1]]

                end_point = self.__findIntersection(point_1, point_2, self.values[point_1], self.values[point_2], const)
                return end_side, end_point
            else:
                end_side = 2
                # calc B and C
                point_1 = self.points[trapezium[1]]
                point_2 = self.points[trapezium[2]]

                end_point = self.__findIntersection(point_1, point_2, self.values[point_1], self.values[point_2], const)
                return end_side, end_point

        # 0 -- 1    1 -- 0
        # |    | or |    |
        # 0 -- 1    1 -- 0
        elif index == 6 or index == 9:
            if start_side == 0:
                end_side = 2
                # calc B and C
                point_1 = self.points[trapezium[1]]
                point_2 = self.points[trapezium[2]]

                end_point = self.__findIntersection(point_1, point_2, self.values[point_1], self.values[point_2], const)
                return end_side, end_point
            else:
                end_side = 0
                # calc A and D
                point_1 = self.points[trapezium[0]]
                point_2 = self.points[trapezium[3]]

                end_point = self.__findIntersection(point_1, point_2, self.values[point_1], self.values[point_2], const)
                return end_side, end_point

        # 0 -- 0    1 -- 1
        # |    | or |    |
        # 1 -- 1    0 -- 0
        elif index == 3 or index == 12:
            if start_side == 1:
                end_side = 3
                # calc A and B
                point_1 = self.points[trapezium[0]]
                point_2 = self.points[trapezium[1]]

                end_point = self.__findIntersection(point_1, point_2, self.values[point_1], self.values[point_2], const)
                return end_side, end_point
            else:
                end_side = 1
                # calc C and D
                point_1 = self.points[trapezium[2]]
                point_2 = self.points[trapezium[3]]

                end_point = self.__findIntersection(point_1, point_2, self.values[point_1], self.values[point_2], const)
                return end_side, end_point

        # 0 -- 1
        # |    |
        # 1 -- 0
        elif index == 5:  # TODO check first and change index
            if start_side == 0:
                end_side = 3
                # calc A and B
                point_1 = self.points[trapezium[0]]
                point_2 = self.points[trapezium[1]]

                end_point = self.__findIntersection(point_1, point_2, self.values[point_1], self.values[point_2], const)
                return end_side, end_point
            elif start_side == 3:
                end_side = 0
                # calc A and D
                point_1 = self.points[trapezium[0]]
                point_2 = self.points[trapezium[3]]

                end_point = self.__findIntersection(point_1, point_2, self.values[point_1], self.values[point_2], const)
                return end_side, end_point
            elif start_side == 1:
                end_side = 2
                # calc B and C
                point_1 = self.points[trapezium[1]]
                point_2 = self.points[trapezium[2]]

                end_point = self.__findIntersection(point_1, point_2, self.values[point_1], self.values[point_2], const)
                return end_side, end_point
            else:  # start_side == 2
                end_side = 1
                # calc C and D
                point_1 = self.points[trapezium[2]]
                point_2 = self.points[trapezium[3]]

                end_point = self.__findIntersection(point_1, point_2, self.values[point_1], self.values[point_2], const)
                return end_side, end_point

        # 0 -- 1
        # |    |
        # 1 -- 0
        elif index == 10:
            if start_side == 0:
                end_side = 1
                # calc C and D
                point_1 = self.points[trapezium[2]]
                point_2 = self.points[trapezium[3]]

                end_point = self.__findIntersection(point_1, point_2, self.values[point_1], self.values[point_2], const)
                return end_side, end_point
            elif start_side == 1:
                end_side = 0
                # calc A and D
                point_1 = self.points[trapezium[0]]
                point_2 = self.points[trapezium[3]]

                end_point = self.__findIntersection(point_1, point_2, self.values[point_1], self.values[point_2], const)
                return end_side, end_point
            elif start_side == 2:
                end_side = 3
                # calc A and B
                point_1 = self.points[trapezium[0]]
                point_2 = self.points[trapezium[1]]

                end_point = self.__findIntersection(point_1, point_2, self.values[point_1], self.values[point_2], const)
                return end_side, end_point
            else:  # start_side == 3
                end_side = 2
                # calc B and C
                point_1 = self.points[trapezium[1]]
                point_2 = self.points[trapezium[2]]

                end_point = self.__findIntersection(point_1, point_2, self.values[point_1], self.values[point_2], const)
                return end_side, end_point

        else:
            print("WARNING: strange segment:", trapezium, index, start_side, const)
            raise False

    def __getTrapeziumIndex(self, trapezium, const):
        index = 0

        # in a clockwise direction
        if (self.values[self.points[trapezium[0]]] > const):  # point A > const
            index += 8  # int('1000', 2)
        if (self.values[self.points[trapezium[3]]] > const):  # point D > const
            index += 4  # int('0100', 2)
        if (self.values[self.points[trapezium[2]]] > const):  # point C > const
            index += 2  # int('0010', 2)
        if (self.values[self.points[trapezium[1]]] > const):  # point B > const
            index += 1  # int('0001', 2)

        return index
