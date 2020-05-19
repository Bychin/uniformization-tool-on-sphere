#ifndef CLASSIC_GRID_HPP
#define CLASSIC_GRID_HPP

#include <array>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <boost/container_hash/hash.hpp>

#include "distributions/angular_gauss.hpp"

class ClassicGrid {
    int div;
    //double (*function_on_grid)(std::array<double, 3>&);

    // TODO std::array<int, 2> as ClassicGridPoint? typedef
    // TODO same for std::array<double, 3>?
    // TODO std::array<std::array<int, 2>, 4> -> Trapezium?

    // points is a map of (i, j) pairs, that are indices of grid points, to
    // their (x, y, z) coordinates.
    std::unordered_map<std::array<int, 2>, std::array<double, 3>, boost::hash<std::array<int, 2>>> points;
    // values is an evaluated function on grid.
    std::unordered_map<std::array<double, 3>, double, boost::hash<std::array<double, 3>>> values;
    // trapeziums are the elementary objects of which the grid consists.
    std::vector<std::vector<std::array<std::array<int, 2>, 4>>> trapeziums;

    void GenerateGridAndEvaluateFunc();
    std::tuple<int, std::array<int, 2>> GetNextTrapeziumIndices(int prev_side, const std::array<int, 2>& prev_trapezium_indices);
    std::array<int, 2> GetUpperTrapeziumIndices(const std::array<int, 2>& prev_trapezium_indices);
    std::array<int, 2> GetLowerTrapeziumIndices(const std::array<int, 2>& prev_trapezium_indices);
    std::array<int, 2> GetLeftTrapeziumIndices(const std::array<int, 2>& prev_trapezium_indices);
    std::array<int, 2> GetRightTrapeziumIndices(const std::array<int, 2>& prev_trapezium_indices);
    std::array<double, 3> FindIntersection(
        const std::array<double, 3>& x1, const std::array<double, 3>& x2, double a, double b, double c);
    std::tuple<int, std::array<double, 3>> ProcessSegment(const std::array<std::array<int, 2>, 4>& trapezium,
                                                          int index,
                                                          int start_side,
                                                          double iso_value);
    int GetTrapeziumIndex(const std::array<std::array<int, 2>, 4>& trapezium, double value);

public:
    AngularGauss* distr;

    ClassicGrid(int, AngularGauss*);
    std::array<double, 3> GetCoordsOfPoint(const std::array<double, 2>& point);
    std::array<double, 2> GetAnglesOfPoint(const std::array<double, 3>& point);
    std::vector<std::array<double, 3>> GetIsolineCoords(double iso_value);
};

#endif // CLASSIC_GRID_HPP
