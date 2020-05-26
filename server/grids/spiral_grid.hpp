#ifndef SPIRAL_GRID_HPP
#define SPIRAL_GRID_HPP

#include <array>
#include <unordered_map>
#include <vector>

#include <boost/container_hash/hash.hpp>

#include "distributions/angular_gauss.hpp"

class SpiralGrid {
    int points_amount;
    double max_value;

    // points is the grid in Cartesian coordinates.
    std::vector<std::array<double, 3>> points;

    // values is an evaluated function on grid, value in i'th position
    // corresponds to the i'th point in points vector
    std::vector<double> values;

    void GenerateGrid();
    void EvaluateFunc();
    double EvaluateFuncRoutine(int lower_bound, int upper_bound);

public:
    AngularGauss* distr;

    SpiralGrid(int, AngularGauss*);
    double CalcIntegralInsideIsoline(double isoline_value);
    const std::vector<double>& Values();
    double MaxValue();
};

#endif // SPIRAL_GRID_HPP
