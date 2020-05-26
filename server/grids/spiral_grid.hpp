#ifndef SPIRAL_GRID_HPP
#define SPIRAL_GRID_HPP

#include <array>
#include <unordered_map>
#include <vector>

#include <boost/container_hash/hash.hpp>

#include "distributions/angular_gauss.hpp"
#include "types/types.hpp"

class SpiralGrid {
    int points_amount;
    double max_value;
    AngularGauss* distr;

    // points is the grid in Cartesian coordinates.
    std::vector<CoordsOfPoint> points;

    // values is an evaluated function on grid, value in i'th position
    // corresponds to the i'th point in points vector
    std::vector<double> values;

    void GenerateGrid(void);
    void EvaluateFunc(void);
    double EvaluateFuncRoutine(int lower_bound, int upper_bound);

public:
    // SpiralGrid creates new spiral grid with N points and function distr on it
    SpiralGrid(int N, AngularGauss* disrt);

    // Func returns function on grid
    const AngularGauss* Func(void) const;
    // MaxValue returns the maximum value of the evaluated function the spiral grid
    double MaxValue(void) const;
    // Values returns vector of evaluated function's values on spiral grid
    const std::vector<double>& Values(void) const;

    // CalcIntegralInsideIsoline calculates surface integral of function on
    // spiral grid
    double CalcIntegralInsideIsoline(double isoline_value);
};

#endif // SPIRAL_GRID_HPP
