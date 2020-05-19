#ifndef SPIRAL_GRID_HPP
#define SPIRAL_GRID_HPP

#include <array>
#include <unordered_map>
#include <vector>

#include <boost/container_hash/hash.hpp>

#include "distributions/angular_gauss.hpp"

class SpiralGrid {
    int points_amount;
    double max_value; // useful only for functions on grid that have positive values
    //double (*function_on_grid)(std::array<double, 3>&);

    // points is the grid in Cartesian coordinates.
    std::vector<std::array<double, 3>> points;
    // values is an evaluated function on grid.
    std::unordered_map<std::array<double, 3>, double, boost::hash<std::array<double, 3>>> values;
    // data is a vector from flattened values map.
    std::vector<double> data;

    void GenerateGrid();
    void EvaluateFunc();

public:
    AngularGauss* distr;

    SpiralGrid(int, AngularGauss*);
    double CalcIntegralInsideIsoline(double isoline_value);
    const std::vector<double>& Data();
    double MaxValue();
};

#endif // SPIRAL_GRID_HPP
