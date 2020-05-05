#ifndef SPIRAL_GRID_HPP
#define SPIRAL_GRID_HPP

#include <array>
#include <unordered_map>
#include <vector>

#include <boost/container_hash/hash.hpp>

#include "distributions/angular_gauss.hpp"

class SpiralGrid {
    int points_amount;
    AngularGauss* distr;
    //double (*function_on_grid)(std::array<double, 3>&);

    // elementary_part_area is an area of an elementary part of a sphere for
    // every grid value
    double elementary_part_area;

    // points is the grid in Cartesian coordinates.
    std::vector<std::array<double, 3>> points;
    // values is an evaluated function on grid.
    std::unordered_map<std::array<double, 3>, double, boost::hash<std::array<double, 3>>> values;
    // data is a vector from flattened values map.
    std::vector<double> data;

    void GenerateGrid();
    void EvaluateFunc();

public:
    SpiralGrid(int, AngularGauss*);
    double CalcIntegralInsideIsoline(double isoline_value);
    const std::vector<double>& Data();
};

#endif // SPIRAL_GRID_HPP
