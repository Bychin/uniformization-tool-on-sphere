#ifndef ISOLINE_HPP
#define ISOLINE_HPP

#include <array>
#include <unordered_map>
#include <vector>

#include <grids/classic_grid.hpp>
#include <grids/spiral_grid.hpp>

class IsolineAPI {
    std::vector<double> ratios;

    ClassicGrid* classic_grid;
    SpiralGrid* spiral_grid;

    double GetIsolineValueByRatio(double ratio);

public:
    IsolineAPI(std::vector<double>& ratios,
               ClassicGrid* classic_grid,
               SpiralGrid* spiral_grid);
    // TODO bool Validate();
    std::unordered_map<double, std::vector<std::array<double, 3>>> GetIsolines();
};

#endif // ISOLINE_HPP
