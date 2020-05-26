#ifndef ISOLINE_HPP
#define ISOLINE_HPP

#include <array>
#include <unordered_map>
#include <vector>

#include <grids/classic_grid.hpp>
#include <grids/spiral_grid.hpp>
#include <types/types.hpp>

class IsolineAPI {
    // ratios is a vector of ratios for the corresponding isolines
    std::vector<double> ratios;

    ClassicGrid* classic_grid;
    SpiralGrid* spiral_grid;

    // GetIsolineValueByRatio uses simple bisection method to find the
    // corresponding function value on classic/spiral grid for a given isoline
    // ratio
    double GetIsolineValueByRatio(double ratio) const;

public:
    IsolineAPI(std::vector<double>& ratios, ClassicGrid* classic_grid, SpiralGrid* spiral_grid);

    void Validate(void) const;
    // GetIsolines returns a map of isoline values to their coordinates
    std::unordered_map<double, IsolineCoords> GetIsolines(void);
};

#endif // ISOLINE_HPP
