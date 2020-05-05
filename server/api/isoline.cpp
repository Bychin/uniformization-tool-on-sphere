#include "isoline.hpp"

#include <array>
#include <cmath>
#include <limits>
#include <vector>

#include <grids/classic_grid.hpp>
#include <grids/spiral_grid.hpp>

IsolineAPI::IsolineAPI(std::array<double, 3> mean,
                       std::array<double, 6> cov,
                       std::vector<double> ratios,
                       ClassicGrid* classic_grid,
                       SpiralGrid* spiral_grid)
    : mean(mean), cov(cov), ratios(ratios), classic_grid(classic_grid), spiral_grid(spiral_grid) {
}

double IsolineAPI::GetIsolineValueByRatio(double ratio) {
    double f1            = 0;
    double f2            = *std::max_element(spiral_grid->Data().begin(), spiral_grid->Data().end());
    double isoline_value = 0;

    double integral = std::numeric_limits<double>::max();
    double eps      = 0.001;

    for (int i = 0; i < 50 || std::fabs(integral - ratio) > eps; ++i) {
        isoline_value = (f1 + f2) / 2;
        integral      = spiral_grid->CalcIntegralInsideIsoline(isoline_value);

        if (integral > ratio)
            f1 = isoline_value;
        else
            f2 = isoline_value;
    }

    return isoline_value;
}

std::unordered_map<double, std::vector<std::array<double, 3>>> IsolineAPI::GetIsolines() {
    std::unordered_map<double, std::vector<std::array<double, 3>>> isolines;
    isolines.reserve(ratios.size());

    for (auto r : ratios) {
        auto value  = GetIsolineValueByRatio(r);
        auto points = classic_grid->GetIsolineCoords(value);
        isolines[r] = points;
    }

    return isolines;
}
