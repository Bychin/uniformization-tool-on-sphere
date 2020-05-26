#include "isoline.hpp"

#include <array>
#include <chrono>
#include <cmath>
#include <exception>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <vector>

#include <grids/classic_grid.hpp>
#include <grids/spiral_grid.hpp>
#include "types/types.hpp"

const int kBisectionMethodMaxIter = 50;

IsolineAPI::IsolineAPI(std::vector<double>& ratios, ClassicGrid* classic_grid, SpiralGrid* spiral_grid)
    : ratios(ratios), classic_grid(classic_grid), spiral_grid(spiral_grid) {
}

double IsolineAPI::GetIsolineValueByRatio(double ratio) const {
    const double eps = 0.0001;

    double f1 = 0;
    double f2 = spiral_grid->MaxValue();
    double isoline_value = 0;
    double integral = std::numeric_limits<double>::max();

    for (int i = 0; i < kBisectionMethodMaxIter && std::fabs(integral - ratio) > eps; ++i) {
        isoline_value = (f1 + f2) / 2;
        integral = spiral_grid->CalcIntegralInsideIsoline(isoline_value);

        if (integral > ratio)
            f1 = isoline_value;
        else
            f2 = isoline_value;
    }

    return isoline_value;
}

void IsolineAPI::Validate(void) const {
    if (classic_grid == nullptr) {
        throw std::invalid_argument("invalid IsolineAPI: classic_grid is nullptr");
    }
    if (spiral_grid == nullptr) {
        throw std::invalid_argument("invalid IsolineAPI: spiral_grid is nullptr");
    }
    if (ratios.size() == 0) {
        throw std::invalid_argument("invalid IsolineAPI: empty ratios");
    }
}

std::unordered_map<double, IsolineCoords> IsolineAPI::GetIsolines(void) {
    std::unordered_map<double, IsolineCoords> isolines;
    isolines.reserve(ratios.size());

    for (const auto r : ratios) {
        auto value  = GetIsolineValueByRatio(r);
        std::cout << "debug: for ratio=" << r << " got value=" << value << ", start getting isoline coords..." << std::endl;
        auto start = std::chrono::high_resolution_clock::now();
        auto points = classic_grid->GetIsolineCoords(value);
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start);
        std::cout << "debug: for value=" << value << " got " << points.size() << " isoline coords in " << duration.count() << "ms" << std::endl;
        isolines[r] = points;
    }

    return isolines;
}
