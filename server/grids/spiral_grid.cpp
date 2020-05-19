#include "spiral_grid.hpp"

#include <cmath>
#include <vector>

#include "distributions/angular_gauss.hpp"

SpiralGrid::SpiralGrid(int N, AngularGauss* distr) : points_amount(N), distr(distr) {
    max_value = 0;
    GenerateGrid();
    EvaluateFunc();
}

void SpiralGrid::GenerateGrid() {
    points.reserve(points_amount);
    double theta = 0.;

    for (int k = 1; k <= points_amount; ++k) {
        double h = -1. + (2. * k - 2.) / (points_amount - 1.);
        double phi = std::acos(h);

        double s_phi = std::sin(phi);
        double c_phi = std::cos(phi);

        if (k != 1 && k != points_amount)
            theta = std::fmod(theta + 3.8 / std::sqrt(points_amount * (1. - h * h)), 2. * M_PI);
        else
            theta = 0;

        double s_theta = std::sin(theta);
        double c_theta = std::cos(theta);

        double x = s_phi * c_theta;
        double y = s_phi * s_theta;
        double z = c_phi;

        points.push_back({x, y, z});
    }
}

// TODO parallel https://solarianprogrammer.com/2011/12/16/cpp-11-thread-tutorial/
void SpiralGrid::EvaluateFunc() {
    values.reserve(points_amount); // actually, we don't need this map
    data.reserve(points_amount);

    for (const auto& it : points) {
        auto value = distr->Calc(it);
        if (value > max_value)
            max_value = value;

        values[it] = value;
        data.push_back(value);
    }
}

double SpiralGrid::CalcIntegralInsideIsoline(double isoline_value) {
    // elementary_part_area is an area of an elementary part on a sphere under
    // every grid value
    double elementary_part_area = 4 * M_PI / points_amount;
    double sum = 0;

    for (const auto it : data)
        if (it > isoline_value)
            sum += it;

    return sum * elementary_part_area;
}

const std::vector<double>& SpiralGrid::Data() {
    return data;
}

double SpiralGrid::MaxValue() {
    return max_value;
}
