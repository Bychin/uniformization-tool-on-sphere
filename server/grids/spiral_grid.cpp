#include "spiral_grid.hpp"

#include <cmath>
#include <vector>

#include "distributions/angular_gauss.hpp"

SpiralGrid::SpiralGrid(int N, AngularGauss* distr) : points_amount(N), distr(distr) {
    //points_amount              = N;
    //function_on_grid = distribution_func;

    elementary_part_area = 4 * M_PI / N;

    GenerateGrid();
    EvaluateFunc();
}

void SpiralGrid::GenerateGrid() {
    points.reserve(points_amount);
    double theta = 0.; // TODO probably you don't need it

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

// TODO parallel
void SpiralGrid::EvaluateFunc() {
    values.reserve(points_amount);
    for (auto& it : points)
        values[it] = distr->Calc(it); //function_on_grid(it);

    data.reserve(points_amount);
    for (auto it = values.begin(); it != values.end(); ++it) {
        data.push_back(it->second);
    }
}

double SpiralGrid::CalcIntegralInsideIsoline(double isoline_value) {
    int counter = 0;
    double sum = 0;

    for (auto it : data) {
        if (it > isoline_value) {
            ++counter;
            sum += it;
        }
    }

    if (!counter) { // if isoline_value is pdf's maximum
        return 0;
    }

    return sum * elementary_part_area;
}

const std::vector<double>& SpiralGrid::Data() {
    return data;
}
