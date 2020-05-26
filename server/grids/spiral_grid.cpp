#include "spiral_grid.hpp"

#include <cmath>
#include <future>
#include <vector>

#include "distributions/angular_gauss.hpp" // AngularGauss
#include "util/cfg.hpp"                    // cfg::kConfig
#include "util/util.hpp"                   // Bounds()

SpiralGrid::SpiralGrid(int N, AngularGauss* distr) : points_amount(N), distr(distr) {
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

double SpiralGrid::EvaluateFuncRoutine(int lower_bound, int upper_bound) {
    double local_max_value = 0;

    for (int i = lower_bound; i < upper_bound; ++i) {
        auto value = distr->Calc(points[i]);

        if (value > local_max_value)
            local_max_value = value;

        values[i] = value;
    }

    return local_max_value;
}

void SpiralGrid::EvaluateFunc() {
    values.resize(points_amount);

    int threads_amount = cfg::kConfig["threads"].get<int>();
    auto bounds = Bounds(threads_amount, points_amount);
    auto futures = new std::future<double>[threads_amount-1];

    for (int i = 0; i < threads_amount-1; ++i)
        futures[i] = std::async(&SpiralGrid::EvaluateFuncRoutine, this, bounds[i], bounds[i+1]);

    max_value = EvaluateFuncRoutine(bounds[threads_amount-1], bounds[threads_amount]);

    for (int i = 0; i < threads_amount-1; ++i) {
        double max_value_from_routine = futures[i].get();
        if (max_value_from_routine > max_value)
            max_value = max_value_from_routine;
    }
    delete[] futures;
}

double SpiralGrid::CalcIntegralInsideIsoline(double isoline_value) {
    // elementary_part_area is an area of an elementary part on a sphere under
    // every grid value
    const double elementary_part_area = 4 * M_PI / points_amount;

    double sum = 0;
    for (const auto it : values)
        if (it > isoline_value)
            sum += it;

    return sum * elementary_part_area;
}

const std::vector<double>& SpiralGrid::Values() {
    return values;
}

double SpiralGrid::MaxValue() {
    return max_value;
}
