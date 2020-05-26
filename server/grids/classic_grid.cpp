#include "classic_grid.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <string>
#include <thread>
#include <unordered_set>
#include <vector>

#include "distributions/angular_gauss.hpp" // AngularGauss
#include "types/types.hpp"
#include "util/cfg.hpp"                    // cfg::kConfig
#include "util/util.hpp"                   // Bounds()

ClassicGrid::ClassicGrid(int grid_div, AngularGauss* distr) : div(grid_div), distr(distr) {
    GenerateGridAndEvaluateFunc();
}

const AngularGauss* ClassicGrid::Func(void) const {
    return distr;
}

void ClassicGrid::EvaluateFuncRoutine(int lower_bound, int upper_bound) {
    for (int i = lower_bound; i < upper_bound; ++i)
        values[i] = distr->Calc(cartesian_points[i]);
}

void ClassicGrid::GenerateGridAndEvaluateFunc(void) {
    double angle = M_PI / div;

    points.reserve(div * (div + 1));
    cartesian_points.reserve(div * (div + 1));

    for (int i = 0; i < div; ++i) { // i goes around the sphere
        double theta = 2 * i * angle;
        double si    = std::sin(theta);
        double ci    = std::cos(theta);

        for (int j = 0; j <= div; ++j) { // j goes from N to S of the unit sphere
            double phi = j * angle;
            double sj  = std::sin(phi);
            double cj  = std::cos(phi);

            CoordsOfPoint point = {ci * sj, si * sj, cj}; // X, Y, Z
            cartesian_points.push_back(point);
            points[{i, j}] = cartesian_points.size()-1;
        }
    }

    values.reserve(div * (div + 1));

    int threads_amount = cfg::kConfig["threads"].get<int>();
    auto bounds = Bounds(threads_amount, cartesian_points.size());
    auto threads = std::vector<std::thread>(threads_amount-1);

    for (int i = 0; i < threads_amount-1; ++i)
        threads[i] = std::thread(&ClassicGrid::EvaluateFuncRoutine, this, bounds[i], bounds[i+1]);

    EvaluateFuncRoutine(bounds[threads_amount-1], bounds[threads_amount]);

    for (int i = 0; i < threads_amount-1; ++i)
        threads[i].join();

    /* trapezium is stored as 4 vertices from A to D anticlockwise:
     *   A -<- D          A (== D)
     *  /       \   or   / \
     * B --->--- C      B - C
     */
    trapeziums.reserve(div);

    for (int i = 0; i < div; ++i) {
        std::vector<Trapezium> trapezium_column;
        trapezium_column.reserve(div);

        for (int j = 0; j < div; ++j) {
            int right_i = i < div - 1 ? i + 1 : 0;

            ClassicGridPoint A = {i, j};
            ClassicGridPoint B = {i, j + 1};
            ClassicGridPoint C = {right_i, j + 1};
            ClassicGridPoint D = {right_i, j};

            trapezium_column.push_back({A, B, C, D});
        }

        trapeziums.push_back(trapezium_column);
    }
}

CoordsOfPoint ClassicGrid::GetCoordsOfPoint(const AnglesOfPoint& point) {
    double theta = point[0], phi = point[1];

    double si = std::sin(theta);
    double ci = std::cos(theta);

    double sj = std::sin(phi);
    double cj = std::cos(phi);

    return {ci * sj, si * sj, cj}; // X, Y, Z
}

AnglesOfPoint ClassicGrid::GetAnglesOfPoint(const CoordsOfPoint& point) {
    double phi = std::atan2(point[1], point[0]); // arctg(y/x)
    if (phi < 0) {
        phi += 2 * M_PI; // due to our classic grid construction logic
    }

    double theta = std::acos(point[2]); // arccos(z)

    return {phi, theta};
}

std::vector<CoordsOfPoint> ClassicGrid::GetIsolineCoords(double iso_value) {
    const size_t initial_isoline_points_amount = 1000;
    std::vector<CoordsOfPoint> isoline_points;
    isoline_points.reserve(initial_isoline_points_amount);

    int initial_trapezium_side;
    TrapeziumIndex initial_trapezium_index;

    bool starting_trapezium_found = false;
    const std::unordered_set<int> marching_indicies_of_bad_to_start_trapeziums = {0, 15, 5, 10};

    for (int i = 0; i < trapeziums.size() && !starting_trapezium_found; ++i) {
        for (int j = 0; j < trapeziums[i].size(); ++j) {
            auto trapezium = trapeziums[i][j];
            auto marching_index = GetMarchingSquareIndex(trapezium, iso_value);

            if (marching_indicies_of_bad_to_start_trapeziums.find(marching_index) == marching_indicies_of_bad_to_start_trapeziums.end()) {
                CoordsOfPoint last_point;
                std::tie(initial_trapezium_side, last_point) = ProcessSegment(trapezium, marching_index, -1, iso_value);
                isoline_points.push_back(last_point);
                initial_trapezium_index = {i, j};

                starting_trapezium_found = true;
                break;
            }
        }
    }

    int side;
    TrapeziumIndex trapezium_index;
    std::tie(side, trapezium_index) = GetNextTrapeziumIndex(initial_trapezium_side, initial_trapezium_index);

    while (trapezium_index != initial_trapezium_index) {
        auto trapezium = trapeziums[trapezium_index[0]][trapezium_index[1]];
        auto marching_index = GetMarchingSquareIndex(trapezium, iso_value);

        CoordsOfPoint point;
        std::tie(side, point) = ProcessSegment(trapezium, marching_index, side, iso_value);
        isoline_points.push_back(point);
        std::tie(side, trapezium_index) = GetNextTrapeziumIndex(side, trapezium_index);
    }

    return isoline_points;
}

std::tuple<int, TrapeziumIndex> ClassicGrid::GetNextTrapeziumIndex(int previous_side, const TrapeziumIndex& previous) const {
    switch (previous_side) {
    case 0: // we must get the upper one
        return std::make_tuple(2, GetUpperTrapeziumIndex(previous));

    case 1: // get the right one
        return std::make_tuple(3, GetRightTrapeziumIndex(previous));

    case 2: // get the lower one
        return std::make_tuple(0, GetLowerTrapeziumIndex(previous));

    case 3: // get the left one
        return std::make_tuple(1, GetLeftTrapeziumIndex(previous));

    default:
        throw std::invalid_argument("GetNextTrapeziumIndices: got invalid prev_side argument: " + std::to_string(previous_side));
    }
}

TrapeziumIndex ClassicGrid::GetUpperTrapeziumIndex(const TrapeziumIndex& previous) const {
    int i = previous[0], j = previous[1];
    assert(j != 0);
    return {i, j - 1};
}

TrapeziumIndex ClassicGrid::GetLowerTrapeziumIndex(const TrapeziumIndex& previous) const {
    int i = previous[0], j = previous[1];
    assert(j != trapeziums[i].size() - 1);
    return {i, j + 1};
}

TrapeziumIndex ClassicGrid::GetLeftTrapeziumIndex(const TrapeziumIndex& previous) const {
    int i = previous[0], j = previous[1];
    if (i == 0) {
        i = trapeziums.size();
    }
    return {i - 1, j};
}

TrapeziumIndex ClassicGrid::GetRightTrapeziumIndex(const TrapeziumIndex& previous) const {
    int i = previous[0], j = previous[1];
    if (i == trapeziums.size() - 1) {
        i = -1;
    }
    return {i + 1, j};
}

CoordsOfPoint ClassicGrid::FindIntersection(const CoordsOfPoint& x1, const CoordsOfPoint& x2, double a, double b, double c) const {
    double k = (c - a) / (b - a);
    if (k >= 1.) {
        std::string err = "FindIntersection, got invalid k: k="+std::to_string(k)+", x1=["+std::to_string(x1[0])+","+std::to_string(x1[1])+","+std::to_string(x1[2])+", x2="+
            std::to_string(x2[0])+","+std::to_string(x2[1])+","+std::to_string(x2[2])+", a="+std::to_string(a)+", b="+std::to_string(b)+", c="+std::to_string(c);
        throw std::invalid_argument(err);
    }

    // coordinates of intersection point
    double x = (1 - k) * x1[0] + k * x2[0];
    double y = (1 - k) * x1[1] + k * x2[1];
    double z = (1 - k) * x1[2] + k * x2[2];

    return {x, y, z};
}

std::tuple<int, CoordsOfPoint> ClassicGrid::ProcessSegment(const Trapezium& trapezium, int index, int start_side, double iso_value) {
    // start_side must be one of: 0 - AD, 1 - DC, 2 - CB, 3 - BA
    // (see trapezium scheme in GenerateGridAndEvaluateFunc())

    int end_side = -1;
    int point_1_index, point_2_index;

    // 0 -- 1    1 -- 0
    // |    | or |    |
    // 1 -- 1    0 -- 0
    if (index == 7 || index == 8) {
        if (start_side == 0) {
            end_side = 3;
            // calc A and B
            point_1_index = points[trapezium[0]];
            point_2_index = points[trapezium[1]];
        } else {
            end_side = 0;
            // calc A and D
            point_1_index = points[trapezium[0]];
            point_2_index = points[trapezium[3]];
        }
    }

    // 0 -- 1    1 -- 0
    // |    | or |    |
    // 0 -- 0    1 -- 1
    else if (index == 4 || index == 11) {
        if (start_side == 0) {
            end_side = 1;
            // calc C and D
            point_1_index = points[trapezium[2]];
            point_2_index = points[trapezium[3]];
        } else {
            end_side = 0;
            // calc A and D
            point_1_index = points[trapezium[0]];
            point_2_index = points[trapezium[3]];
        }
    }

    // 0 -- 0    1 -- 1
    // |    | or |    |
    // 0 -- 1    1 -- 0
    else if (index == 2 || index == 13) {
        if (start_side == 1) {
            end_side = 2;
            // calc B and C
            point_1_index = points[trapezium[1]];
            point_2_index = points[trapezium[2]];
        } else {
            end_side = 1;
            // calc C and D
            point_1_index = points[trapezium[2]];
            point_2_index = points[trapezium[3]];
        }
    }

    // 0 -- 0    1 -- 1
    // |    | or |    |
    // 1 -- 0    0 -- 1
    else if (index == 1 || index == 14) {
        if (start_side == 2) {
            end_side = 3;
            // calc A and B
            point_1_index = points[trapezium[0]];
            point_2_index = points[trapezium[1]];
        } else {
            end_side = 2;
            // calc B and C
            point_1_index = points[trapezium[1]];
            point_2_index = points[trapezium[2]];
        }
    }

    // 0 -- 1    1 -- 0
    // |    | or |    |
    // 0 -- 1    1 -- 0
    else if (index == 6 || index == 9) {
        if (start_side == 0) {
            end_side = 2;
            // calc B and C
            point_1_index = points[trapezium[1]];
            point_2_index = points[trapezium[2]];
        } else {
            end_side = 0;
            // calc A and D
            point_1_index = points[trapezium[0]];
            point_2_index = points[trapezium[3]];
        }
    }

    // 0 -- 0    1 -- 1
    // |    | or |    |
    // 1 -- 1    0 -- 0
    else if (index == 3 || index == 12) {
        if (start_side == 1) {
            end_side = 3;
            // calc A and B
            point_1_index = points[trapezium[0]];
            point_2_index = points[trapezium[1]];
        } else {
            end_side = 1;
            // calc C and D
            point_1_index = points[trapezium[2]];
            point_2_index = points[trapezium[3]];
        }
    }

    // 0 -- 1
    // |    |
    // 1 -- 0
    // TODO think about it one more time (going to the wrong side second time?)
    else if (index == 5) { // TODO check first and change index
        if (start_side == 0) {
            end_side = 3;
            // calc A and B
            point_1_index = points[trapezium[0]];
            point_2_index = points[trapezium[1]];
        } else if (start_side == 3) {
            end_side = 0;
            // calc A and D
            point_1_index = points[trapezium[0]];
            point_2_index = points[trapezium[3]];
        } else if (start_side == 1) {
            end_side = 2;
            // calc B and C
            point_1_index = points[trapezium[1]];
            point_2_index = points[trapezium[2]];
        } else { // start_side == 2
            end_side = 1;
            // calc C and D
            point_1_index = points[trapezium[2]];
            point_2_index = points[trapezium[3]];
        }
    }

    // 0 -- 1
    // |    |
    // 1 -- 0
    else if (index == 10) {
        if (start_side == 0) {
            end_side = 1;
            // calc C and D
            point_1_index = points[trapezium[2]];
            point_2_index = points[trapezium[3]];
        } else if (start_side == 1) {
            end_side = 0;
            // calc A and D
            point_1_index = points[trapezium[0]];
            point_2_index = points[trapezium[3]];
        } else if (start_side == 2) {
            end_side = 3;
            // calc A and B
            point_1_index = points[trapezium[0]];
            point_2_index = points[trapezium[1]];
        } else { // start_side == 3
            end_side = 2;
            // calc B and C
            point_1_index = points[trapezium[1]];
            point_2_index = points[trapezium[2]];
        }
    }

    else {
        std::string err = "ProcessSegment: got invalid index/start_side, index="+std::to_string(index)+", start_side="+std::to_string(start_side)+", iso_value="+std::to_string(start_side);
        throw std::invalid_argument(err);
    }

    auto point_1 = cartesian_points[point_1_index];
    auto point_2 = cartesian_points[point_2_index];

    auto end_point = FindIntersection(point_1, point_2, values[point_1_index], values[point_2_index], iso_value);
    return std::make_tuple(end_side, end_point);
}

int ClassicGrid::GetMarchingSquareIndex(const Trapezium& trapezium, double value) {
    int index = 0;

    // in a clockwise direction
    if (values[points[trapezium[0]]] > value) // value of point A > value
        index += 8;                           // 0b1000
    if (values[points[trapezium[3]]] > value) // value of point D > value
        index += 4;                           // 0b0100
    if (values[points[trapezium[2]]] > value) // value of point C > value
        index += 2;                           // 0b0010
    if (values[points[trapezium[1]]] > value) // value of point B > value
        index += 1;                           // 0b0001

    return index;
}
