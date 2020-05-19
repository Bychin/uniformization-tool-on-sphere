#include "classic_grid.hpp"

#include <array>
#include <cmath>
#include <string>
#include <unordered_set>

#include "distributions/angular_gauss.hpp"

ClassicGrid::ClassicGrid(int grid_div, AngularGauss* distr) : div(grid_div), distr(distr) { //double (*distribution_calc)(std::array<double, 3>&)) {
    //div              = grid_div;
    //distr = distr;
    //function_on_grid = distribution_calc;

    GenerateGridAndEvaluateFunc();
}

// TODO parallel
void ClassicGrid::GenerateGridAndEvaluateFunc() {
    double angle = M_PI / div;

    points.reserve(div * (div + 1));

    for (int i = 0; i < div; ++i) { // i goes around the sphere
        double theta = 2 * i * angle;
        double si    = std::sin(theta);
        double ci    = std::cos(theta);

        for (int j = 0; j <= div; ++j) { // j goes from N to S of the unit sphere
            double phi = j * angle;
            double sj  = std::sin(phi);
            double cj  = std::cos(phi);

            std::array<double, 3> point = {ci * sj, si * sj, cj}; // X, Y, Z
            points[{i, j}] = point;
        }
    }

    values.reserve(div * (div + 1));
    for (const auto& it : points) {
        values[it.second] =  distr->Calc(it.second); //function_on_grid(it->second);
    }

    /* trapezium is stored as 4 vertices from A to D anticlockwise:
     *   A -<- D          A (== D)
     *  /       \   or   / \
     * B --->--- C      B - C
     */
    // trapeziums = std::vector<std::vector<std::array<std::array<int, 2>, 4>>>(div);
    trapeziums.reserve(div);

    for (int i = 0; i < div; ++i) {
        std::vector<std::array<std::array<int, 2>, 4>> trapezium_column;
        trapezium_column.reserve(div);

        for (int j = 0; j < div; ++j) {
            int right_i = i < div - 1 ? i + 1 : 0;

            std::array<int, 2> A = {i, j};
            std::array<int, 2> B = {i, j + 1};
            std::array<int, 2> C = {right_i, j + 1};
            std::array<int, 2> D = {right_i, j};

            trapezium_column.push_back({A, B, C, D});
        }

        trapeziums.push_back(trapezium_column);
    }
}

// GetCoordsOfPoint converts spherical coordinates of a point to Cartesian.
// TODO @classmethod or static
std::array<double, 3> ClassicGrid::GetCoordsOfPoint(const std::array<double, 2>& point) {
    double theta = point[0], phi = point[1];

    double si = std::sin(theta);
    double ci = std::cos(theta);

    double sj = std::sin(phi);
    double cj = std::cos(phi);

    return {ci * sj, si * sj, cj}; // X, Y, Z
}

// GetAnglesOfPoint converts Cartesian coordinates of a point to spherical.
// TODO @classmethod or static
std::array<double, 2> ClassicGrid::GetAnglesOfPoint(const std::array<double, 3>& point) {
    double theta = std::atan2(point[1], point[0]); // arctg(y/x)
    if (theta < 0) {
        theta += 2 * M_PI; // due to our classic grid construction logic
    }
    double phi = std::acos(point[2]); // arccos(z)

    return {theta, phi};
}

// GetIsolineCoords returns points' coordinates for isoline with value 'iso_value'.
// It uses marching squares algorithm (https://en.wikipedia.org/wiki/Marching_squares).
std::vector<std::array<double, 3>> ClassicGrid::GetIsolineCoords(double iso_value) {
    const size_t initial_isoline_points_amount = 1000; // TODO move somewhere
    std::vector<std::array<double, 3>> isoline_points;
    isoline_points.reserve(initial_isoline_points_amount);

    int end_side; // TODO rename end
    std::array<double, 3> end_point;
    std::array<int, 2> initial_trapezium_indices;

    bool starting_trapezium_found                                     = false;
    const std::unordered_set<int> indicies_of_bad_to_start_trapeziums = {0, 15, 5, 10};

    for (int i = 0; i < trapeziums.size() && !starting_trapezium_found; ++i) {
        for (int j = 0; j < trapeziums[i].size(); ++j) {
            auto trapezium = trapeziums[i][j];
            auto index     = GetTrapeziumIndex(trapezium, iso_value);

            if (indicies_of_bad_to_start_trapeziums.find(index) == indicies_of_bad_to_start_trapeziums.end()) {
                std::tie(end_side, end_point) = ProcessSegment(trapezium, index, -1, iso_value);
                isoline_points.push_back(end_point);
                initial_trapezium_indices = {i, j};

                starting_trapezium_found = true;
                break;
            }
        }
    }

    int side;
    std::array<int, 2> trapezium_indices;
    std::tie(side, trapezium_indices) = GetNextTrapeziumIndices(end_side, initial_trapezium_indices);

    while (trapezium_indices != initial_trapezium_indices) {
        auto trapezium = trapeziums[trapezium_indices[0]][trapezium_indices[1]];
        auto index     = GetTrapeziumIndex(trapezium, iso_value);

        std::array<double, 3> point;
        std::tie(side, point) = ProcessSegment(trapezium, index, side, iso_value);
        isoline_points.push_back(point);
        std::tie(side, trapezium_indices) = GetNextTrapeziumIndices(side, trapezium_indices);
    }

    return isoline_points;
}

std::tuple<int, std::array<int, 2>> ClassicGrid::GetNextTrapeziumIndices(int prev_side,
                                                                         const std::array<int, 2>& prev_trapezium_indices) {
    switch (prev_side) {
    case 0: // we must get the one above
        return std::make_tuple(2, GetUpperTrapeziumIndices(prev_trapezium_indices));

    case 1: // get the right one
        return std::make_tuple(3, GetRightTrapeziumIndices(prev_trapezium_indices));

    case 2: // get the lower one
        return std::make_tuple(0, GetLowerTrapeziumIndices(prev_trapezium_indices));

    case 3: // get the left one
        return std::make_tuple(1, GetLeftTrapeziumIndices(prev_trapezium_indices));

    default:
        throw std::invalid_argument("GetNextTrapeziumIndices: got invalid prev_side argument: " +
                                    std::to_string(prev_side));
    }
}

std::array<int, 2> ClassicGrid::GetUpperTrapeziumIndices(const std::array<int, 2>& prev_trapezium_indices) {
    int i = prev_trapezium_indices[0], j = prev_trapezium_indices[1];
    assert(j != 0);
    return {i, j - 1};
}

std::array<int, 2> ClassicGrid::GetLowerTrapeziumIndices(const std::array<int, 2>& prev_trapezium_indices) {
    int i = prev_trapezium_indices[0], j = prev_trapezium_indices[1];
    assert(j != trapeziums[i].size() - 1);
    return {i, j + 1};
}

std::array<int, 2> ClassicGrid::GetLeftTrapeziumIndices(const std::array<int, 2>& prev_trapezium_indices) {
    int i = prev_trapezium_indices[0], j = prev_trapezium_indices[1];
    if (i == 0) {
        i = trapeziums.size();
    }
    return {i - 1, j};
}

std::array<int, 2> ClassicGrid::GetRightTrapeziumIndices(const std::array<int, 2>& prev_trapezium_indices) {
    int i = prev_trapezium_indices[0], j = prev_trapezium_indices[1];
    if (i == trapeziums.size() - 1) {
        i = -1;
    }
    return {i + 1, j};
}

// Linear interpolation to find intersection point with value c between vertices with values a and b
// where x1 = f^(-1)(a); x2 = f^(-1)(b); sign((a - c) * (b - c)) should be -1.
// Note that x1 and x2 - Cartesian coordinates of a point.
std::array<double, 3> ClassicGrid::FindIntersection(
        const std::array<double, 3>& x1, const std::array<double, 3>& x2, double a, double b, double c) {
    double k = (c - a) / (b - a);
    if (k >= 1.) {
        std::cerr << "FindIntersection, k=" << k << ", x1=" << x1.data() << ", x2=" << x2.data()
                  << ", a=" << a << ", b=" << b << ", c=" << c << std::endl;
        throw std::invalid_argument("FindIntersection: got invalid k from arguments: " + std::to_string(k));
    }

    // coordinates of intersection point
    double x = (1 - k) * x1[0] + k * x2[0];
    double y = (1 - k) * x1[1] + k * x2[1];
    double z = (1 - k) * x1[2] + k * x2[2];

    return {x, y, z};
}

// start_side is 0 - AD, 1 - DC, 2 - CB, 3 - BA
// end_side, end_point are returned
std::tuple<int, std::array<double, 3>> ClassicGrid::ProcessSegment(
        const std::array<std::array<int, 2>, 4>& trapezium, int index, int start_side, double iso_value) {
    int end_side = -1;
    std::array<double, 3> point_1, point_2;

    // 0 -- 1    1 -- 0
    // |    | or |    |
    // 1 -- 1    0 -- 0
    if (index == 7 || index == 8) {
        if (start_side == 0) {
            end_side = 3;
            // calc A and B
            point_1 = points[trapezium[0]];
            point_2 = points[trapezium[1]];
        } else {
            end_side = 0;
            // calc A and D
            point_1 = points[trapezium[0]];
            point_2 = points[trapezium[3]];
        }
    }

    // 0 -- 1    1 -- 0
    // |    | or |    |
    // 0 -- 0    1 -- 1
    else if (index == 4 || index == 11) {
        if (start_side == 0) {
            end_side = 1;
            // calc C and D
            point_1 = points[trapezium[2]];
            point_2 = points[trapezium[3]];
        } else {
            end_side = 0;
            // calc A and D
            point_1 = points[trapezium[0]];
            point_2 = points[trapezium[3]];
        }
    }

    // 0 -- 0    1 -- 1
    // |    | or |    |
    // 0 -- 1    1 -- 0
    else if (index == 2 || index == 13) {
        if (start_side == 1) {
            end_side = 2;
            // calc B and C
            point_1 = points[trapezium[1]];
            point_2 = points[trapezium[2]];
        } else {
            end_side = 1;
            // calc C and D
            point_1 = points[trapezium[2]];
            point_2 = points[trapezium[3]];
        }
    }

    // 0 -- 0    1 -- 1
    // |    | or |    |
    // 1 -- 0    0 -- 1
    else if (index == 1 || index == 14) {
        if (start_side == 2) {
            end_side = 3;
            // calc A and B
            point_1 = points[trapezium[0]];
            point_2 = points[trapezium[1]];
        } else {
            end_side = 2;
            // calc B and C
            point_1 = points[trapezium[1]];
            point_2 = points[trapezium[2]];
        }
    }

    // 0 -- 1    1 -- 0
    // |    | or |    |
    // 0 -- 1    1 -- 0
    else if (index == 6 || index == 9) {
        if (start_side == 0) {
            end_side = 2;
            // calc B and C
            point_1 = points[trapezium[1]];
            point_2 = points[trapezium[2]];
        } else {
            end_side = 0;
            // calc A and D
            point_1 = points[trapezium[0]];
            point_2 = points[trapezium[3]];
        }
    }

    // 0 -- 0    1 -- 1
    // |    | or |    |
    // 1 -- 1    0 -- 0
    else if (index == 3 || index == 12) {
        if (start_side == 1) {
            end_side = 3;
            // calc A and B
            point_1 = points[trapezium[0]];
            point_2 = points[trapezium[1]];
        } else {
            end_side = 1;
            // calc C and D
            point_1 = points[trapezium[2]];
            point_2 = points[trapezium[3]];
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
            point_1 = points[trapezium[0]];
            point_2 = points[trapezium[1]];
        } else if (start_side == 3) {
            end_side = 0;
            // calc A and D
            point_1 = points[trapezium[0]];
            point_2 = points[trapezium[3]];
        } else if (start_side == 1) {
            end_side = 2;
            // calc B and C
            point_1 = points[trapezium[1]];
            point_2 = points[trapezium[2]];
        } else { // start_side == 2
            end_side = 1;
            // calc C and D
            point_1 = points[trapezium[2]];
            point_2 = points[trapezium[3]];
        }
    }

    // 0 -- 1
    // |    |
    // 1 -- 0
    else if (index == 10) {
        if (start_side == 0) {
            end_side = 1;
            // calc C and D
            point_1 = points[trapezium[2]];
            point_2 = points[trapezium[3]];
        } else if (start_side == 1) {
            end_side = 0;
            // calc A and D
            point_1 = points[trapezium[0]];
            point_2 = points[trapezium[3]];
        } else if (start_side == 2) {
            end_side = 3;
            // calc A and B
            point_1 = points[trapezium[0]];
            point_2 = points[trapezium[1]];
        } else { // start_side == 3
            end_side = 2;
            // calc B and C
            point_1 = points[trapezium[1]];
            point_2 = points[trapezium[2]];
        }
    }

    else {
        std::cerr << "ProcessSegment, trapezium=" << trapezium.data() << ", index=" << index
                  << ", start_side=" << start_side << ", iso_value=" << iso_value << std::endl;
        throw std::invalid_argument("ProcessSegment: got invalid index/start_side: " + std::to_string(index) +
                                    " " + std::to_string(start_side));
    }

    auto end_point = FindIntersection(point_1, point_2, values[point_1], values[point_2], iso_value);
    return std::make_tuple(end_side, end_point);
}

// GetTrapeziumIndex return the trapezium index needed for marching squares algorithm.
int ClassicGrid::GetTrapeziumIndex(const std::array<std::array<int, 2>, 4>& trapezium, double value) {
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
