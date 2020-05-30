#include "classic_grid.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <string>
#include <thread>
#include <unordered_set>
#include <vector>

#include <iostream> // TODO
extern "C" {
#include "s2kit/cospml.h" // TODO
#include "s2kit/FST_semi_memo.h" // TODO
}

#include "distributions/angular_gauss.hpp"    // AngularGauss
#include "distributions/von_mises_fisher.hpp" // VonMisesFisher
#include "types/types.hpp"
#include "util/cfg.hpp"  // cfg::kConfig
#include "util/util.hpp" // Bounds()

ClassicGrid::ClassicGrid(int grid_div, AngularGauss* distr) : div(grid_div), distr(distr), fisher_distr(nullptr) {
    assert(div % 2 == 0);
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

    points.reserve(div * div + 2);
    cartesian_points.reserve(div * div + 2);

    // north pole
    cartesian_points.push_back({0, 0, 1});
    points[{-1, -1}] = cartesian_points.size()-1;

    for (int i = 0; i < div; ++i) { // i goes from N to S
        double theta = (i + 0.5) * angle;
        double si    = std::sin(theta);
        double ci    = std::cos(theta);

        for (int j = 0; j < div; ++j) { // j goes around the sphere
            double phi = 2 * j * angle;
            double sj  = std::sin(phi);
            double cj  = std::cos(phi);

            CoordsOfPoint point = {cj * si, sj * si, ci}; // X, Y, Z
            cartesian_points.push_back(point);
            points[{i, j}] = cartesian_points.size()-1;
        }
    }

    // south pole
    cartesian_points.push_back({0, 0, -1});
    points[{div, div}] = cartesian_points.size()-1;

    values.resize(div * div + 2);

    int threads_amount = cfg::kConfig["threads"].get<int>();
    auto bounds = Bounds(threads_amount, cartesian_points.size());
    auto threads = std::vector<std::thread>(threads_amount-1);

    for (int i = 0; i < threads_amount-1; ++i)
        threads[i] = std::thread(&ClassicGrid::EvaluateFuncRoutine, this, bounds[i], bounds[i+1]);

    EvaluateFuncRoutine(bounds[threads_amount-1], bounds[threads_amount]);

    for (int i = 0; i < threads_amount-1; ++i)
        threads[i].join();

    /*
     * trapezium is stored as 4 vertices from A to D anticlockwise:
     *    A -<- D
     *   /       \
     *  B --->--- C
     * 
     * each trapezium in ring on north pole is:
     *    A (== D)
     *   / \
     *  B - C
     * where A (D) has indices {-1, -1} (north pole)
     * 
     * each trapezium in ring on south pole is:
     *  A - D
     *   \ /
     *    B (== C)
     * where B (C) has indices {div, div} (south pole)
     */

    trapeziums.reserve(div + 1);

    std::vector<Trapezium> north_trapezium_ring;
    north_trapezium_ring.reserve(div);
    for (int j = 0; j < div; ++j) { // ring on north pole
        int right_j = j < div - 1 ? j + 1 : 0;

        ClassicGridPoint A = {-1, -1};
        ClassicGridPoint B = {0, j};
        ClassicGridPoint C = {0, right_j};
        ClassicGridPoint D = {-1, -1};

        north_trapezium_ring.push_back({A, B, C, D});
    }
    trapeziums.push_back(north_trapezium_ring);

    for (int i = 0; i < div-1; ++i) {
        std::vector<Trapezium> trapezium_ring;
        trapezium_ring.reserve(div);

        for (int j = 0; j < div; ++j) {
            int right_j = j < div - 1 ? j + 1 : 0;

            ClassicGridPoint A = {i, j};
            ClassicGridPoint B = {i + 1, j};
            ClassicGridPoint C = {i + 1, right_j};
            ClassicGridPoint D = {i, right_j};

            trapezium_ring.push_back({A, B, C, D});
        }

        trapeziums.push_back(trapezium_ring);
    }

    std::vector<Trapezium> south_trapezium_ring;
    south_trapezium_ring.reserve(div);
    for (int j = 0; j < div; ++j) { // ring on south pole
        int right_j = j < div - 1 ? j + 1 : 0;

        ClassicGridPoint A = {div-1, j};
        ClassicGridPoint B = {div, div};
        ClassicGridPoint C = {div, div};
        ClassicGridPoint D = {div-1, right_j};

        south_trapezium_ring.push_back({A, B, C, D});
    }

    trapeziums.push_back(south_trapezium_ring);
}

void ClassicGrid::EvaluateFisherFuncRoutine(int lower_bound, int upper_bound) {
    for (int i = lower_bound; i < upper_bound; ++i)
        fisher_values[i] = fisher_distr->Calc(cartesian_points[i]);
}

void ClassicGrid::SetFisherKappa(double kappa) {
    fisher_distr = new VonMisesFisher(kappa);

    fisher_values.resize(div * div + 2);

    int threads_amount = cfg::kConfig["threads"].get<int>();
    auto bounds = Bounds(threads_amount, cartesian_points.size());
    auto threads = std::vector<std::thread>(threads_amount-1);

    for (int i = 0; i < threads_amount-1; ++i)
        threads[i] = std::thread(&ClassicGrid::EvaluateFisherFuncRoutine, this, bounds[i], bounds[i+1]);

    EvaluateFisherFuncRoutine(bounds[threads_amount-1], bounds[threads_amount]);

    for (int i = 0; i < threads_amount-1; ++i)
        threads[i].join();

    int values_size = div * div;

    double* values_start = &values[1]; // skip north pole extra value
    std::vector<double> i_values(values_size, 0.);
    double* i_values_start = &i_values[0];

    double* fisher_values_start = &fisher_values[1]; // skip north pole extra value
    std::vector<double> i_fisher_values(values_size, 0.);
    double* i_fisher_values_start = &i_fisher_values[0];

    values_after_convolution.resize(values_size + 2);
    double* values_after_convolution_start = &values_after_convolution[1]; // skip north pole extra value
    std::vector<double> i_values_after_convolution(values_size, 0.);
    double* i_values_after_convolution_start = &i_values_after_convolution[0];

    int bw = div / 2;
    int cutoff = bw; // seminaive all orders
    int legendre_size = Reduced_Naive_TableSize(bw, cutoff) + Reduced_SpharmonicTableSize(bw, cutoff);
    int workspace_size = 2 * legendre_size + 12 * bw * bw + 12 * bw;
    std::vector<double> workspace(workspace_size, 0.);
    double* workspace_start = &workspace[0];

    ConvOn2SphereSemiMemo(values_start, i_values_start, fisher_values_start, i_fisher_values_start,
        values_after_convolution_start, i_values_after_convolution_start, bw, workspace_start);

    // add value on north and south poles as simple mean of neighbour values

    double north_value = 0;
    for (int i = 0; i < trapeziums[0].size(); ++i) {
        ClassicGridPoint B = trapeziums[0][i][1];
        ClassicGridPoint C = trapeziums[0][i][2];
        north_value += values_after_convolution[points[B]] + values_after_convolution[points[C]];
    }
    north_value /= 2 * trapeziums[0].size();
    values_after_convolution[0] = north_value;

    double south_value = 0;
    for (int i = 0; i < trapeziums[div].size(); ++i) {
        ClassicGridPoint A = trapeziums[0][i][0];
        ClassicGridPoint D = trapeziums[0][i][3];
        south_value += values_after_convolution[points[A]] + values_after_convolution[points[D]];
    }
    south_value /= 2 * trapeziums[div].size();
    values_after_convolution[values_after_convolution.size()-1] = south_value;
}

// assumes that point is normed due to GetTrapeziumIndex() // TODO really?
double ClassicGrid::CalcFunc(const CoordsOfPoint& u) {
    if (fisher_distr == nullptr || values_after_convolution.size() != values.size()) {
        // std::cout << "WITHOUT CONV\n"; // TODO debug
        return distr->Calc(u);
    }

    // we must find the trapezium to which u belongs and interpolate function
    // value using values_after_convolution

    auto trapezium_index = GetTrapeziumIndex(u);
    auto trapezium = trapeziums[trapezium_index[0]][trapezium_index[1]];

    // maybe u is one of the grids points (one of trapezium's points)
    for (int i = 0; i < trapezium.size(); ++i)
        if (u == cartesian_points[points[trapezium[i]]])
            return values_after_convolution[points[trapezium[i]]];

    if (trapezium_index[0] == 0) { // north ring
        return InterpolateInsideABC(u, trapezium[0], trapezium[1], trapezium[2]);
    } else if (trapezium_index[0] == trapeziums.size()-1) { // south ring
        return InterpolateInsideABC(u, trapezium[0], trapezium[3], trapezium[2]);
    }

    auto A_coords = cartesian_points[points[trapezium[0]]];
    auto C_coords = cartesian_points[points[trapezium[2]]];

    // we should detect to which triangle inside trapezium u belongs:
    // ABC or ADC
    Vector AU = GetVector(A_coords, u);
    Vector AC = GetVector(A_coords, C_coords);
    Vector OA = GetVector({0, 0, 0}, A_coords);
    Vector Across; // cross product of AU and AC vectors
    Across[0] = AU[1] * AC[2] - AU[2] * AC[1];
    Across[1] = AU[2] * AC[0] - AU[0] * AC[2];
    Across[2] = AU[0] * AC[1] - AU[1] * AC[0];

    double Across_projection_on_OA = 0; // dot product of Across and OA
    for (int i = 0; i < 3; ++i)
        Across_projection_on_OA += Across[i] * OA[i];
    
    if (Across_projection_on_OA > 0) // then u is inside ABC
        return InterpolateInsideABC(u, trapezium[0], trapezium[1], trapezium[2]);
    
    // u is inside ADC
    return InterpolateInsideABC(u, trapezium[0], trapezium[3], trapezium[2]);
}

TrapeziumIndex ClassicGrid::GetTrapeziumIndex(const CoordsOfPoint& u) const {
    auto angles_of_point = GetAnglesOfPoint(u);
    double phi = angles_of_point[0], theta = angles_of_point[1];
    double angle = M_PI / div;

    int i = static_cast<int>(std::floor(theta / angle - 0.5)) + 1;
    int j = static_cast<int>(std::floor(phi / (2 * angle))) % div;

    return {i, j};
}

CoordsOfPoint ClassicGrid::GetCoordsOfPoint(const AnglesOfPoint& point) {
    double phi = point[0], theta = point[1];

    double si = std::sin(phi);
    double ci = std::cos(phi);

    double sj = std::sin(theta);
    double cj = std::cos(theta);

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
    const size_t initial_isoline_points_amount = 2 * div;
    std::vector<CoordsOfPoint> isoline_points;
    isoline_points.reserve(initial_isoline_points_amount);

    // TODO code duplication can be reduced

    auto initial_trapezium_index = FindTrapeziumWithIsolineInside(iso_value);
    auto initial_trapezium = trapeziums[initial_trapezium_index[0]][initial_trapezium_index[1]];
    auto initial_marching_index = GetMarchingSquareIndex(initial_trapezium, iso_value);
    auto [initial_trapezium_side, initial_point] = ProcessSegment(initial_trapezium, initial_marching_index, -1, iso_value);
    isoline_points.push_back(initial_point);

    auto [side, trapezium_index] = GetNextTrapeziumIndex(initial_trapezium_side, initial_trapezium_index);

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

TrapeziumIndex ClassicGrid::FindTrapeziumWithIsolineInside(double iso_value) {
    int initial_trapezium_side;
    TrapeziumIndex initial_trapezium_index;

    const std::unordered_set<int> marching_indicies_of_bad_to_start_trapeziums = {0, 15, 5, 10};

    // TODO probably better solution for unimodal functions is to use point
    // with actual pdf maximum instead of mean vector

    Vector mean = distr->Mean();
    double squared_sum = 0.;
    for (int i = 0; i < 3; ++i)
        squared_sum += mean[i] * mean[i];

    double mean_len = std::sqrt(squared_sum);
    Vector normed_mean;
    for (int i = 0; i < 3; ++i)
        normed_mean[i] = mean[i] / mean_len;

    auto trapezium_with_mean_index = GetTrapeziumIndex(normed_mean);
    int i_with_mean = trapezium_with_mean_index[0]; // polar
    int j_with_mean = trapezium_with_mean_index[1]; // azimuthal

    for (int i = i_with_mean; i < trapeziums.size(); ++i) {
        auto trapezium = trapeziums[i][j_with_mean];
        auto marching_index = GetMarchingSquareIndex(trapezium, iso_value);

        if (marching_indicies_of_bad_to_start_trapeziums.find(marching_index) == marching_indicies_of_bad_to_start_trapeziums.end())
            return {i, j_with_mean};
    }

    // I have explicitly split one loop into two, since starting search from
    // i_with_mean to trapeziums.size() may return needed trapezium faster
    for (int i = 0; i < i_with_mean; ++i) {
        auto trapezium = trapeziums[i][j_with_mean];
        auto marching_index = GetMarchingSquareIndex(trapezium, iso_value);

        if (marching_indicies_of_bad_to_start_trapeziums.find(marching_index) == marching_indicies_of_bad_to_start_trapeziums.end())
            return {i, j_with_mean};
    }

    // We are unlucky and have to check the opposite azimuth
    int j_opposite_to_mean = (j_with_mean + div/2) % div;
    for (int i = 0; i < trapeziums.size(); ++i) {
        auto trapezium = trapeziums[i][j_opposite_to_mean];
        auto marching_index = GetMarchingSquareIndex(trapezium, iso_value);

        if (marching_indicies_of_bad_to_start_trapeziums.find(marching_index) == marching_indicies_of_bad_to_start_trapeziums.end())
            return {i, j_with_mean};
    }

    // TODO maybe should leave full scan of all trapeziums?

    throw std::invalid_argument("FindTrapeziumWithIsolineInside: could not find trapezium with isovalue=" + std::to_string(iso_value));

    /*
    for (int i = 0; i < trapeziums.size(); ++i) {
        for (int j = 0; j < trapeziums[i].size(); ++j) {
            auto trapezium = trapeziums[i][j];
            auto marching_index = GetMarchingSquareIndex(trapezium, iso_value);

            if (marching_indicies_of_bad_to_start_trapeziums.find(marching_index) == marching_indicies_of_bad_to_start_trapeziums.end())
                return {i, j};
        }
    }
    */
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
    assert(i != 0);
    return {i - 1, j};
}

TrapeziumIndex ClassicGrid::GetLowerTrapeziumIndex(const TrapeziumIndex& previous) const {
    int i = previous[0], j = previous[1];
    assert(i != trapeziums.size()-1);
    return {i + 1, j};
}

TrapeziumIndex ClassicGrid::GetLeftTrapeziumIndex(const TrapeziumIndex& previous) const {
    int i = previous[0], j = previous[1];
    if (j == 0)
        return {i, static_cast<int>(trapeziums[i].size()-1)};

    return {i, j - 1};
}

TrapeziumIndex ClassicGrid::GetRightTrapeziumIndex(const TrapeziumIndex& previous) const {
    int i = previous[0], j = previous[1];
    if (j == trapeziums[i].size()-1)
        return {i, 0};

    return {i, j + 1};
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


// TODO should move to util
Vector ClassicGrid::GetVector(const CoordsOfPoint& x1, const CoordsOfPoint& x2) const {
    return {x2[0] - x1[0], x2[1] - x1[1], x2[2] - x1[2]};
}

double ClassicGrid::InterpolateInsideABC(const CoordsOfPoint& u, const ClassicGridPoint& A, const ClassicGridPoint& B, const ClassicGridPoint& C) {
    std::cout << "enter InterpolateInsideABC\n";

    auto A_coords = cartesian_points[points[A]];
    auto B_coords = cartesian_points[points[B]];
    auto C_coords = cartesian_points[points[C]];

    Vector AB = GetVector(A_coords, B_coords);
    Vector AC = GetVector(A_coords, C_coords);
    Vector AU = GetVector(A_coords, u);

    double AB_len = std::sqrt(AB[0]*AB[0] + AB[1]*AB[1] + AB[2]*AB[2]);

    double AC_square = AC[0]*AC[0] + AC[1]*AC[1] + AC[2]*AC[2];
    double AU_square = AU[0]*AU[0] + AU[1]*AU[1] + AU[2]*AU[2];

    // recoordinate A, B, C and u, since they are in the same plane

    std::array<double, 2> p1 = {0, 0}; // A
    std::array<double, 2> p2 = {AB_len, 0}; // B

    double p3x = (AC[0]*AB[0] + AC[1]*AB[1] + AC[2]*AB[2]) / AB_len; // Cx
    double p3y = std::sqrt(AC_square - p3x*p3x); // Cy

    double ux = (AU[0]*AB[0] + AU[1]*AB[1] + AU[2]*AB[2]) / AB_len;
    double uy = std::sqrt(AU_square - ux*ux);

    double det = (p2[0] - p1[0])*(p3y - p1[1]) - (p3x - p1[0])*(p2[1] - p1[1]);
    double det1 = (p2[0]*p3y - p3x*p2[1]) + ux*(p2[1] - p3y) + uy*(p3x - p2[0]);
    double det2 = (p3x*p1[1] - p1[0]*p3y) + ux*(p3y - p1[1]) + uy*(p1[0] - p3x);

    // barycentric coordinates of the point u
    double alpha = det1 / det;
    double beta = det2 / det;
    double gamma = 1 - alpha - beta;

    // TODO DECREASES PERFORMANCE!
    auto values = this->values_after_convolution;
    if (fisher_distr == nullptr || this->values_after_convolution.size() != this->values.size()) {
        // std::cout << "InterpolateInsideABC: WITHOUT CONV\n";
        values = this->values;
    }

    return alpha*values[points[A]] + beta*values[points[B]] + gamma*values[points[C]];
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

    // TODO DECREASES PERFORMANCE!
    //auto values = this->values_after_convolution;
    //if (fisher_distr == nullptr || this->values_after_convolution.size() != this->values.size()) {
    //    // std::cout << "ProcessSegment: WITHOUT CONV\n";
    //    values = this->values;
    //}

    auto end_point = FindIntersection(point_1, point_2, values[point_1_index], values[point_2_index], iso_value);
    return std::make_tuple(end_side, end_point);
}

int ClassicGrid::GetMarchingSquareIndex(const Trapezium& trapezium, double value) {
    int index = 0;

    // TODO DECREASES PERFORMANCE!
    //auto values = this->values_after_convolution;
    //if (fisher_distr == nullptr || this->values_after_convolution.size() != this->values.size()) {
    //    // std::cout << "GetMarchingSquareIndex: WITHOUT CONV\n";
    //    values = this->values;
    //}

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
