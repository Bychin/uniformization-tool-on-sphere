#include "stats.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <limits>
#include <vector>

#include "types/types.hpp"

StatsAPI::StatsAPI(std::vector<CoordsOfPoint>& points, ClassicGrid* classic_grid, SpiralGrid* spiral_grid)
    : points(points), classic_grid(classic_grid), spiral_grid(spiral_grid) {
}

double StatsAPI::DistanceBetweenPoints(const CoordsOfPoint& point1, const CoordsOfPoint& point2) const {
    double sum = 0.;
    for (int i = 0; i < 3; ++i) {
        double diff = point1[i] - point2[i];
        sum += diff * diff;
    }
    
    return std::sqrt(sum);
}

int StatsAPI::InsertPointIntoIsoline(const CoordsOfPoint& point, IsolineCoords& isoline_points) const {
    int point_index = -1;

    if (isoline_points.size() < 2)
        throw std::logic_error("InsertPointIntoIsoline: isoline_points length is less than 2");

    // Calculate distance between point and each isoline point.
    // The smallest two distances will indicate the place, where the point must
    // be inserted into isoline_points array.

    double first_smallest_distance = DistanceBetweenPoints(isoline_points[0], point);
    double second_smallest_distance = DistanceBetweenPoints(isoline_points[1], point);
    int index_of_first_smallest = 0;
    int index_of_second_smallest = 1;

    if (first_smallest_distance > second_smallest_distance) {
        std::swap(first_smallest_distance, second_smallest_distance);
        std::swap(index_of_first_smallest, index_of_second_smallest);
    }

    for (int i = 2; i < isoline_points.size(); ++i) {
        double distance = DistanceBetweenPoints(isoline_points[i], point);
        if (distance < first_smallest_distance) {
            second_smallest_distance = first_smallest_distance;
            first_smallest_distance = distance;

            index_of_second_smallest = index_of_first_smallest;
            index_of_first_smallest = i;

        } else if (distance < second_smallest_distance) {
            second_smallest_distance = distance;
            index_of_second_smallest = i;
        }
    }

    int indices_diff = std::abs(index_of_first_smallest - index_of_second_smallest);
    if (indices_diff == isoline_points.size() - 1)
        // it is the first and the last points of isoline_points
        point_index = isoline_points.size();
    else if (indices_diff == 1)
        // it is the two adjacent points of isoline_points
        point_index = std::max(index_of_first_smallest, index_of_second_smallest);
    
    if (point_index == -1)
        throw std::logic_error("InsertPointIntoIsoline: the smallest two distances are not nearby");

    isoline_points.insert(isoline_points.begin() + point_index, point);

    return point_index;
}

Vector StatsAPI::GetNormalizedVector(const CoordsOfPoint& point1, const CoordsOfPoint& point2) const {
    Vector vec;
    double squared_sum = 0.;

    for (int i = 0; i < 3; ++i) {
        double value = point2[i] - point1[i];

        vec[i] = value;
        squared_sum += value * value;
    }

    double vec_len = std::sqrt(squared_sum);
    for (int i = 0; i < 3; ++i)
        vec[i] /= vec_len;

    return vec;
}

int StatsAPI::GetIndexOfPrevIsolinePoint(int point_index, const IsolineCoords& isoline_points) const {
    if (point_index == 0)
        return isoline_points.size() - 1;
    return point_index - 1;
}

int StatsAPI::GetIndexOfNextIsolinePoint(int point_index, const IsolineCoords& isoline_points) const {
    if (point_index == isoline_points.size() - 1)
        return 0;
    return point_index + 1;
}

bool StatsAPI::CheckClockwiseDirection(const Vector& M, const Vector& I, const Vector& A) const {
    Vector AI; // a cross product of A and I vectors
    AI[0] = A[1] * I[2] - A[2] * I[1]; 
    AI[1] = A[2] * I[0] - A[0] * I[2]; 
    AI[2] = A[0] * I[1] - A[1] * I[0];

    double AI_projection_on_M = 0; // a dot product of AI and M
    for (int i = 0; i < 3; ++i)
        AI_projection_on_M += AI[i] * M[i];
    
    if (AI_projection_on_M > 0)
        return true;
    return false;
}

double StatsAPI::CalcIntegralOnInfinitesimalCurve(const CoordsOfPoint& point1, const CoordsOfPoint& point2) const {
    CoordsOfPoint center_point;
    for (int i = 0; i < 3; ++i)
        center_point[i] = (point1[i] + point2[i]) / 2;

    double gradient_modulus = 1; // TODO np.linalg.norm(self.grad_of_d(middle_point)) // self.grad_of_d = grad(self.distribution.calcForGradient)
    // Градиент считается приближенно, через величины функции на сетке --- взять соседние точки и решить линейное уравнение, например.

    return DistanceBetweenPoints(point1, point2) / gradient_modulus;
}

double StatsAPI::CalcIntegralOnFullIsoline(const IsolineCoords& isoline_points) const {
    double sum = 0;
    for (int i = 0; i < isoline_points.size()-1; ++i)
        sum += CalcIntegralOnInfinitesimalCurve(isoline_points[i], isoline_points[i+1]);

    return sum + CalcIntegralOnInfinitesimalCurve(isoline_points[isoline_points.size()-1], isoline_points[0]);
}

double StatsAPI::CalcIntegralOnLeftCurveOnIsoline(const IsolineCoords& isoline_points, int start_index, int end_index) const {
    if (start_index == end_index)
        throw std::invalid_argument("CalcIntegralOnLeftCurveOnIsoline: start_index equals end_index, start_index="+std::to_string(start_index));

    double sum = 0;

    if (end_index < start_index) {
        for (int i = end_index; i < start_index; ++i)
            sum += CalcIntegralOnInfinitesimalCurve(isoline_points[i], isoline_points[i+1]);
        return sum;
    }

    for (int i = 0; i < start_index; ++i)
        sum += CalcIntegralOnInfinitesimalCurve(isoline_points[i], isoline_points[i+1]);
    for (int i = end_index; i < isoline_points.size()-1; ++i)
        sum += CalcIntegralOnInfinitesimalCurve(isoline_points[i], isoline_points[i+1]);
    return sum + CalcIntegralOnInfinitesimalCurve(isoline_points[isoline_points.size()-1], isoline_points[0]);
}

double StatsAPI::CalcIntegralOnRightCurveOnIsoline(const IsolineCoords& isoline_points, int start_index, int end_index) const {
    if (start_index == end_index)
        throw std::invalid_argument("CalcIntegralOnRightCurveOnIsoline: start_index equals end_index, start_index="+std::to_string(start_index));

    double sum = 0;

    if (start_index < end_index) {
        for (int i = start_index; i < end_index; ++i)
            sum += CalcIntegralOnInfinitesimalCurve(isoline_points[i], isoline_points[i+1]);
        return sum;
    }

    for (int i = 0; i < end_index; ++i)
        sum += CalcIntegralOnInfinitesimalCurve(isoline_points[i], isoline_points[i+1]);
    for (int i = start_index; i < isoline_points.size()-1; ++i)
        sum += CalcIntegralOnInfinitesimalCurve(isoline_points[i], isoline_points[i+1]);
    return sum + CalcIntegralOnInfinitesimalCurve(isoline_points[isoline_points.size()-1], isoline_points[0]);
}

double StatsAPI::CalculateSStat(const CoordsOfPoint& point, double value) const {
    auto isoline_points = classic_grid->GetIsolineCoords(value);
    int point_index = InsertPointIntoIsoline(point, isoline_points);

    CoordsOfPoint zero_point = {0, 0, 0};
    auto mean = classic_grid->Func()->Mean();
    auto normed_mean = GetNormalizedVector(zero_point, mean);
    auto angles_of_mean = classic_grid->GetAnglesOfPoint(normed_mean); // TODO should I change distribution mean with the actual PDF's maximum?

    double mean_phi = angles_of_mean[0]; // azimuthal angle
    double mean_theta = angles_of_mean[1]; // polar angle

    if (normed_mean[0] == 0 and normed_mean[1] == 0) // if mean vector points to one of the poles, correct it's polar angle
        mean_theta = 0;

    // see my graduation work in order to understand the logic below

    std::vector<std::array<double, 2>> angles_of_isoline_points;
    angles_of_isoline_points.reserve(isoline_points.size());
    for (const auto& p : isoline_points)
        angles_of_isoline_points.push_back(classic_grid->GetAnglesOfPoint(p));
    
    std::vector<int> indices_between_zero_and_theta_mean;
    indices_between_zero_and_theta_mean.reserve(isoline_points.size()); // not optimal but at least something
    for (int i = 0; i < angles_of_isoline_points.size(); ++i)
        if (0 <= angles_of_isoline_points[i][1] && angles_of_isoline_points[i][1] <= mean_theta)
            indices_between_zero_and_theta_mean.push_back(i);


    std::vector<int> indices_between_zero_and_pi_minus_theta_mean;
    indices_between_zero_and_pi_minus_theta_mean.reserve(isoline_points.size()); // not optimal but at least something
    for (int i = 0; i < angles_of_isoline_points.size(); ++i)
        if (0 <= angles_of_isoline_points[i][1] && angles_of_isoline_points[i][1] <= (M_PI - mean_theta))
            indices_between_zero_and_pi_minus_theta_mean.push_back(i);

    std::vector<std::array<double, 2>> intersection_points;
    intersection_points.reserve(2);

    double smallest_diff_to_theta_mean = -1; // maybe should use intersection_points.size() for check instead

    // WARNING: all differencies below are calculated as |a-b|, since a and b must be non-negative angles!

    if (indices_between_zero_and_theta_mean.size() >= 2) { // TODO what if only 1 point?
        // find the smallest two differences to theta mean

        double first_smallest_diff = std::fabs(angles_of_isoline_points[indices_between_zero_and_theta_mean[0]][0] - angles_of_mean[0]);
        double second_smallest_diff = std::fabs(angles_of_isoline_points[indices_between_zero_and_theta_mean[1]][0] - angles_of_mean[0]);

        int index_of_first_smallest = 0; // in indices_between_zero_and_theta_mean vector
        int index_of_second_smallest = 1;

        if (first_smallest_diff > second_smallest_diff) {
            std::swap(first_smallest_diff, second_smallest_diff);
            std::swap(index_of_first_smallest, index_of_second_smallest);
        }

        for (int i = 2; i < indices_between_zero_and_theta_mean.size(); ++i) {
            double diff = std::fabs(angles_of_isoline_points[indices_between_zero_and_theta_mean[i]][0] - angles_of_mean[0]);
            if (diff < first_smallest_diff) {
                second_smallest_diff = first_smallest_diff;
                first_smallest_diff = diff;

                index_of_second_smallest = index_of_first_smallest;
                index_of_first_smallest = i;

            } else if (diff < second_smallest_diff) {
                second_smallest_diff = diff;
                index_of_second_smallest = i;
            }
        }

        smallest_diff_to_theta_mean = first_smallest_diff;

        intersection_points.push_back(angles_of_isoline_points[indices_between_zero_and_theta_mean[index_of_first_smallest]]);
        intersection_points.push_back(angles_of_isoline_points[indices_between_zero_and_theta_mean[index_of_second_smallest]]);
    }

    if (indices_between_zero_and_pi_minus_theta_mean.size() >= 2) { // TODO what if only 1 point?
        double theta_mean_plus_pi = std::fmod(angles_of_mean[0] + M_PI, 2 * M_PI);
        // find the smallest two differences to (theta mean + PI)

        double first_smallest_diff = std::fabs(angles_of_isoline_points[indices_between_zero_and_pi_minus_theta_mean[0]][0] - theta_mean_plus_pi);
        double second_smallest_diff = std::fabs(angles_of_isoline_points[indices_between_zero_and_pi_minus_theta_mean[1]][0] - theta_mean_plus_pi);

        int index_of_first_smallest = 0; // in indices_between_zero_and_pi_minus_theta_mean vector
        int index_of_second_smallest = 1;

        if (first_smallest_diff > second_smallest_diff) {
            std::swap(first_smallest_diff, second_smallest_diff);
            std::swap(index_of_first_smallest, index_of_second_smallest);
        }

        for (int i = 2; i < indices_between_zero_and_pi_minus_theta_mean.size(); ++i) {
            double diff = std::fabs(angles_of_isoline_points[indices_between_zero_and_pi_minus_theta_mean[i]][0] - theta_mean_plus_pi);
            if (diff < first_smallest_diff) {
                second_smallest_diff = first_smallest_diff;
                first_smallest_diff = diff;

                index_of_second_smallest = index_of_first_smallest;
                index_of_first_smallest = i;

            } else if (diff <= second_smallest_diff) {
                second_smallest_diff = diff;
                index_of_second_smallest = i;
            }
        }

        if (smallest_diff_to_theta_mean != -1) {
            if (smallest_diff_to_theta_mean > first_smallest_diff) {
                intersection_points.clear();
                intersection_points.push_back(angles_of_isoline_points[indices_between_zero_and_pi_minus_theta_mean[index_of_first_smallest]]);
                intersection_points.push_back(angles_of_isoline_points[indices_between_zero_and_pi_minus_theta_mean[index_of_second_smallest]]);
            }
        } else {
            intersection_points.push_back(angles_of_isoline_points[indices_between_zero_and_pi_minus_theta_mean[index_of_first_smallest]]);
            intersection_points.push_back(angles_of_isoline_points[indices_between_zero_and_pi_minus_theta_mean[index_of_second_smallest]]);
        }
    }

    std::array<double, 2> intersection_point = {(intersection_points[0][0] + intersection_points[1][0]) / 2, (intersection_points[0][1] + intersection_points[1][1]) / 2};
    auto intersection_point_coords = classic_grid->GetCoordsOfPoint(intersection_point);
    int intersection_point_index = InsertPointIntoIsoline(intersection_point_coords, isoline_points);

    auto point_A = isoline_points[GetIndexOfNextIsolinePoint(intersection_point_index, isoline_points)];
    auto vec_A = GetNormalizedVector(normed_mean, point_A);
    auto vec_I = GetNormalizedVector(normed_mean, intersection_point_coords);

    double isoline_integral = CalcIntegralOnFullIsoline(isoline_points);
    double curve_integral = 0;
    bool clockwiseDirection = false;
    if (CheckClockwiseDirection(normed_mean, vec_I, vec_A)) {
        clockwiseDirection = true;
        curve_integral = CalcIntegralOnRightCurveOnIsoline(isoline_points, intersection_point_index, point_index);
    } else {
        curve_integral = CalcIntegralOnLeftCurveOnIsoline(isoline_points, intersection_point_index, point_index);
    }

    // return clockwiseDirection too
    // return {"isoline": isoline_points, "point": intersection_point_coords, "S": curve_integral / isoline_integral}
    return curve_integral / isoline_integral;
}

double StatsAPI::CalculateTStat(double value) const {
    return spiral_grid->CalcIntegralInsideIsoline(value);
}

const std::vector<CoordsOfPoint> StatsAPI::Points(void) const {
    return points;
}

void StatsAPI::Validate(void) const {
    if (classic_grid == nullptr) {
        throw std::invalid_argument("invalid StatsAPI: classic_grid is nullptr");
    }
    if (spiral_grid == nullptr) {
        throw std::invalid_argument("invalid StatsAPI: spiral_grid is nullptr");
    }
    if (points.size() == 0) {
        throw std::invalid_argument("invalid StatsAPI: empty points");
    }
}

double StatsAPI::CalculateFunc(const CoordsOfPoint& point) const {
    return spiral_grid->Func()->Calc(point);
}
