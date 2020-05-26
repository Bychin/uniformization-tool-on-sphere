#ifndef CLASSIC_GRID_HPP
#define CLASSIC_GRID_HPP

#include <array>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <boost/container_hash/hash.hpp>

#include "distributions/angular_gauss.hpp"

typedef std::array<int, 2> ClassicGridPoint;
typedef std::array<ClassicGridPoint, 4> Trapezium;

typedef std::array<double, 2> AnglesOfPoint;
typedef std::array<double, 3> CoordsOfPoint;

class ClassicGrid {
    int div;

    // points is a map of (i, j) pairs, that are indices of grid points, to
    // the same point in Cartesian coordinates' indicies
    std::unordered_map<ClassicGridPoint, int, boost::hash<ClassicGridPoint>> points;

    // cartesian_points is a vector of points in Cartesian coordinates
    std::vector<CoordsOfPoint> cartesian_points;

    // values is an evaluated function on grid, value in i'th position
    // corresponds to the i'th point in cartesian_points vector
    std::vector<double> values;

    // trapeziums are the elementary objects of which the grid consists. They
    // are stored as vectors of columns of trapeziums on a sphere.
    std::vector<std::vector<Trapezium>> trapeziums;

    void GenerateGridAndEvaluateFunc();
    void EvaluateFuncRoutine(int lower_bound, int upper_bound);

    std::tuple<int, ClassicGridPoint> GetNextTrapeziumIndices(int prev_side, const ClassicGridPoint& prev_trapezium_indices); // TODO TrapeziumIndices -> next point and add explanation about moving using trapeziums
    ClassicGridPoint GetUpperTrapeziumIndices(const ClassicGridPoint& prev_trapezium_indices);
    ClassicGridPoint GetLowerTrapeziumIndices(const ClassicGridPoint& prev_trapezium_indices);
    ClassicGridPoint GetLeftTrapeziumIndices(const ClassicGridPoint& prev_trapezium_indices);
    ClassicGridPoint GetRightTrapeziumIndices(const ClassicGridPoint& prev_trapezium_indices);
    CoordsOfPoint FindIntersection(
        const CoordsOfPoint& x1, const CoordsOfPoint& x2, double a, double b, double c);
    std::tuple<int, CoordsOfPoint> ProcessSegment(const Trapezium& trapezium,
                                                  int index,
                                                  int start_side,
                                                  double iso_value);
    int GetTrapeziumIndex(const Trapezium& trapezium, double value);

public:
    AngularGauss* distr;

    ClassicGrid(int, AngularGauss*);
    CoordsOfPoint GetCoordsOfPoint(const AnglesOfPoint& point);
    AnglesOfPoint GetAnglesOfPoint(const CoordsOfPoint& point);
    std::vector<CoordsOfPoint> GetIsolineCoords(double iso_value);
};

#endif // CLASSIC_GRID_HPP
