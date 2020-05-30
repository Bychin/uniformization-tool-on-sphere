#ifndef CLASSIC_GRID_HPP
#define CLASSIC_GRID_HPP

#include <array>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <boost/container_hash/hash.hpp>

#include "distributions/angular_gauss.hpp"
#include "distributions/von_mises_fisher.hpp"
#include "types/types.hpp"

class ClassicGrid {
    int div;
    AngularGauss* distr;
    VonMisesFisher* fisher_distr;

    // points is a map of (i, j) pairs, that are indices of grid points, to
    // the indicies of the same grid's points in Cartesian coordinates'.
    //
    // Note, that these (i, j) pairs could NOT be used directly for
    // spherical coordinates.
    std::unordered_map<ClassicGridPoint, int, boost::hash<ClassicGridPoint>> points;

    // cartesian_points is a vector of grid's points in Cartesian coordinates
    std::vector<CoordsOfPoint> cartesian_points;

    // values is an evaluated function on grid, value in i'th position
    // corresponds to the function's value of i'th point in cartesian_points
    // vector
    std::vector<double> values;
    // TODO
    std::vector<double> fisher_values;
    // TODO
    std::vector<double> values_after_convolution;

    // trapeziums are the elementary objects of which the grid consists. They
    // are stored as vectors of columns of trapeziums on a sphere.
    std::vector<std::vector<Trapezium>> trapeziums;

    void GenerateGridAndEvaluateFunc(void);
    void EvaluateFuncRoutine(int lower_bound, int upper_bound);
    void EvaluateFisherFuncRoutine(int lower_bound, int upper_bound);

    TrapeziumIndex GetTrapeziumIndex(const CoordsOfPoint& u) const;
    TrapeziumIndex FindTrapeziumWithIsolineInside(double iso_value);

    std::tuple<int, TrapeziumIndex> GetNextTrapeziumIndex(int previous_side, const TrapeziumIndex& previous) const; // TODO remove TrapeziumIndex
    TrapeziumIndex GetUpperTrapeziumIndex(const TrapeziumIndex& previous) const;
    TrapeziumIndex GetLowerTrapeziumIndex(const TrapeziumIndex& previous) const;
    TrapeziumIndex GetLeftTrapeziumIndex(const TrapeziumIndex& previous) const;
    TrapeziumIndex GetRightTrapeziumIndex(const TrapeziumIndex& previous) const;

    // FindIntersection uses linear interpolation to find an intersection point
    // with value c between vertices with values a and b, where x1 = f^(-1)(a),
    // x2 = f^(-1)(b), sign((a - c) * (b - c)) must be -1.
    CoordsOfPoint FindIntersection(const CoordsOfPoint& x1, const CoordsOfPoint& x2, double a, double b, double c) const;
    // GetVector creates new vector from x1 to x2
    Vector GetVector(const CoordsOfPoint& x1, const CoordsOfPoint& x2) const;
    // TODO
    double InterpolateInsideABC(const CoordsOfPoint& u, const ClassicGridPoint& A, const ClassicGridPoint& B, const ClassicGridPoint& C);

    // ProcessSegment processes trapezium using marching square algorithm and
    // returns last processed trapezium side and point.
    std::tuple<int, CoordsOfPoint> ProcessSegment(const Trapezium& trapezium, int index, int start_side, double iso_value);

    // GetMarchingSquareIndex returns the trapezium index needed for marching squares algorithm.
    // See Lookup Table here: https://en.wikipedia.org/wiki/Marching_squares#Basic_algorithm
    int GetMarchingSquareIndex(const Trapezium& trapezium, double value);

public:
    // ClassicGrid creates new classic grid with bandwidth grid_div and
    // function distr on it
    ClassicGrid(int grid_div, AngularGauss* distr);

    // TODO
    void SetFisherKappa(double kappa);
    double CalcFunc(const CoordsOfPoint& u);

    // Func returns function on grid
    const AngularGauss* Func(void) const; // TODO rename? // TODO need only for Mean()

    // GetCoordsOfPoint converts spherical coordinates of a point to Cartesian
    static CoordsOfPoint GetCoordsOfPoint(const AnglesOfPoint& point);
    // GetAnglesOfPoint converts Cartesian coordinates of a point to spherical
    static AnglesOfPoint GetAnglesOfPoint(const CoordsOfPoint& point);

    // GetIsolineCoords returns the coordinates of the points for the function
    // isoline with a value equal to iso_value.
    // It uses marching squares algorithm, see more info:
    // https://en.wikipedia.org/wiki/Marching_squares
    std::vector<CoordsOfPoint> GetIsolineCoords(double iso_value);
};

#endif // CLASSIC_GRID_HPP
