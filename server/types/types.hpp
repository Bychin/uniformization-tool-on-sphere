#ifndef SERVER_TYPES_HPP
#define SERVER_TYPES_HPP

#include <array>
#include <tuple>
#include <vector>

#include "nlohmann/json.hpp"

// ClassicGridPoint represents (i, j) pair of point indices on a classic grid
typedef std::array<int, 2> ClassicGridPoint;
// Trapezium is a simple element that makes up a classic grid
typedef std::array<ClassicGridPoint, 4> Trapezium;
// TrapeziumIndex represents (i, j) pair as index of trapezium on a classic grid
typedef std::array<int, 2> TrapeziumIndex;

// AnglesOfPoint is a 3D point in spherical coordinates,
// (phi, theta) pair, where phi - azimuthal angle, theta - polar angle
typedef std::array<double, 2> AnglesOfPoint;
// CoordsOfPoint is a 3D point in Cartesian coordinates
typedef std::array<double, 3> CoordsOfPoint;
// Vector in 3D, starting from origin
typedef std::array<double, 3> Vector;

// List of points in Cartesian coordinates making up isoline coordinates
typedef std::vector<CoordsOfPoint> IsolineCoords;

// t and s statistics for directional data
typedef std::tuple<std::vector<double>, std::vector<double>> Stats;
// t, s statistics and extra info for s statistic for directional data
typedef std::tuple<std::vector<double>, std::vector<double>, nlohmann::json> StatsWithDebugInfo;

#endif // SERVER_TYPES_HPP
