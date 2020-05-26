#ifndef SERVER_TYPES_HPP
#define SERVER_TYPES_HPP

#include <array>
#include <vector>

// ClassicGridPoint represents (i, j) pair of point indices on a classic grid
typedef std::array<int, 2> ClassicGridPoint;
// Trapezium is a simple element that makes up a classic grid
typedef std::array<ClassicGridPoint, 4> Trapezium;

// AnglesOfPoint is a 3D point in spherical coordinates
typedef std::array<double, 2> AnglesOfPoint;
// CoordsOfPoint is a 3D point in Cartesian coordinates
typedef std::array<double, 3> CoordsOfPoint;
// Vector in 3D, starting from origin
typedef std::array<double, 3> Vector;

// List of points in Cartesian coordinates making up isoline coordinates
typedef std::vector<CoordsOfPoint> IsolineCoords;

#endif // SERVER_TYPES_HPP
