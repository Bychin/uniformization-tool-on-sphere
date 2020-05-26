#ifndef STATS_HPP
#define STATS_HPP

#include <array>
#include <vector>

#include <grids/classic_grid.hpp>
#include <grids/spiral_grid.hpp>
#include <types/types.hpp>

class StatsAPI {
    std::vector<CoordsOfPoint> points;
    ClassicGrid* classic_grid;
    SpiralGrid* spiral_grid;

    double DistanceBetweenPoints(const CoordsOfPoint&, const CoordsOfPoint&) const;
    int InsertPointIntoIsoline(const CoordsOfPoint&, IsolineCoords&) const;
    Vector GetNormalizedVector(const CoordsOfPoint&, const CoordsOfPoint&) const;
    int GetIndexOfPrevIsolinePoint(int, const IsolineCoords&) const;
    int GetIndexOfNextIsolinePoint(int, const IsolineCoords&) const;
    bool CheckClockwiseDirection(const Vector&, const Vector&, const Vector&) const;

    // Calculates line integral on infinitesimal straight line between two points
    double CalcIntegralOnInfinitesimalCurve(const CoordsOfPoint&, const CoordsOfPoint&) const;
    double CalcIntegralOnFullIsoline(const IsolineCoords&) const;
    double CalcIntegralOnLeftCurveOnIsoline(const IsolineCoords&, int, int) const;
    double CalcIntegralOnRightCurveOnIsoline(const IsolineCoords&, int, int) const;

public:
    StatsAPI(std::vector<CoordsOfPoint>& points, ClassicGrid* classic_grid, SpiralGrid* spiral_grid);

    const std::vector<CoordsOfPoint> Points(void) const;

    void Validate(void) const;
    double CalculateFunc(const CoordsOfPoint& point) const;
    double CalculateTStat(double value) const;
    double CalculateSStat(const CoordsOfPoint& point, double value) const;
};

#endif // STATS_HPP
