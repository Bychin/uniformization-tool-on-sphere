#ifndef STATS_HPP
#define STATS_HPP

#include <array>
#include <vector>

#include "nlohmann/json.hpp"

#include <grids/classic_grid.hpp>
#include <grids/spiral_grid.hpp>
#include <types/types.hpp>

using json = nlohmann::json;

class StatsAPI {
    ClassicGrid* classic_grid;
    SpiralGrid* spiral_grid;

    std::vector<CoordsOfPoint> points;
    std::vector<double> alpha95; // i'th coeff correspons to i'th point

    double DistanceBetweenPoints(const CoordsOfPoint&, const CoordsOfPoint&) const;
    int InsertPointIntoIsoline(const CoordsOfPoint&, IsolineCoords&) const;
    Vector GetNormalizedVector(const CoordsOfPoint&, const CoordsOfPoint&) const;
    int GetIndexOfPrevIsolinePoint(int, const IsolineCoords&) const;
    int GetIndexOfNextIsolinePoint(int, const IsolineCoords&) const;
    bool CheckClockwiseDirection(const Vector&, const Vector&, const Vector&) const;

    double CalcGradientModulus(const CoordsOfPoint&, const CoordsOfPoint&) const;
    // Calculates line integral on infinitesimal straight line between two points
    double CalcIntegralOnInfinitesimalCurve(const CoordsOfPoint&, const CoordsOfPoint&) const;
    double CalcIntegralOnFullIsoline(const IsolineCoords&) const;
    double CalcIntegralOnLeftCurveOnIsoline(const IsolineCoords&, int, int) const;
    double CalcIntegralOnRightCurveOnIsoline(const IsolineCoords&, int, int) const;

public:
    StatsAPI(std::vector<CoordsOfPoint>& points, ClassicGrid* classic_grid, SpiralGrid* spiral_grid);
    StatsAPI(std::vector<CoordsOfPoint>& points, std::vector<double>& alpha95Coeffs, ClassicGrid* classic_grid, SpiralGrid* spiral_grid);

    const std::vector<CoordsOfPoint> Points(void) const;
    const std::vector<double> Alpha95Coeffs(void) const;

    void Validate(void) const;
    Stats CalculateStats();
    StatsWithDebugInfo CalculateStatsWithDebugInfo();

    // TODO make private
    double CalculateFunc(const CoordsOfPoint& point) const; // TODO remove?
    double CalculateTStat(double value) const;
    double CalculateSStat(const CoordsOfPoint& point, double value) const;
    json CalculateSStatWithDebugInfo(const CoordsOfPoint& point, double value) const;
};

#endif // STATS_HPP
