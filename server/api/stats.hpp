#ifndef STATS_HPP
#define STATS_HPP

#include <array>
#include <vector>

#include <grids/classic_grid.hpp>
#include <grids/spiral_grid.hpp>

class StatsAPI {
    double DistanceBetweenPoints(std::array<double, 3>&, std::array<double, 3>&);
    int InsertPointIntoIsoline(std::array<double, 3>&, std::vector<std::array<double, 3>>&);
    std::array<double, 3> GetNormalizedVector(std::array<double, 3>&, std::array<double, 3>&);
    int GetIndexOfPrevIsolinePoint(int, std::vector<std::array<double, 3>>&);
    int GetIndexOfNextIsolinePoint(int, std::vector<std::array<double, 3>>&);
    bool CheckClockwiseDirection(std::array<double, 3>&, std::array<double, 3>&, std::array<double, 3>&);
    double CalcIntegralOnInfinitesimalCurve(std::array<double, 3>&, std::array<double, 3>&);
    double CalcIntegralOnFullIsoline(std::vector<std::array<double, 3>>&);
    double CalcIntegralOnLeftCurveOnIsoline(std::vector<std::array<double, 3>>&, int, int);
    double CalcIntegralOnRightCurveOnIsoline(std::vector<std::array<double, 3>>&, int, int);

public:
    std::vector<std::array<double, 3>> points;
    ClassicGrid* classic_grid;
    SpiralGrid* spiral_grid;

    StatsAPI(std::vector<std::array<double, 3>>& points,
             ClassicGrid* classic_grid,
             SpiralGrid* spiral_grid);
    // TODO bool Validate();
    double CalculateTStat(double);
    double CalculateSStat(std::array<double, 3>&, double);
};

#endif // STATS_HPP
