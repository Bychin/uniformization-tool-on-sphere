#include "util.hpp"

#include <cassert>
#include <cmath>
#include <vector>

std::vector<int> Bounds(int parts_amount, int points_amount) {
    assert(parts_amount <= points_amount);

    int part_max_size = std::round(double(points_amount) / double(parts_amount));

    auto bounds = std::vector<int>(parts_amount+1);
    bounds[0] = 0;
    for (int i = 1; i < parts_amount+1; ++i)
        bounds[i] = bounds[i-1] + part_max_size;

    bounds[parts_amount] = points_amount;

    return bounds;
}
