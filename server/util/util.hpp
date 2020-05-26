#ifndef UTIL_HPP
#define UTIL_HPP

#include <vector>

// Bounds splits points_amount into parts_amount and returns a list of bounds,
// for example: [b1, b2, b3], where b2 is the end of the first bound and the
// start of second bound.
std::vector<int> Bounds(int parts_amount, int points_amount);

#endif // UTIL_HPP
