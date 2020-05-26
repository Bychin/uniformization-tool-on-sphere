#ifndef PSVT_STATISTICS_HPP
#define PSVT_STATISTICS_HPP

#include <vector>

// UnifTestResults contains the results of two uniformity tests:
// Kolmogorov-Smirnov test and Anderson-Darling test
struct UnifTestResults {
    double KS_measure;
    double KS_estimate;
    double AD_measure;
    double AD_estimate;
};

// UniformityTests performs for sample KS and AD tests for uniformity
UnifTestResults UniformityTests(const std::vector<double>& sample);

#endif // PSVT_STATISTICS_HPP
