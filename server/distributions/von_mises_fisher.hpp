#ifndef VON_MISES_FISHER_HPP
#define VON_MISES_FISHER_HPP

#include "types/types.hpp"

const Vector kMeanDir = {0, 0, 1};

// von Mises–Fisher distribution
class VonMisesFisher {
    Vector mean;
    double kappa;

    double calcCoeff;

public:
    // VonMisesFisher creates new von Mises–Fisher distribution with mean
    // direction in kMeanDir and concentration parameter kappa
    VonMisesFisher(double kappa);

    // Calc calculates the von Mises–Fisher PDF at u
    double Calc(const CoordsOfPoint& u) const;
};

#endif // VON_MISES_FISHER_HPP
