#include "von_mises_fisher.hpp"

#include "cmath"

#include "types/types.hpp"

VonMisesFisher::VonMisesFisher(double kappa) : kappa(kappa), mean(kMeanDir) {
    calcCoeff = kappa / (2 * M_PI * (std::exp(kappa) - std::exp(-kappa)));
};

double VonMisesFisher::Calc(const CoordsOfPoint& u) const {
    // actually it can be reduced to mean[0] * u[0] due to kMeanDir
    double dot_product = mean[0] * u[0] + mean[1] * u[1] + mean[2] * u[2];

    return calcCoeff * std::exp(kappa * dot_product);
}
