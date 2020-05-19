#ifndef ANGULAR_GAUSS_HPP
#define ANGULAR_GAUSS_HPP

#include <array>

#include <boost/numeric/ublas/matrix.hpp>

namespace ublas = boost::numeric::ublas;

// Angular Gaussian distribution
class AngularGauss {
    std::array<double, 3> mean;
    double mean_norm;

    ublas::matrix<double> lambda;
    double det_lambda;

    void InvMatrix(const ublas::matrix<double>&, ublas::matrix<double>&);
    double Det(const ublas::matrix<double>&);
    double InnerProduct(const std::array<double, 3>&, const std::array<double, 3>&);

public:
    AngularGauss(const std::array<double, 3>&, const ublas::matrix<double>&);
    double Calc(const std::array<double, 3>&);
    std::array<double, 3> Mean();
};

#endif // ANGULAR_GAUSS_HPP
