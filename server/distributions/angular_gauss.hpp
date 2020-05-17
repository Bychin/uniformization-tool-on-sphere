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

    void InvMatrix(ublas::matrix<double>&, ublas::matrix<double>&);
    double Det(ublas::matrix<double>&);
    double InnerProduct(std::array<double, 3>&, std::array<double, 3>&);

public:
    AngularGauss(std::array<double, 3>&, ublas::matrix<double>&);
    double Calc(std::array<double, 3>&);
    std::array<double, 3> Mean();
};

#endif // ANGULAR_GAUSS_HPP
