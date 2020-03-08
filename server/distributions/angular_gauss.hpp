#ifndef ANG_GAUSSIAN_DISTRIBUTION_HPP
#define ANG_GAUSSIAN_DISTRIBUTION_HPP

#include <boost/numeric/ublas/matrix.hpp>

namespace ublas = boost::numeric::ublas;

// Angular Gaussian distribution
class AngularGauss {
    ublas::vector<double> mean;
    double mean_norm;

    ublas::matrix<double> lambda;
    double det_lambda;

    void InvMatrix(ublas::matrix<double>&, ublas::matrix<double>&);
    double Det(ublas::matrix<double>&);
    double InnerProduct(ublas::vector<double>&, ublas::vector<double>&);

public:
    AngularGauss(ublas::vector<double>&, ublas::matrix<double>&);
    double Calc(ublas::vector<double>&);
};

#endif // ANG_GAUSSIAN_DISTRIBUTION_HPP
