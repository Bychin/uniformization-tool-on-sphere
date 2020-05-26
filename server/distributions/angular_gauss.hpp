#ifndef ANGULAR_GAUSS_HPP
#define ANGULAR_GAUSS_HPP

#include <array>

#include <boost/numeric/ublas/matrix.hpp>

#include "types/types.hpp"

namespace ublas = boost::numeric::ublas;

// Angular Gaussian distribution
class AngularGauss {
    Vector mean;
    double mean_norm;

    ublas::matrix<double> lambda;
    double det_lambda;

    // InvMatrix calculates inverse matrix for input_mat and stores it in
    // inverse_mat
    void InvMatrix(const ublas::matrix<double>& input_mat, ublas::matrix<double>& inverse_mat);
    // Det calculates m's determinant
    double Det(const ublas::matrix<double>& m) const;
    // InnerProduct calculates inner product of x and y using lambda
    double InnerProduct(const Vector& x, const Vector& y) const;

public:
    // AngularGauss creates new Angular Gaussian distribution with mean using
    // mean_vec and covariance using cov_mat
    AngularGauss(const Vector& mean_vec, const ublas::matrix<double>& cov_mat);

    // Mean returns mean vector for Angular Gaussian distribution
    const Vector Mean(void) const;

    // Calc calculates the Angular Gaussian PDF at u
    double Calc(const CoordsOfPoint& u) const;
};

#endif // ANGULAR_GAUSS_HPP
