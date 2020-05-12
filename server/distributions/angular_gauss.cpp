#include "angular_gauss.hpp"

#include <array>
#include <cmath>

#include <boost/math/special_functions/erf.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/matrix.hpp>

namespace ublas = boost::numeric::ublas;

void AngularGauss::InvMatrix(ublas::matrix<double>& input_mat, ublas::matrix<double>& inverse_mat) {
    ublas::matrix<double> A(input_mat);

    ublas::permutation_matrix<std::size_t> pm(A.size1());
    int factorization_status = ublas::lu_factorize(A, pm);
    assert(factorization_status == 0);

    inverse_mat.assign(ublas::identity_matrix<double>(A.size1()));
    lu_substitute(A, pm, inverse_mat);
}

double AngularGauss::Det(ublas::matrix<double>& m) {
    const double a = m(0, 0);
    const double b = m(0, 1);
    const double c = m(0, 2);
    const double d = m(1, 0);
    const double e = m(1, 1);
    const double f = m(1, 2);
    const double g = m(2, 0);
    const double h = m(2, 1);
    const double k = m(2, 2);

    const double determinant =
        (a * ((e * k) - (f * h))) - (b * ((k * d) - (f * g))) + (c * ((d * h) - (e * g)));

    return determinant;
}

double AngularGauss::InnerProduct(std::array<double, 3>& x, std::array<double, 3>& y) {
    double result = 0;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            result += lambda(i, j) * x.at(i) * y.at(j);

    return result;
}

AngularGauss::AngularGauss(std::array<double, 3>& mean_vec, ublas::matrix<double>& cov_mat) {
    assert(cov_mat.size1() == 3 && cov_mat.size2() == 3);

    std::copy(mean_vec.begin(), mean_vec.end(), mean.begin());
    lambda = ublas::matrix<double>(3,3);
    InvMatrix(cov_mat, lambda);
    det_lambda = Det(lambda);
    mean_norm  = std::sqrt(InnerProduct(mean, mean));
}

double AngularGauss::Calc(std::array<double, 3>& u) {
    double u_norm = std::sqrt(InnerProduct(u, u));
    double z      = InnerProduct(mean, u) / u_norm;

    double coeff = std::exp(-0.5 * mean_norm * mean_norm) * std::sqrt(det_lambda) /
                   (4 * M_PI * u_norm * u_norm * u_norm);

    return coeff * (z * std::sqrt(M_2_PI) +
                    std::exp(0.5 * z * z) * (1 + z * z) * (1 + boost::math::erf(z / M_SQRT2)));
}
