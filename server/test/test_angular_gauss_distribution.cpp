#include <array>
#include <boost/numeric/ublas/matrix.hpp>
#include <iostream>

#include "distributions/angular_gauss.hpp"

namespace ublas = boost::numeric::ublas;

// Simple check function. If parameter is not true, test fails.
#define IS_TRUE(x)                                                                    \
    {                                                                                 \
        if (!x)                                                                       \
            std::cout << __FUNCTION__ << " failed on line " << __LINE__ << std::endl; \
    }

// Function to test
bool function1() {
    std::array<double, 3> mean_vec = {1., 1., 1.};

    auto cov_values = {1.5, 1., 1., 0., 1., 0., 0., 0., 1.};
    ublas::unbounded_array<double> cov_storage(3 * 3);
    std::copy(cov_values.begin(), cov_values.end(), cov_storage.begin());
    ublas::matrix<double> cov_mat(3, 3, cov_storage);

    auto g = AngularGauss(mean_vec, cov_mat);

    std::array<double, 3> u_vec = {1., 1., 1.};
    std::cout << g.Calc(u_vec) << std::endl;

    return true;
}

// Test for function1()
// You would need to write these even when using a framework
void test_function1() {
    IS_TRUE(function1());
    IS_TRUE(function1());
    IS_TRUE(function1());
}

int main(void) {
    // Call all tests. Using a test framework would simplify this.
    test_function1();
}
