#ifndef BESSEL_HPP
#define BESSEL_HPP

#include <cmath>                    // For std::sph_bessel and other math functions
#include <boost/math/special_functions/factorials.hpp> // For Boost double factorial
#include <source_location>          // For enhanced debugging (C++23)

// Computes the double factorial using Boost
inline double DoubleFactorial(int l) {
    return boost::math::double_factorial<double>(l);
}

// Computes spherical Bessel function using std::sph_bessel (C++17+)
constexpr double bessel_j(int l, double x) {
    return std::sph_bessel(l, x);
}

// Computes spherical Bessel function using downward recursion
inline double Bessel_j(int l, double x) {
    if (x > l) return std::sph_bessel(l, x); // Use std::sph_bessel for large x
    if (x < 1e-20) return 0; // Avoid division by zero
    int lstart = l + static_cast<int>(sqrt(40.0 * l) / 2.0); // Estimate starting point
    double j2 = 0, j1 = 1, j0, jl, x1 = 1 / x; // Store 1/x for performance
    for (int i = lstart; i >= 0; i--) {
        j0 = (2 * i + 3.0) * x1 * j1 - j2;
        if (i == l) jl = j0;
        j2 = j1;
        j1 = j0;
    }
    double true_j0 = std::sph_bessel(0, x);
    return jl * true_j0 / j0; // Renormalize result
}

// Computes the derivative of the spherical Bessel function
inline double dbessel_j(int l, double x) {
    if (l == 0 && fabs(x) < 1e-20) return 0;
    if (fabs(x) < 1e-5) return l * pow(x, l - 1) / DoubleFactorial(2 * l + 1);
    return l * std::sph_bessel(l, x) / x - std::sph_bessel(l + 1, x);
}

// Computes x * d/dx log(j_l(x)) and j_l(x)
inline void dlog_bessel_j(int l, double x, double &dlogj, double &jl, double &jldlogj) {
    if (fabs(x) < 1e-5) {
        dlogj = l;
        jl = pow(x, l) / DoubleFactorial(2 * l + 1);
        jldlogj = jl * dlogj;
        return;
    }
    jl = std::sph_bessel(l, x);
    jldlogj = l * jl - x * std::sph_bessel(l + 1, x);
    dlogj = jldlogj / jl;
}

// Computes the Legendre polynomial P_l(x) using std::legendre (C++17+)
constexpr double Legendre(int l, double x) {
    return std::legendre(l, x);
}

#endif // BESSEL_HPP
