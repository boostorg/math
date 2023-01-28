//  (C) Copyright Nick Thompson, John Maddock 2023.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#include <random>
#include <numeric>
#include <boost/math/tools/estrin.hpp>
#include "math_unit_test.hpp"
#if BOOST_MATH_TEST_FLOAT128
#include <boost/multiprecision/float128.hpp>
#endif

using boost::math::tools::estrin;

template<typename Real, typename Complex>
void test_polynomial_properties(size_t seed = 0) {
   if (seed == 0) {
     std::random_device rd;
     seed = rd();
   }
   std::mt19937_64 gen(seed);
   std::uniform_real_distribution<Real> dis(-1, 1);
   size_t n = seed % 128;
   std::vector<Real> coeffs(n);
   for (auto& c : coeffs) {
        c = dis(gen);
   }
   // Test evaluation at zero:
   auto p0 = estrin(coeffs, Complex(0));
   if(!CHECK_ULP_CLOSE(p0.real(), coeffs[0], 0)) {
       std::cerr << "  p(0) != c0 with seed " << seed << "\n";  
   }
   // Just given the structure of the calculation, I believe the complex part should be identically zero:
   if (!CHECK_ULP_CLOSE(p0.imag(), Real(0), 0)) {
       std::cerr << "  p(0) != c0 with seed " << seed << "\n";
   }

   auto p1_computed = estrin(coeffs, Complex(1));
   // Might want to consider using Kahan summation so the expected value is super accurate:
   auto p1_expected = std::reduce(coeffs.begin(), coeffs.end());
   if (!CHECK_ULP_CLOSE(p1_expected, p1_computed.real(), coeffs.size())) {
       std::cerr << " p(1) != sum(coeffs) with seed " << seed << "\n";
   }
   if (!CHECK_ULP_CLOSE(Real(0), p1_computed.imag(), coeffs.size())) {
       std::cerr << " p(1) != sum(coeffs) with seed " << seed << "\n";
   }
}

template<typename Real>
void test_std_array_overload(size_t seed = 0) {
   if (seed == 0) {
     std::random_device rd;
     seed = rd();
   }
   std::mt19937_64 gen(seed);
   std::uniform_real_distribution<Real> dis(-1, 1);
   std::array<Real, 12> coeffs;
   for (auto& c : coeffs) {
        c = dis(gen);
   }
   // Test evaluation at zero:
   auto p0 = estrin(coeffs, Real(0));
   if(!CHECK_ULP_CLOSE(p0, coeffs[0], 0)) {
       std::cerr << "  p(0) != c0 with seed " << seed << "\n";  
   }
   auto p0_real = estrin(coeffs, Real(0));
   if(!CHECK_ULP_CLOSE(p0_real, coeffs[0], 0)) {
       std::cerr << "  p(0) != c0 with seed " << seed << "\n";  
   }

   auto p1_computed = estrin(coeffs, Real(1));
   auto p1_expected = std::reduce(coeffs.begin(), coeffs.end());
   if (!CHECK_ULP_CLOSE(p1_expected, p1_computed, coeffs.size())) {
       std::cerr << " p(1) != sum(coeffs) with seed " << seed << "\n";
   }
}


int main() {
   test_polynomial_properties<float, std::complex<float>>();
   test_polynomial_properties<double, std::complex<double>>();
   test_std_array_overload<float>();
   test_std_array_overload<double>();
#if BOOST_MATH_TEST_FLOAT128
   test_std_array_overload<boost::multiprecision::float128>();
#endif
   return boost::math::test::report_errors();
}
