/*
 * Copyright Nick Thompson, 2026
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 * Module build check for <boost/math/tools/fft.hpp>. It exercises fft and real_fft through
 * the boost.math module with std::vector and std::array of std::complex<double> (and of double
 * for real_fft), on the three code paths (powers of 2, mixed radix, Bluestein).
 * Multiprecision does not work in module mode, so the multiprecision cases of test_fft.cpp are
 * dropped and the reference is a naive DFT in long double.
 */

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/tools/fft.hpp>
#else
import boost.math;
#endif

#ifndef BOOST_MATH_BUILD_MODULE
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <stdexcept>
#include <vector>
#endif
#include "math_unit_test.hpp"

using boost::math::tools::fft;
using boost::math::tools::real_fft;
using complex = std::complex<double>;

std::vector<std::complex<long double>> naive_dft(std::vector<complex> const& x)
{
   std::size_t n = x.size();
   std::vector<std::complex<long double>> X(n);
   for (std::size_t k = 0; k < n; ++k)
   {
      std::complex<long double> s(0, 0);
      for (std::size_t j = 0; j < n; ++j)
      {
         long double theta = -2 * 3.14159265358979323846264338327950288L * static_cast<long double>((j * k) % n) / static_cast<long double>(n);
         s += std::complex<long double>(x[j].real(), x[j].imag()) * std::complex<long double>(std::cos(theta), std::sin(theta));
      }
      X[k] = s;
   }
   return X;
}

template <class Container>
void check_spectrum(Container const& y, std::vector<std::complex<long double>> const& X, std::size_t count, std::size_t n)
{
   long double err2 = 0;
   long double norm2 = 0;
   for (std::size_t k = 0; k < count; ++k)
   {
      err2 += std::norm(std::complex<long double>(y[k].real(), y[k].imag()) - X[k]);
      norm2 += std::norm(X[k]);
   }
   double L = n < 2 ? 1.0 : std::log2(double(n));
   // Allow for the reference's own error where long double is no wider than double.
   double reference_error = double(4 * std::numeric_limits<long double>::epsilon()) * std::sqrt(double(n));
   CHECK_LE(double(std::sqrt(err2 / norm2)), 4 * std::numeric_limits<double>::epsilon() * std::sqrt(L) + reference_error);
}

template <class Container>
void check_round_trip(Container const& y, std::vector<complex> const& x)
{
   double L = x.size() < 2 ? 1.0 : std::log2(double(x.size()));
   for (std::size_t j = 0; j < x.size(); ++j)
   {
      CHECK_ABSOLUTE_ERROR(x[j].real(), y[j].real(), 8 * std::numeric_limits<double>::epsilon() * L);
      CHECK_ABSOLUTE_ERROR(x[j].imag(), y[j].imag(), 8 * std::numeric_limits<double>::epsilon() * L);
   }
}

template <std::size_t N>
void test_containers()
{
   std::vector<complex> x(N);
   for (std::size_t j = 0; j < N; ++j)
   {
      x[j] = complex(std::cos(1.7 * double(j)), std::sin(0.3 * double(j * j)));
   }
   auto X = naive_dft(x);

   fft<complex> f(N);
   std::vector<complex> v(x);
   std::array<complex, N> a;
   std::copy(x.begin(), x.end(), a.begin());
   f.transform(v);
   f.transform(a);
   check_spectrum(v, X, N, N);
   check_spectrum(a, X, N, N);
   f.inverse_transform(v);
   f.inverse_transform(a);
   check_round_trip(v, x);
   check_round_trip(a, x);

   std::vector<complex> xr(N);
   std::vector<double> rv(N);
   std::array<double, N> ra;
   for (std::size_t j = 0; j < N; ++j)
   {
      xr[j] = complex(x[j].real(), 0);
      rv[j] = x[j].real();
      ra[j] = rv[j];
   }
   auto XR = naive_dft(xr);
   constexpr std::size_t M = N / 2 + 1;
   real_fft<complex> rf(N);
   std::vector<complex> out_v(M);
   std::array<complex, M> out_a;
   rf.transform(rv, out_v);
   check_spectrum(out_v, XR, M, N);
   rf.transform(ra, out_a);
   check_spectrum(out_a, XR, M, N);
   std::array<complex, M + 1> wrong;
   CHECK_THROW(rf.transform(rv, wrong), std::domain_error);
}

int main()
{
   test_containers<1>();
   test_containers<2>();
   test_containers<16>();
   test_containers<12>();
   test_containers<45>();
   test_containers<17>();
   test_containers<38>();
   test_containers<97>();
   CHECK_THROW(fft<complex>(0), std::domain_error);
   return boost::math::test::report_errors();
}
