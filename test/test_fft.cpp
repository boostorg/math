/*
 * Copyright Nick Thompson, 2026
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#include <boost/math/tools/fft.hpp>
#include "math_unit_test.hpp"
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_complex.hpp>
#include <array>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <new>
#include <random>
#include <stdexcept>
#include <vector>

// Counts every allocation, so that the tests can check that transforms don't allocate after construction.
static std::size_t allocations = 0;

void* operator new(std::size_t size)
{
   ++allocations;
   if (void* p = std::malloc(size == 0 ? 1 : size))
   {
      return p;
   }
   throw std::bad_alloc();
}

void operator delete(void* p) noexcept
{
   std::free(p);
}

void operator delete(void* p, std::size_t) noexcept
{
   std::free(p);
}

using boost::math::tools::fft;
using boost::math::tools::real_fft;
using boost::multiprecision::cpp_bin_float_50;
using boost::multiprecision::cpp_bin_float_100;
using boost::multiprecision::cpp_complex_50;
using boost::multiprecision::cpp_complex_100;
using boost::multiprecision::cpp_complex_quad;

template <class Complex>
using real_t = typename fft<Complex>::real_type;

template <class Complex>
real_t<Complex> squared_modulus(Complex const& z)
{
   return z.real() * z.real() + z.imag() * z.imag();
}

inline double log2_length(std::size_t n)
{
   return n < 2 ? 1.0 : std::log2(double(n));
}

// O(n^2) DFT in the precision of RefComplex, with twiddles computed directly from the reduced angle 2*pi*(jk mod n)/n.
template <class RefComplex>
std::vector<RefComplex> naive_dft(std::vector<RefComplex> const& x)
{
   using Ref = real_t<RefComplex>;
   using std::cos;
   using std::sin;
   std::size_t n = x.size();
   std::vector<Ref> wr(n), wi(n);
   Ref two_pi = 2 * boost::math::constants::pi<Ref>();
   for (std::size_t k = 0; k < n; ++k)
   {
      Ref theta = two_pi * Ref(k) / Ref(n);
      wr[k] = cos(theta);
      wi[k] = -sin(theta);
   }
   std::vector<RefComplex> X(n);
   for (std::size_t k = 0; k < n; ++k)
   {
      Ref re = 0;
      Ref im = 0;
      for (std::size_t j = 0; j < n; ++j)
      {
         std::size_t i = (j * k) % n;
         re += x[j].real() * wr[i] - x[j].imag() * wi[i];
         im += x[j].real() * wi[i] + x[j].imag() * wr[i];
      }
      X[k] = RefComplex(re, im);
   }
   return X;
}

// Random input whose components are exactly representable in float, so every type sees the same data.
std::vector<std::complex<double>> random_input(std::size_t n, std::mt19937_64& gen)
{
   std::uniform_real_distribution<float> dis(-1, 1);
   std::vector<std::complex<double>> x(n);
   for (auto& z : x)
   {
      z = std::complex<double>(dis(gen), dis(gen));
   }
   return x;
}

template <class Complex>
Complex convert(std::complex<double> const& z)
{
   using Real = real_t<Complex>;
   return Complex(Real(z.real()), Real(z.imag()));
}

template <class RefComplex>
std::vector<RefComplex> reference(std::vector<std::complex<double>> const& x)
{
   std::vector<RefComplex> xr(x.size());
   for (std::size_t j = 0; j < x.size(); ++j)
   {
      xr[j] = convert<RefComplex>(x[j]);
   }
   return naive_dft(xr);
}

// Normwise relative error O(eps*sqrt(log n)); largest single error relative to the rms of X, O(eps*log n).
template <class Complex, class ComplexContainer, class RefComplex>
bool check_spectrum(ComplexContainer const& y, std::vector<RefComplex> const& X, std::size_t count, std::size_t n)
{
   using Real = real_t<Complex>;
   using Ref = real_t<RefComplex>;
   using std::sqrt;
   Ref err2 = 0;
   Ref norm2 = 0;
   Ref max_err2 = 0;
   for (std::size_t k = 0; k < count; ++k)
   {
      RefComplex d(Ref(y[k].real()) - X[k].real(), Ref(y[k].imag()) - X[k].imag());
      Ref e2 = squared_modulus(d);
      err2 += e2;
      norm2 += squared_modulus(X[k]);
      if (e2 > max_err2)
      {
         max_err2 = e2;
      }
   }
   Real eps = std::numeric_limits<Real>::epsilon();
   Real L = Real(log2_length(n));
   Real rms = Real(sqrt(err2 / norm2));
   Real worst = Real(sqrt(max_err2 / (norm2 / Ref(count))));
   bool ok = CHECK_LE(rms, 2 * eps * sqrt(L));
   ok = CHECK_LE(worst, 3 * eps * L) && ok;
   if (!ok)
   {
      std::cerr << "  n = " << n << "\n";
   }
   return ok;
}

template <class Complex, class ComplexContainer>
void check_round_trip(ComplexContainer const& y, std::vector<std::complex<double>> const& x)
{
   using Real = real_t<Complex>;
   using std::sqrt;
   std::size_t n = x.size();
   Real err2 = 0;
   Real norm2 = 0;
   for (std::size_t j = 0; j < n; ++j)
   {
      Complex xj = convert<Complex>(x[j]);
      err2 += squared_modulus(Complex(y[j].real() - xj.real(), y[j].imag() - xj.imag()));
      norm2 += squared_modulus(xj);
   }
   Real eps = std::numeric_limits<Real>::epsilon();
   if (!CHECK_LE(Real(sqrt(err2 / norm2)), 4 * eps * Real(sqrt(Real(log2_length(n))))))
   {
      std::cerr << "  n = " << n << "\n";
   }
}

template <class Complex, class RefComplex>
void check_against_reference(std::vector<std::complex<double>> const& x, std::vector<RefComplex> const& X)
{
   std::size_t n = x.size();
   std::vector<Complex> y(n);
   for (std::size_t j = 0; j < n; ++j)
   {
      y[j] = convert<Complex>(x[j]);
   }
   fft<Complex> f(n);
   CHECK_EQUAL(f.size(), n);
   f.transform(y);
   check_spectrum<Complex>(y, X, n, n);
   f.inverse_transform(y);
   check_round_trip<Complex>(y, x);
}

void test_against_naive_dft()
{
   std::mt19937_64 gen(17);
   // Powers of 2; mixed-radix lengths, whose prime factors are at most 97 (including 11-, 13- and 17-smooth ones
   // and the primes 89 and 97 just below the cutoff); and Bluestein lengths, with a prime factor above 97
   // (101 and 103 just above it, 127, 509, and composites such as 2*101).
   std::vector<std::size_t> lengths{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 21, 22, 25, 27, 30,
                                    31, 32, 45, 49, 60, 64, 89, 97, 98, 100, 101, 103, 120, 125, 127, 128, 143, 169,
                                    187, 194, 202, 210, 243, 256, 289, 343, 360, 509, 512, 735, 1000, 1001, 1024};
   for (std::size_t n : lengths)
   {
      auto x = random_input(n, gen);
      auto X = reference<cpp_complex_50>(x);
      check_against_reference<std::complex<float>>(x, X);
      check_against_reference<std::complex<double>>(x, X);
      check_against_reference<std::complex<long double>>(x, X);
   }
   for (std::size_t n : {1, 2, 3, 5, 6, 7, 8, 11, 12, 13, 16, 17, 30, 49, 64, 89, 97, 100, 101, 128, 143})
   {
      auto x = random_input(n, gen);
      auto X = reference<cpp_complex_100>(x);
      check_against_reference<cpp_complex_50>(x, X);
      check_against_reference<cpp_complex_quad>(x, X);
   }
}

// Lengths too long for the O(n^2) reference: compare against the FFT in cpp_complex_50,
// whose own accuracy is established above. This checks the log n growth of the error bound.
void test_long_lengths()
{
   std::mt19937_64 gen(29);
   for (std::size_t n : {2018, 2401, 3072, 3125, 4096, 4099, 6208, 6464, 10000, 11011, 12297, 16384, 65536, 100000})
   {
      auto x = random_input(n, gen);
      std::vector<cpp_complex_50> X(n);
      for (std::size_t j = 0; j < n; ++j)
      {
         X[j] = convert<cpp_complex_50>(x[j]);
      }
      fft<cpp_complex_50>(n).transform(X);
      check_against_reference<std::complex<float>>(x, X);
      check_against_reference<std::complex<double>>(x, X);
   }
}

// The container matrix: std::vector and std::array, of std::complex<double> and of cpp_complex_50,
// for C2C, and the same with real inputs for R2C, at a power of 2, a smooth even and a smooth odd composite,
// and an odd and an even length with a large prime factor.
template <class Complex, class RefComplex, std::size_t N>
void test_containers()
{
   using Real = real_t<Complex>;
   std::mt19937_64 gen(N);
   auto x = random_input(N, gen);
   auto X = reference<RefComplex>(x);

   // C2C, std::vector and std::array.
   fft<Complex> f(N);
   std::vector<Complex> v(N);
   std::array<Complex, N> a;
   for (std::size_t j = 0; j < N; ++j)
   {
      v[j] = convert<Complex>(x[j]);
      a[j] = v[j];
   }
   f.transform(v);
   f.transform(a);
   check_spectrum<Complex>(v, X, N, N);
   check_spectrum<Complex>(a, X, N, N);
   f.inverse_transform(v);
   f.inverse_transform(a);
   check_round_trip<Complex>(v, x);
   check_round_trip<Complex>(a, x);

   // R2C: the real parts of x as input, in a vector or an array, into a vector or an array of n/2 + 1.
   std::vector<std::complex<double>> xr(N);
   for (std::size_t j = 0; j < N; ++j)
   {
      xr[j] = std::complex<double>(x[j].real(), 0);
   }
   auto XR = reference<RefComplex>(xr);
   std::vector<Real> rv(N);
   std::array<Real, N> ra;
   for (std::size_t j = 0; j < N; ++j)
   {
      rv[j] = Real(x[j].real());
      ra[j] = rv[j];
   }
   real_fft<Complex> rf(N);
   CHECK_EQUAL(rf.size(), N);
   constexpr std::size_t M = N / 2 + 1;
   std::vector<Complex> out_v(M);
   std::array<Complex, M> out_a;
   rf.transform(rv, out_v);
   check_spectrum<Complex>(out_v, XR, M, N);
   rf.transform(ra, out_a);
   check_spectrum<Complex>(out_a, XR, M, N);
   rf.transform(rv, out_a);
   check_spectrum<Complex>(out_a, XR, M, N);
   rf.transform(ra, out_v);
   check_spectrum<Complex>(out_v, XR, M, N);

   // Round trip through the complex inverse, after filling in X[n-k] = conj(X[k]).
   std::vector<Complex> full(N);
   for (std::size_t k = 0; k < N; ++k)
   {
      Complex z = out_v[k < M ? k : N - k];
      full[k] = k < M ? z : Complex(z.real(), -z.imag());
   }
   f.inverse_transform(full);
   check_round_trip<Complex>(full, xr);

   // Wrong sizes throw.
   std::vector<Complex> too_long(M + 1);
   std::array<Complex, M + 1> too_long_array;
   std::vector<Complex> too_short(N - 1);
   std::vector<Real> wrong_input(N + 1);
   CHECK_THROW(rf.transform(rv, too_long), std::domain_error);
   CHECK_THROW(rf.transform(ra, too_long_array), std::domain_error);
   CHECK_THROW(rf.transform(wrong_input, out_v), std::domain_error);
   CHECK_THROW(f.transform(too_short), std::domain_error);
   CHECK_THROW(f.inverse_transform(too_short), std::domain_error);
}

template <class Complex, class RefComplex>
void test_container_matrix()
{
   test_containers<Complex, RefComplex, 1>();
   test_containers<Complex, RefComplex, 2>();
   test_containers<Complex, RefComplex, 16>();
   test_containers<Complex, RefComplex, 12>();
   test_containers<Complex, RefComplex, 45>();
   test_containers<Complex, RefComplex, 17>();
   test_containers<Complex, RefComplex, 38>();
}

// exp(2*pi*i*p/n), with the argument reduced exactly.
template <class Complex>
Complex unit_root(std::size_t p, std::size_t n)
{
   using Real = real_t<Complex>;
   Real x = Real(2 * (p % n)) / Real(n);
   return Complex(boost::math::cos_pi(x), boost::math::sin_pi(x));
}

template <class Complex>
void test_known_signals()
{
   using Real = real_t<Complex>;
   Real eps = std::numeric_limits<Real>::epsilon();
   for (std::size_t n : {1, 2, 3, 4, 7, 8, 12, 16, 31, 64, 100})
   {
      fft<Complex> f(n);
      Real L = Real(log2_length(n));
      // Delta at 0 -> all ones.
      std::vector<Complex> v(n, Complex(0, 0));
      v[0] = Complex(1, 0);
      f.transform(v);
      for (std::size_t k = 0; k < n; ++k)
      {
         CHECK_ABSOLUTE_ERROR(Real(1), Real(v[k].real()), 4 * eps * L);
         CHECK_ABSOLUTE_ERROR(Real(0), Real(v[k].imag()), 4 * eps * L);
      }
      // Constant -> n at k = 0.
      std::fill(v.begin(), v.end(), Complex(1, 0));
      f.transform(v);
      CHECK_ABSOLUTE_ERROR(Real(n), Real(v[0].real()), 4 * n * eps * L);
      for (std::size_t k = 1; k < n; ++k)
      {
         CHECK_ABSOLUTE_ERROR(Real(0), Real(v[k].real()), 4 * n * eps * L);
         CHECK_ABSOLUTE_ERROR(Real(0), Real(v[k].imag()), 4 * n * eps * L);
      }
      // Single tone exp(2*pi*i*m*j/n) -> n at k = m: this checks the sign convention.
      std::size_t m = n / 3;
      for (std::size_t j = 0; j < n; ++j)
      {
         v[j] = unit_root<Complex>(m * j, n);
      }
      f.transform(v);
      for (std::size_t k = 0; k < n; ++k)
      {
         CHECK_ABSOLUTE_ERROR(Real(k == m ? n : 0), Real(v[k].real()), 4 * n * eps * L);
         CHECK_ABSOLUTE_ERROR(Real(0), Real(v[k].imag()), 4 * n * eps * L);
      }
      // A delta at m transforms to exp(-2*pi*i*m*k/n).
      std::fill(v.begin(), v.end(), Complex(0, 0));
      v[m] = Complex(1, 0);
      f.transform(v);
      for (std::size_t k = 0; k < n; ++k)
      {
         Complex expected = unit_root<Complex>(m * k, n);
         CHECK_ABSOLUTE_ERROR(Real(expected.real()), Real(v[k].real()), 4 * eps * L);
         CHECK_ABSOLUTE_ERROR(Real(-expected.imag()), Real(v[k].imag()), 4 * eps * L);
      }
   }
   // Power-of-two deltas are exact.
   fft<Complex> f(1024);
   std::vector<Complex> v(1024, Complex(0, 0));
   v[0] = Complex(1, 0);
   f.transform(v);
   for (auto const& z : v)
   {
      CHECK_EQUAL(Real(z.real()), Real(1));
      CHECK_EQUAL(Real(z.imag()), Real(0));
   }
}

template <class Complex>
void test_parseval_and_linearity()
{
   using Real = real_t<Complex>;
   using std::sqrt;
   std::mt19937_64 gen(5);
   std::uniform_real_distribution<Real> dis(-1, 1);
   Real eps = std::numeric_limits<Real>::epsilon();
   for (std::size_t n : {3, 16, 100, 127, 256, 1000, 4096, 6000})
   {
      fft<Complex> f(n);
      Real L = Real(log2_length(n));
      std::vector<Complex> x(n), y(n), s(n);
      Real a = dis(gen);
      Real b = dis(gen);
      for (std::size_t j = 0; j < n; ++j)
      {
         x[j] = Complex(dis(gen), dis(gen));
         y[j] = Complex(dis(gen), dis(gen));
         s[j] = a * x[j] + b * y[j];
      }
      Real energy = 0;
      for (auto const& z : x)
      {
         energy += squared_modulus(z);
      }
      f.transform(x);
      f.transform(y);
      f.transform(s);
      // Parseval: sum |X|^2 = n sum |x|^2.
      Real spectral_energy = 0;
      for (auto const& z : x)
      {
         spectral_energy += squared_modulus(z);
      }
      CHECK_ULP_CLOSE(energy, spectral_energy / Real(n), static_cast<std::size_t>(8 * sqrt(Real(n)) * L));
      // Linearity, relative to the size of the coefficients, which is about sqrt(n).
      Real worst = 0;
      for (std::size_t k = 0; k < n; ++k)
      {
         worst = (std::max)(worst, Real(sqrt(squared_modulus(Complex(s[k] - (a * x[k] + b * y[k]))))));
      }
      CHECK_LE(worst, 16 * eps * L * sqrt(Real(n)));
   }
}

template <class Complex>
void test_real_fft()
{
   using Real = real_t<Complex>;
   using std::sqrt;
   std::mt19937_64 gen(11);
   std::uniform_real_distribution<float> dis(-1, 1);
   Real eps = std::numeric_limits<Real>::epsilon();
   for (std::size_t n : {1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 15, 16, 64, 100, 127, 1024, 1025})
   {
      std::vector<Real> x(n);
      std::vector<Complex> z(n);
      for (std::size_t j = 0; j < n; ++j)
      {
         x[j] = dis(gen);
         z[j] = Complex(x[j], 0);
      }
      fft<Complex>(n).transform(z);
      real_fft<Complex> rf(n);
      std::vector<Complex> X(n / 2 + 1);
      rf.transform(x, X);
      Real scale = sqrt(Real(n));
      Real L = Real(log2_length(n));
      for (std::size_t k = 0; k < X.size(); ++k)
      {
         CHECK_ABSOLUTE_ERROR(Real(z[k].real()), Real(X[k].real()), 8 * eps * L * scale);
         CHECK_ABSOLUTE_ERROR(Real(z[k].imag()), Real(X[k].imag()), 8 * eps * L * scale);
      }
   }
}

template <class Complex>
void test_interfaces()
{
   using Real = real_t<Complex>;
   CHECK_THROW(fft<Complex>(0), std::domain_error);
   CHECK_THROW(real_fft<Complex>(0), std::domain_error);
   // n = 1 is the identity.
   fft<Complex> one(1);
   std::vector<Complex> w{Complex(3, 4)};
   one.transform(w);
   CHECK_EQUAL(Real(w[0].real()), Real(3));
   CHECK_EQUAL(Real(w[0].imag()), Real(4));
   one.inverse_transform(w);
   CHECK_EQUAL(Real(w[0].real()), Real(3));
   CHECK_EQUAL(Real(w[0].imag()), Real(4));
}

template <class Complex>
bool identical(std::vector<Complex> const& a, std::vector<Complex> const& b)
{
   if (a.size() != b.size())
   {
      return false;
   }
   for (std::size_t k = 0; k < a.size(); ++k)
   {
      if (!(a[k].real() == b[k].real() && a[k].imag() == b[k].imag()))
      {
         return false;
      }
   }
   return true;
}

// The out-of-place overloads agree with the in-place ones to 0 ulps, leave their input alone, and repeated
// transforms (which reuse the object's workspace) and a second object give the same bits.
// After construction, no transform allocates.
template <class Complex>
void test_overloads(std::vector<std::size_t> const& lengths)
{
   using Real = real_t<Complex>;
   std::mt19937_64 gen(17);
   std::uniform_real_distribution<double> dis(-1, 1);
   for (std::size_t n : lengths)
   {
      std::vector<Complex> x(n);
      for (auto& z : x)
      {
         z = Complex(Real(dis(gen)), Real(dis(gen)));
      }
      std::vector<Complex> const x0 = x;
      fft<Complex> f(n);
      fft<Complex> g(n);
      std::vector<Complex> in_place = x;
      std::vector<Complex> again = x;
      std::vector<Complex> other = x;
      std::vector<Complex> out(n);
      std::vector<Complex> inv_in_place = x;
      std::vector<Complex> inv_out(n);
      std::size_t before = allocations;
      f.transform(in_place);
      f.transform(x, out);
      f.transform(again);
      g.transform(other);
      f.inverse_transform(inv_in_place);
      f.inverse_transform(x, inv_out);
      std::size_t after = allocations;
      if (!CHECK_EQUAL(after - before, std::size_t(0)))
      {
         std::cerr << "  fft of length " << n << " allocated.\n";
      }
      bool ok = CHECK_EQUAL(identical(x, x0), true);
      ok = CHECK_EQUAL(identical(out, in_place), true) && ok;
      ok = CHECK_EQUAL(identical(again, in_place), true) && ok;
      ok = CHECK_EQUAL(identical(other, in_place), true) && ok;
      ok = CHECK_EQUAL(identical(inv_out, inv_in_place), true) && ok;
      if (!ok)
      {
         std::cerr << "  Overloads of the fft of length " << n << " disagree.\n";
      }

      // real_fft: repeated calls reuse its workspace (and its odd-n buffer).
      std::vector<Real> r(n);
      for (auto& t : r)
      {
         t = Real(dis(gen));
      }
      real_fft<Complex> rf(n);
      real_fft<Complex> rg(n);
      std::vector<Complex> R1(n / 2 + 1);
      std::vector<Complex> R2(n / 2 + 1);
      std::vector<Complex> R3(n / 2 + 1);
      before = allocations;
      rf.transform(r, R1);
      rf.transform(r, R2);
      rg.transform(r, R3);
      after = allocations;
      if (!CHECK_EQUAL(after - before, std::size_t(0)))
      {
         std::cerr << "  real_fft of length " << n << " allocated.\n";
      }
      ok = CHECK_EQUAL(identical(R2, R1), true);
      ok = CHECK_EQUAL(identical(R3, R1), true) && ok;
      if (!ok)
      {
         std::cerr << "  Repeated real_fft transforms of length " << n << " disagree.\n";
      }
   }
}

int main()
{
   test_against_naive_dft();
   test_long_lengths();

   test_container_matrix<std::complex<double>, cpp_complex_50>();
   test_container_matrix<cpp_complex_50, cpp_complex_100>();

   test_known_signals<std::complex<float>>();
   test_known_signals<std::complex<double>>();
   test_known_signals<std::complex<long double>>();
   test_known_signals<cpp_complex_50>();

   test_parseval_and_linearity<std::complex<float>>();
   test_parseval_and_linearity<std::complex<double>>();
   test_parseval_and_linearity<std::complex<long double>>();

   test_real_fft<std::complex<float>>();
   test_real_fft<std::complex<double>>();
   test_real_fft<std::complex<long double>>();
   test_real_fft<cpp_complex_50>();

   test_interfaces<std::complex<float>>();
   test_interfaces<std::complex<double>>();
   test_interfaces<cpp_complex_50>();

   // Mixed-radix lengths (including 11-, 13- and 17-smooth ones, and 89, 97 and 97*64 just below the cutoff),
   // and Bluestein ones (101, 103 and 101*64 just above it, 2018, 4099).
   std::vector<std::size_t> const lengths{1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 16, 17, 25, 27, 30, 49, 64, 89,
                                          97, 100, 101, 103, 128, 143, 169, 210, 243, 289, 360, 1000, 1001, 1024, 1156,
                                          2018, 4099, 6208, 6464};
   test_overloads<std::complex<float>>(lengths);
   test_overloads<std::complex<double>>(lengths);
   test_overloads<cpp_complex_50>({1, 2, 3, 4, 5, 6, 7, 8, 11, 12, 13, 16, 17, 30, 64, 89, 97, 100, 101, 143, 289, 360});

   return boost::math::test::report_errors();
}
