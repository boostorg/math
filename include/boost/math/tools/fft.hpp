//  Copyright Nick Thompson 2026.
//  Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at
//  https://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_FFT_HPP
#define BOOST_MATH_TOOLS_FFT_HPP

#include <boost/math/tools/config.hpp>
#include <boost/math/tools/detail/complex_from_real.hpp>
#include <boost/math/tools/assert.hpp>
#include <boost/math/special_functions/sin_pi.hpp>
#include <boost/math/special_functions/cos_pi.hpp>

#ifndef BOOST_MATH_BUILD_MODULE
#include <algorithm>
#include <cstddef>
#include <initializer_list>
#include <iterator>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>
#endif

BOOST_MATH_NAMESPACE_BEGIN namespace tools {

namespace fft_detail {

// The real type of a complex type: std::complex<T> gives T, and a Boost.Multiprecision complex
// number gives its component type. Both have a real() member returning it.
template <class Complex>
using real_type_t = BOOST_MATH_NAMESPACE::tools::detail::real_from_complex_t<Complex>;

template <class Complex>
inline Complex conj(Complex const& z)
{
   return Complex(z.real(), -z.imag());
}

// Explicit products: std::complex operator* may call a C99 Annex G helper (e.g. __muldc3)
// that is much slower than the four multiplies we need.
template <class Complex>
inline Complex mul(Complex const& a, Complex const& b)
{
   return Complex(a.real() * b.real() - a.imag() * b.imag(),
                  a.real() * b.imag() + a.imag() * b.real());
}

// conj(a)*conj(b)
template <class Complex>
inline Complex mul_conj(Complex const& a, Complex const& b)
{
   return Complex(a.real() * b.real() - a.imag() * b.imag(),
                  -(a.real() * b.imag() + a.imag() * b.real()));
}

// exp(i*pi*p/q) for integers p >= 0, q > 0.
// The reduction to the first octant is done in integer arithmetic, so the only rounding
// in the argument is that of a single quotient in [0, 1/4]; each twiddle is then within
// a couple of ulps of the exact root of unity, with no accumulation from a recurrence.
template <class Complex>
Complex cis_pi(std::size_t p, std::size_t q)
{
   using Real = real_type_t<Complex>;
   using BOOST_MATH_NAMESPACE::sin_pi;
   using BOOST_MATH_NAMESPACE::cos_pi;
   p %= 2 * q;
   bool negate_cos = false;
   bool negate_sin = false;
   if (p >= q)
   {
      // theta -> theta - pi
      p -= q;
      negate_cos = true;
      negate_sin = true;
   }
   if (2 * p > q)
   {
      // theta -> pi - theta
      p = q - p;
      negate_cos = !negate_cos;
   }
   // Now theta = pi*p/q lies in [0, pi/2].
   Real c;
   Real s;
   if (4 * p > q)
   {
      // theta = pi/2 - pi*x with x in [0, 1/4).
      Real x = Real(q - 2 * p) / Real(2 * q);
      c = sin_pi(x);
      s = cos_pi(x);
   }
   else
   {
      Real x = Real(p) / Real(q);
      c = cos_pi(x);
      s = sin_pi(x);
   }
   if (negate_cos)
   {
      c = -c;
   }
   if (negate_sin)
   {
      s = -s;
   }
   return Complex(c, s);
}

// The largest prime factor the mixed-radix FFT takes; a length with a larger one uses Bluestein.
// A radix-p stage costs O(p) per point, Bluestein a constant factor of about 3 to 4 over a 5-smooth FFT;
// measured (complex<double>, n = p 2^k and 1000 p), Bluestein starts to win between p = 97 and p = 127,
// and the mixed-radix FFT is the more accurate up to about the same point.
constexpr std::size_t max_prime_radix = 97;

// The prime factors of n, in increasing order, with multiplicity.
inline std::vector<std::size_t> prime_factors(std::size_t n)
{
   std::vector<std::size_t> f;
   for (std::size_t p = 2; p * p <= n; ++p)
   {
      while (n % p == 0)
      {
         f.push_back(p);
         n /= p;
      }
   }
   if (n > 1)
   {
      f.push_back(n);
   }
   return f;
}

// Whether every prime factor of n is at most max_prime_radix; other lengths use Bluestein.
inline bool is_smooth(std::size_t n)
{
   std::vector<std::size_t> f = prime_factors(n);
   return f.empty() || f.back() <= max_prime_radix;
}

inline std::size_t next_5_smooth(std::size_t n)
{
   for (std::size_t m = n;; ++m)
   {
      std::size_t r = m;
      for (std::size_t p : {2, 3, 5})
      {
         while (r % p == 0)
         {
            r /= p;
         }
      }
      if (r == 1)
      {
         return m;
      }
   }
}

// The length of Bluestein's convolution: the smallest 2^a 3^b 5^c >= n, unless the next power of 2
// is within 20% of it, since radix-4 stages are cheaper than radix-3 and radix-5 ones.
inline std::size_t convolution_length(std::size_t n)
{
   std::size_t smooth = next_5_smooth(n);
   std::size_t power_of_two = 1;
   while (power_of_two < n)
   {
      power_of_two <<= 1;
   }
   return 5 * power_of_two <= 6 * smooth ? power_of_two : smooth;
}

// Radix-2 and radix-4 butterflies, and a butterfly for any odd prime radix R that pairs t with R - t:
//    y[u] = a[0] + sum_{t=1}^{(R-1)/2} (a[t] + a[R-t]) cos(2 pi t u/R) - i (a[t] - a[R-t]) sin(2 pi t u/R).
template <class Complex, std::size_t R>
struct butterfly
{
   using Real = real_type_t<Complex>;
   static void apply(Complex* a, Real const* c, Real const* s)
   {
      constexpr std::size_t h = (R - 1) / 2;
      Complex sum[h + 1];
      Complex diff[h + 1];
      Complex y0 = a[0];
      for (std::size_t t = 1; t <= h; ++t)
      {
         sum[t] = a[t] + a[R - t];
         diff[t] = a[t] - a[R - t];
         y0 += sum[t];
      }
      for (std::size_t u = 1; u <= h; ++u)
      {
         Real even_re = a[0].real();
         Real even_im = a[0].imag();
         Real odd_re = 0;
         Real odd_im = 0;
         for (std::size_t t = 1; t <= h; ++t)
         {
            std::size_t k = (t * u) % R;
            even_re += sum[t].real() * c[k];
            even_im += sum[t].imag() * c[k];
            odd_re += diff[t].real() * s[k];
            odd_im += diff[t].imag() * s[k];
         }
         // even -/+ i*odd
         a[u] = Complex(even_re + odd_im, even_im - odd_re);
         a[R - u] = Complex(even_re - odd_im, even_im + odd_re);
      }
      a[0] = y0;
   }
};

template <class Complex>
struct butterfly<Complex, 2>
{
   using Real = real_type_t<Complex>;
   static void apply(Complex* a, Real const*, Real const*)
   {
      Complex t = a[0];
      a[0] = t + a[1];
      a[1] = t - a[1];
   }
};

template <class Complex>
struct butterfly<Complex, 4>
{
   using Real = real_type_t<Complex>;
   static void apply(Complex* a, Real const*, Real const*)
   {
      Complex s02 = a[0] + a[2];
      Complex d02 = a[0] - a[2];
      Complex s13 = a[1] + a[3];
      Complex d13 = a[1] - a[3];
      // -i*d13
      Complex r(d13.imag(), -d13.real());
      a[0] = s02 + s13;
      a[1] = d02 + r;
      a[2] = s02 - s13;
      a[3] = d02 - r;
   }
};

// One Stockham autosort stage (decimation in frequency). The current sub-transforms have length L = R*m,
// and there are s of them, interleaved with stride s. Each writes R sub-transforms of length m:
//    y[q + s*(R*p + u)] = w_L^{p*u} * sum_t x[q + s*(p + t*m)] w_R^{t*u}.
template <class Complex, std::size_t R>
void stockham_stage(Complex const* x, Complex* y, std::size_t m, std::size_t s, Complex const* twiddles,
                    real_type_t<Complex> const* c, real_type_t<Complex> const* sn)
{
   Complex a[R];
   if (s == 1)
   {
      // The first stage: one sub-transform per p, so make p the inner (vectorizable) loop.
      // The p = 0 twiddles are exactly 1.
      for (std::size_t p = 0; p < m; ++p)
      {
         for (std::size_t t = 0; t < R; ++t)
         {
            a[t] = x[p + t * m];
         }
         butterfly<Complex, R>::apply(a, c, sn);
         Complex const* w = twiddles + p * (R - 1);
         y[R * p] = a[0];
         for (std::size_t u = 1; u < R; ++u)
         {
            y[R * p + u] = mul(a[u], w[u - 1]);
         }
      }
      return;
   }
   for (std::size_t q = 0; q < s; ++q)
   {
      for (std::size_t t = 0; t < R; ++t)
      {
         a[t] = x[q + s * t * m];
      }
      butterfly<Complex, R>::apply(a, c, sn);
      for (std::size_t u = 0; u < R; ++u)
      {
         y[q + s * u] = a[u];
      }
   }
   for (std::size_t p = 1; p < m; ++p)
   {
      Complex const* w = twiddles + p * (R - 1);
      for (std::size_t q = 0; q < s; ++q)
      {
         for (std::size_t t = 0; t < R; ++t)
         {
            a[t] = x[q + s * (p + t * m)];
         }
         butterfly<Complex, R>::apply(a, c, sn);
         y[q + s * R * p] = a[0];
         for (std::size_t u = 1; u < R; ++u)
         {
            y[q + s * (R * p + u)] = mul(a[u], w[u - 1]);
         }
      }
   }
}

// The odd-prime butterfly of butterfly<Complex, R> for a radix R known only at run time.
// sum and diff hold (R + 1)/2 values each. c[k] = cos(2 pi k/R), s[k] = sin(2 pi k/R).
template <class Complex>
void generic_butterfly(Complex* a, std::size_t R, real_type_t<Complex> const* c, real_type_t<Complex> const* sn,
                       Complex* sum, Complex* diff)
{
   using Real = real_type_t<Complex>;
   std::size_t h = (R - 1) / 2;
   Complex y0 = a[0];
   for (std::size_t t = 1; t <= h; ++t)
   {
      sum[t] = a[t] + a[R - t];
      diff[t] = a[t] - a[R - t];
      y0 += sum[t];
   }
   for (std::size_t u = 1; u <= h; ++u)
   {
      Real even_re = a[0].real();
      Real even_im = a[0].imag();
      Real odd_re = 0;
      Real odd_im = 0;
      // k = t*u mod R, updated without a division.
      std::size_t k = 0;
      for (std::size_t t = 1; t <= h; ++t)
      {
         k += u;
         if (k >= R)
         {
            k -= R;
         }
         even_re += sum[t].real() * c[k];
         even_im += sum[t].imag() * c[k];
         odd_re += diff[t].real() * sn[k];
         odd_im += diff[t].imag() * sn[k];
      }
      a[u] = Complex(even_re + odd_im, even_im - odd_re);
      a[R - u] = Complex(even_re - odd_im, even_im + odd_re);
   }
   a[0] = y0;
}

// stockham_stage for an odd prime radix R known only at run time; work holds 2R + 2 values.
template <class Complex>
void generic_stage(Complex const* x, Complex* y, std::size_t R, std::size_t m, std::size_t s, Complex const* twiddles,
                   real_type_t<Complex> const* c, real_type_t<Complex> const* sn, Complex* work)
{
   Complex* a = work;
   Complex* sum = work + R;
   Complex* diff = sum + (R + 1) / 2;
   for (std::size_t p = 0; p < m; ++p)
   {
      Complex const* w = twiddles + p * (R - 1);
      for (std::size_t q = 0; q < s; ++q)
      {
         for (std::size_t t = 0; t < R; ++t)
         {
            a[t] = x[q + s * (p + t * m)];
         }
         generic_butterfly(a, R, c, sn, sum, diff);
         y[q + s * R * p] = a[0];
         for (std::size_t u = 1; u < R; ++u)
         {
            y[q + s * (R * p + u)] = p == 0 ? a[u] : mul(a[u], w[u - 1]);
         }
      }
   }
}

// A forward FFT of a length whose prime factors are all at most max_prime_radix, done by Stockham
// stages that ping-pong between the data and a scratch buffer, so no digit reversal is needed.
template <class Complex>
class smooth_fft
{
public:
   using Real = real_type_t<Complex>;

   smooth_fft() = default;

   explicit smooth_fft(std::size_t n) : m_n(n)
   {
      // Every twiddle is a power of w = exp(-2*pi*i/n); compute the powers we need directly,
      // using exp(-2*pi*i*(n-k)/n) = conj(exp(-2*pi*i*k/n)), which is exact.
      std::vector<Complex> roots(n);
      for (std::size_t k = 0; 2 * k <= n; ++k)
      {
         roots[k] = conj(cis_pi<Complex>(2 * k, n));
         if (k != 0 && 2 * k != n)
         {
            roots[n - k] = conj(roots[k]);
         }
      }
      // Radices: 4s, then a 2, then the odd primes in increasing order.
      std::vector<std::size_t> f = prime_factors(n);
      std::size_t twos = static_cast<std::size_t>(std::count(f.begin(), f.end(), std::size_t(2)));
      for (std::size_t i = 0; i < twos / 2; ++i)
      {
         m_radices.push_back(4);
      }
      if (twos % 2 != 0)
      {
         m_radices.push_back(2);
      }
      for (std::size_t p : f)
      {
         BOOST_MATH_ASSERT(p <= max_prime_radix);
         if (p != 2)
         {
            m_radices.push_back(p);
         }
      }
      // Stage i has s = product of the earlier radices, L = n/s and m = L/R; its twiddle for (p, u) is
      // w_L^{p*u} = roots[p*u*s], with p*u*s < n.
      std::size_t s = 1;
      for (std::size_t R : m_radices)
      {
         std::size_t m = n / (s * R);
         m_offsets.push_back(m_twiddles.size());
         for (std::size_t p = 0; p < m; ++p)
         {
            for (std::size_t u = 1; u < R; ++u)
            {
               m_twiddles.push_back(roots[p * u * s]);
            }
         }
         s *= R;
      }
      // cos(2*pi*k/R) and sin(2*pi*k/R) for each odd radix R, one table per distinct radix.
      for (std::size_t i = 0; i < m_radices.size(); ++i)
      {
         std::size_t R = m_radices[i];
         std::size_t offset = 0;
         if (R % 2 != 0)
         {
            auto earlier = std::find(m_radices.begin(), m_radices.begin() + static_cast<std::ptrdiff_t>(i), R);
            if (earlier != m_radices.begin() + static_cast<std::ptrdiff_t>(i))
            {
               offset = m_trig_offsets[static_cast<std::size_t>(earlier - m_radices.begin())];
            }
            else
            {
               offset = m_cos.size();
               for (std::size_t k = 0; k < R; ++k)
               {
                  Complex z = cis_pi<Complex>(2 * k, R);
                  m_cos.push_back(z.real());
                  m_sin.push_back(z.imag());
               }
            }
            if (R > 7)
            {
               m_generic_work = (std::max)(m_generic_work, 2 * R + 2);
            }
         }
         m_trig_offsets.push_back(offset);
      }
   }

   std::size_t size() const noexcept { return m_n; }

   // The scratch space execute needs: m_n values (none for n = 1), and room for a generic butterfly.
   std::size_t scratch_size() const noexcept { return m_radices.empty() ? 0 : m_n + m_generic_work; }

   // Transforms the m_n contiguous values at data in place, using scratch_size() values at scratch,
   // which must not overlap data.
   void execute(Complex* data, Complex* scratch) const
   {
      if (m_radices.empty())
      {
         return;
      }
      Complex* x = data;
      Complex* y = scratch;
      std::size_t s = 1;
      for (std::size_t i = 0; i < m_radices.size(); ++i)
      {
         std::size_t R = m_radices[i];
         stage(i, x, y, R, m_n / (s * R), s, scratch + m_n);
         std::swap(x, y);
         s *= R;
      }
      if (x != data)
      {
         std::copy(x, x + m_n, data);
      }
   }

   // Transforms in into out, leaving in unchanged; in, out and scratch must not overlap.
   // The stages alternate between out and scratch so that the last one writes out.
   void execute(Complex const* in, Complex* out, Complex* scratch) const
   {
      std::size_t stages = m_radices.size();
      if (stages == 0)
      {
         std::copy(in, in + m_n, out);
         return;
      }
      Complex const* x = in;
      std::size_t s = 1;
      for (std::size_t i = 0; i < stages; ++i)
      {
         std::size_t R = m_radices[i];
         Complex* y = (stages - 1 - i) % 2 == 0 ? out : scratch;
         stage(i, x, y, R, m_n / (s * R), s, scratch + m_n);
         x = y;
         s *= R;
      }
   }

private:
   // Stage i, of radix R; work is room for a generic butterfly.
   void stage(std::size_t i, Complex const* x, Complex* y, std::size_t R, std::size_t m, std::size_t s, Complex* work) const
   {
      Complex const* w = m_twiddles.data() + m_offsets[i];
      Real const* c = m_cos.data() + m_trig_offsets[i];
      Real const* sn = m_sin.data() + m_trig_offsets[i];
      switch (R)
      {
      case 2:
         stockham_stage<Complex, 2>(x, y, m, s, w, c, sn);
         break;
      case 3:
         stockham_stage<Complex, 3>(x, y, m, s, w, c, sn);
         break;
      case 4:
         stockham_stage<Complex, 4>(x, y, m, s, w, c, sn);
         break;
      case 5:
         stockham_stage<Complex, 5>(x, y, m, s, w, c, sn);
         break;
      case 7:
         stockham_stage<Complex, 7>(x, y, m, s, w, c, sn);
         break;
      // Compile-time butterflies for 11 and 13 are 5 to 15% faster than the generic one.
      case 11:
         stockham_stage<Complex, 11>(x, y, m, s, w, c, sn);
         break;
      case 13:
         stockham_stage<Complex, 13>(x, y, m, s, w, c, sn);
         break;
      default:
         generic_stage(x, y, R, m, s, w, c, sn, work);
         break;
      }
   }

   std::size_t m_n = 0;
   std::vector<std::size_t> m_radices;
   std::vector<std::size_t> m_offsets;
   std::vector<Complex> m_twiddles;
   std::vector<Real> m_cos;
   std::vector<Real> m_sin;
   // Where stage i's cos/sin table starts in m_cos and m_sin.
   std::vector<std::size_t> m_trig_offsets;
   // Workspace for a generic (radix > 7) butterfly, after the ping-pong buffer.
   std::size_t m_generic_work = 0;
};

// The documented requirement is a size() member, not iterators:
template <class Container>
inline std::size_t container_size(Container const& v)
{
   return static_cast<std::size_t>(v.size());
}

} // namespace fft_detail

BOOST_MATH_EXPORT template <class Complex>
class real_fft;

// Discrete Fourier transform of a fixed length n, templated on the complex type
// (std::complex<float>, std::complex<double>, std::complex<long double>, or a Boost.Multiprecision
// complex type such as cpp_complex_50):
//
//    X[k] = sum_{j=0}^{n-1} x[j] exp(-2*pi*i*j*k/n),
//
// with the inverse normalised by 1/n, so that inverse_transform(transform(x)) == x.
// Lengths whose prime factors are all at most 97 use a mixed-radix Stockham FFT (radix 4 where possible).
// Every other length uses Bluestein's algorithm, with a convolution whose length is 2^a 3^b 5^c.
// The constructor builds all the tables and allocates all the workspace, so transforms never allocate.
// Like an FFTW plan with its buffers, an object is used by one thread at a time: give each thread its own.
BOOST_MATH_EXPORT template <class Complex>
class fft
{
public:
   using value_type = Complex;
   using real_type = fft_detail::real_type_t<Complex>;

   explicit fft(std::size_t n) : m_n(n)
   {
      if (n == 0)
      {
         throw std::domain_error("The length of the transform must be positive.");
      }
      if (fft_detail::is_smooth(n))
      {
         m_plan = fft_detail::smooth_fft<Complex>(n);
         m_work.resize(m_plan.scratch_size());
      }
      else
      {
         m_plan = fft_detail::smooth_fft<Complex>(fft_detail::convolution_length(2 * n - 1));
         // The convolution buffer, then the scratch space of its FFTs.
         m_work.resize(m_plan.size() + m_plan.scratch_size());
         build_chirp();
      }
   }

   std::size_t size() const noexcept { return m_n; }

   // Transforms v in place. v is a contiguous container of n values of type Complex,
   // such as std::vector<Complex> or std::array<Complex, n>.
   template <class ContiguousContainer>
   void transform(ContiguousContainer& v)
   {
      Complex* x = data_of(v);
      forward(x, x);
   }

   template <class ContiguousContainer>
   void inverse_transform(ContiguousContainer& v)
   {
      Complex* x = data_of(v);
      inverse(x, x);
   }

   // Out of place: out receives the transform of in, which is left unchanged.
   // Both hold n values of type Complex contiguously, and they must not overlap.
   template <class InputContainer, class OutputContainer>
   void transform(InputContainer const& in, OutputContainer& out)
   {
      forward(data_of(in), data_of(out));
   }

   template <class InputContainer, class OutputContainer>
   void inverse_transform(InputContainer const& in, OutputContainer& out)
   {
      inverse(data_of(in), data_of(out));
   }

private:
   template <class> friend class real_fft;

   template <class ContiguousContainer>
   auto data_of(ContiguousContainer& v) const -> decltype(v.data())
   {
      static_assert(std::is_same<typename std::decay<decltype(*v.data())>::type, Complex>::value,
                    "The container's elements must be of the transform's complex type.");
      if (fft_detail::container_size(v) != m_n)
      {
         throw std::domain_error("The length of the data does not match the length of the transform.");
      }
      return v.data();
   }

   // Bluestein: with c[j] = exp(-i*pi*j^2/n), jk = (j^2 + k^2 - (k-j)^2)/2 gives
   //    X[k] = c[k] * sum_j (x[j] c[j]) conj(c[k-j]),
   // a linear convolution that is computed as a cyclic one of length m >= 2n - 1.
   void build_chirp()
   {
      std::size_t len = m_plan.size();
      m_chirp.resize(m_n);
      // j^2 mod 2n, updated exactly so that j^2 never overflows.
      std::size_t two_n = 2 * m_n;
      std::size_t r = 0;
      for (std::size_t j = 0; j < m_n; ++j)
      {
         m_chirp[j] = fft_detail::conj(fft_detail::cis_pi<Complex>(r, m_n));
         // (j+1)^2 = j^2 + 2j + 1
         r = (r + (2 * j + 1) % two_n) % two_n;
      }
      std::vector<Complex> b(len, Complex(0, 0));
      for (std::size_t j = 0; j < m_n; ++j)
      {
         b[j] = fft_detail::conj(m_chirp[j]);
         if (j != 0)
         {
            b[len - j] = b[j];
         }
      }
      m_plan.execute(b.data(), m_work.data() + len);
      // Store B/m: the convolution's inverse FFT is done as conj(F(conj(A)*conj(B/m))).
      real_type m = real_type(len);
      for (std::size_t k = 0; k < len; ++k)
      {
         b[k] = Complex(b[k].real() / m, b[k].imag() / m);
      }
      m_kernel = std::move(b);
   }

   // out = F(in). in and out are either the same array or don't overlap.
   void forward(Complex const* in, Complex* out)
   {
      if (m_chirp.empty())
      {
         if (in == out)
         {
            m_plan.execute(out, m_work.data());
         }
         else
         {
            m_plan.execute(in, out, m_work.data());
         }
         return;
      }
      std::size_t m = m_plan.size();
      Complex* a = m_work.data();
      Complex* scratch = a + m;
      for (std::size_t j = 0; j < m_n; ++j)
      {
         a[j] = fft_detail::mul(in[j], m_chirp[j]);
      }
      std::fill(a + m_n, a + m, Complex(0, 0));
      m_plan.execute(a, scratch);
      for (std::size_t k = 0; k < m; ++k)
      {
         a[k] = fft_detail::mul_conj(a[k], m_kernel[k]);
      }
      m_plan.execute(a, scratch);
      for (std::size_t k = 0; k < m_n; ++k)
      {
         out[k] = fft_detail::mul(fft_detail::conj(a[k]), m_chirp[k]);
      }
   }

   // out = conj(F(conj(in)))/n; in and out are either the same array or don't overlap.
   void inverse(Complex const* in, Complex* out)
   {
      for (std::size_t j = 0; j < m_n; ++j)
      {
         out[j] = fft_detail::conj(in[j]);
      }
      forward(out, out);
      real_type n = real_type(m_n);
      for (std::size_t j = 0; j < m_n; ++j)
      {
         out[j] = Complex(out[j].real() / n, -out[j].imag() / n);
      }
   }

   std::size_t m_n;
   // The FFT that does the work: of length n itself, or of the Bluestein convolution length.
   fft_detail::smooth_fft<Complex> m_plan;
   std::vector<Complex> m_chirp;
   std::vector<Complex> m_kernel;
   // Workspace: the Stockham ping-pong buffer, or Bluestein's convolution buffer followed by that.
   std::vector<Complex> m_work;
};

// Transform of n real samples, templated on the complex type of the result. The input holds n values
// of the real type of Complex; the output must hold exactly n/2 + 1 values of type Complex, and receives
// X[0], ..., X[n/2]. The rest follow from X[n-k] = conj(X[k]). Even n packs the input into a complex FFT
// of length n/2; odd n uses a complex FFT of length n.
BOOST_MATH_EXPORT template <class Complex>
class real_fft
{
public:
   using value_type = Complex;
   using real_type = fft_detail::real_type_t<Complex>;

   explicit real_fft(std::size_t n) : m_n(n), m_fft(n % 2 == 0 && n > 0 ? n / 2 : n)
   {
      if (n % 2 != 0)
      {
         m_z.resize(n);
      }
      else
      {
         // exp(-2*pi*i*k/n) for k in [0, n/2].
         m_twiddles.resize(n / 2 + 1);
         for (std::size_t k = 0; k <= n / 2; ++k)
         {
            m_twiddles[k] = fft_detail::conj(fft_detail::cis_pi<Complex>(2 * k, n));
         }
      }
   }

   std::size_t size() const noexcept { return m_n; }

   // input: n real values; output: exactly n/2 + 1 complex values, which are overwritten.
   // Any container with operator[] works for the input; the output must be contiguous (data() and size()).
   // As with fft, the constructor allocates all the workspace, and an object is used by one thread at a time.
   template <class RealContainer, class ComplexContainer>
   void transform(RealContainer const& input, ComplexContainer& output)
   {
      static_assert(std::is_same<typename std::decay<decltype(*output.data())>::type, Complex>::value,
                    "The output container's elements must be of the transform's complex type.");
      if (fft_detail::container_size(input) != m_n)
      {
         throw std::domain_error("The length of the input does not match the length of the transform.");
      }
      if (fft_detail::container_size(output) != m_n / 2 + 1)
      {
         throw std::domain_error("The output must hold n/2 + 1 values.");
      }
      Complex* out = output.data();
      if (m_n % 2 != 0)
      {
         Complex* z = m_z.data();
         for (std::size_t j = 0; j < m_n; ++j)
         {
            z[j] = Complex(real_type(input[j]), real_type(0));
         }
         m_fft.forward(z, z);
         std::copy(z, z + (m_n / 2 + 1), out);
         return;
      }
      // Pack the input into the first n/2 output slots and transform it there.
      std::size_t h = m_n / 2;
      for (std::size_t j = 0; j < h; ++j)
      {
         out[j] = Complex(real_type(input[2 * j]), real_type(input[2 * j + 1]));
      }
      m_fft.forward(out, out);
      // With Z = DFT(z): E[k] = (Z[k] + conj(Z[h-k]))/2 is the transform of the even samples,
      // O[k] = (Z[k] - conj(Z[h-k]))/(2i) that of the odd ones, and X[k] = E[k] + exp(-2*pi*i*k/n) O[k].
      // X[k] and X[h-k] use the same two values of Z, so compute them together, in place.
      Complex z0 = out[0];
      out[0] = combine(z0, z0, 0);
      out[h] = combine(z0, z0, h);
      for (std::size_t k = 1; 2 * k < h; ++k)
      {
         Complex zk = out[k];
         Complex zh = out[h - k];
         out[k] = combine(zk, zh, k);
         out[h - k] = combine(zh, zk, h - k);
      }
      if (h % 2 == 0)
      {
         out[h / 2] = combine(out[h / 2], out[h / 2], h / 2);
      }
   }

private:
   // X[k] from zk = Z[k] and zh = Z[h-k].
   Complex combine(Complex const& zk, Complex const& zh, std::size_t k) const
   {
      Complex zc = fft_detail::conj(zh);
      Complex e((zk.real() + zc.real()) / 2, (zk.imag() + zc.imag()) / 2);
      Complex d((zk.real() - zc.real()) / 2, (zk.imag() - zc.imag()) / 2);
      Complex o(d.imag(), -d.real());
      Complex t = fft_detail::mul(m_twiddles[k], o);
      return Complex(e.real() + t.real(), e.imag() + t.imag());
   }

   std::size_t m_n;
   fft<Complex> m_fft;
   std::vector<Complex> m_twiddles;
   // Odd n: the input as complex values, transformed in place.
   std::vector<Complex> m_z;
};

} BOOST_MATH_NAMESPACE_END

#endif
