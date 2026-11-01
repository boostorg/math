//  Copyright Nick Thompson, 2026
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Temporary diagnostic: on one CI machine, long double hypergeometric_1F1 is wrong by 1e-7 to 1e-4 for
// large negative a, where A&S 13.3.7 (Tricomi) is used. This computes each step of that method in long
// double and in a cpp_bin_float with the same 64-bit significand, which uses no libm, and prints the
// machine, the x87 and SSE control words, and every step when a final value differs.

#include <cfloat>
#include <iostream>

#if LDBL_MANT_DIG == 64

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/math/special_functions/hypergeometric_1F1.hpp>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#ifdef __GLIBC__
#include <gnu/libc-version.h>
#endif

using E = boost::multiprecision::number<boost::multiprecision::cpp_bin_float<64, boost::multiprecision::digit_base_2, void, std::int16_t, -16382, 16383>, boost::multiprecision::et_off>;
using pol_t = boost::math::policies::policy<>;

template <class T>
using steps_t = std::vector<std::pair<std::string, T>>;

template <class T, class F>
void step(steps_t<T>& s, const std::string& name, F f)
{
   T v;
#ifndef BOOST_NO_EXCEPTIONS
   try
#endif
   {
      v = f();
   }
#ifndef BOOST_NO_EXCEPTIONS
   catch (const std::exception& e)
   {
      std::cout << "  " << name << " threw: " << e.what() << "\n";
      v = std::numeric_limits<T>::quiet_NaN();
   }
#endif
   s.emplace_back(name, v);
}

// The steps of A&S 13.3.7 as hypergeometric_1F1_AS_13_3_7_tricomi takes them; the first one is the result:
template <class T>
steps_t<T> steps(const char* as, const char* bs, const char* zs)
{
   using std::exp; using std::log; using std::pow; using std::sqrt; using std::sin; using std::cos; using std::fabs;
   const T a = T(E(as)), b = T(E(bs)), z = T(E(zs));
   const T arg = (b / 2 - a) * z;
   const T x = 2 * sqrt(arg);
   steps_t<T> s;
   step(s, "hypergeometric_1F1(a, b, z)", [&] { return boost::math::hypergeometric_1F1(a, b, z); });
   step(s, "tricomi(a, b, z) * exp(log_scale)", [&] { long long ls = 0; T r = boost::math::detail::hypergeometric_1F1_AS_13_3_7_tricomi(a, b, z, pol_t(), ls); return T(r * exp(T(ls))); });
   step(s, "divergent_fallback(a, b, z) * exp(log_scale)", [&] { long long ls = 0; T r = boost::math::detail::hypergeometric_1F1_divergent_fallback(a, b, z, pol_t(), ls); return T(r * exp(T(ls))); });
   step(s, "tgamma(b)", [&] { return boost::math::tgamma(b); });
   step(s, "lgamma(b)", [&] { return boost::math::lgamma(b); });
   step(s, "exp(z / 2)", [&] { return T(exp(z / 2)); });
   step(s, "bessel_arg = (b / 2 - a) z", [&] { return arg; });
   step(s, "x = 2 sqrt(bessel_arg)", [&] { return x; });
   step(s, "log(bessel_arg)", [&] { return T(log(arg)); });
   step(s, "sin(x)", [&] { return T(sin(x)); });
   step(s, "cos(x)", [&] { return T(cos(x)); });
   step(s, "pow(bessel_arg, (b - 1) / 2)", [&] { return T(pow(arg, (b - 1) / 2)); });
   for (int k = 0; k < 4; ++k)
      step(s, "cyl_bessel_j(b - 2 + " + std::to_string(k) + ", x)", [&] { return boost::math::cyl_bessel_j(b - 2 + k, x); });
   step(s, "cyl_bessel_j(b + 30, x)", [&] { return boost::math::cyl_bessel_j(b + 30, x); });
   {
      boost::math::detail::hypergeometric_1F1_AS_13_3_7_tricomi_series<T, pol_t> series(a, b, z, pol_t());
      step(s, "series.scale()", [&] { return T(series.scale()); });
      for (int k = 0; k < 12; ++k)
         step(s, "series term " + std::to_string(k), [&] { return series(); });
   }
   {
      boost::math::detail::hypergeometric_1F1_AS_13_3_7_tricomi_series<T, pol_t> series(a, b, z, pol_t());
      std::uintmax_t max_iter = boost::math::policies::get_max_series_iterations<pol_t>();
      T norm = 0, sum = 0;
      step(s, "checked_sum_series", [&] { sum = boost::math::tools::checked_sum_series(series, boost::math::policies::get_epsilon<T, pol_t>(), max_iter, T(0), norm); return sum; });
      step(s, "norm / |sum| (retries above 1 / root_epsilon)", [&] { return T(norm / fabs(sum)); });
      step(s, "terms summed", [&] { return T(max_iter); });
   }
   return s;
}

void print_machine()
{
#ifdef __GLIBC__
   std::cout << "glibc " << gnu_get_libc_version() << "\n";
#endif
   std::ifstream ver("/proc/version");
   std::string line;
   if (std::getline(ver, line))
      std::cout << line << "\n";
   std::ifstream cpu("/proc/cpuinfo");
   while (std::getline(cpu, line) && !line.empty())
   {
      if (line.rfind("vendor_id", 0) == 0 || line.rfind("model", 0) == 0 || line.rfind("cpu family", 0) == 0 || line.rfind("stepping", 0) == 0 || line.rfind("microcode", 0) == 0 || line.rfind("flags", 0) == 0 || line.rfind("bugs", 0) == 0)
         std::cout << line << "\n";
   }
#if defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__))
   unsigned short cw = 0;
   unsigned int mxcsr = 0;
   __asm__ __volatile__("fnstcw %0" : "=m"(cw));
   __asm__ __volatile__("stmxcsr %0" : "=m"(mxcsr));
   std::cout << std::hex << "x87 control word 0x" << cw << " (precision control " << ((cw >> 8) & 3) << ", rounding " << ((cw >> 10) & 3) << "), MXCSR 0x" << mxcsr << std::dec << "\n";
#endif
}

int main()
{
   const char* rows[][3] = {
      { "-814723.75", "-1.45330195664428174496e-05", "1.21108460007235407829e-05" },
      { "-13547.703125", "-1.45330195664428174496e-05", "15.87335205078125" },
      { "-835.008544921875", "19.0480194091796875", "15.87335205078125" },
   };
   std::ostringstream report;
   bool differs = false;
   for (const auto& r : rows)
   {
      steps_t<long double> n = steps<long double>(r[0], r[1], r[2]);
      steps_t<E> e = steps<E>(r[0], r[1], r[2]);
      report << "\nhypergeometric_1F1(" << r[0] << ", " << r[1] << ", " << r[2] << "), long double vs 64-bit cpp_bin_float:\n";
      for (std::size_t i = 0; i < n.size(); ++i)
      {
         long double ref = static_cast<long double>(e[i].second);
         long double rel = (n[i].second == ref) ? 0 : std::fabs((n[i].second - ref) / ref);
         report << std::setprecision(21) << "  " << n[i].first << ":\n    " << n[i].second << "\n    " << ref << std::setprecision(3) << "   relative difference " << rel << (rel > 1e-15L ? "   <--" : "") << "\n";
      }
      long double ref = static_cast<long double>(e[0].second);
      if (!(std::fabs((n[0].second - ref) / ref) <= 1e-10L))
         differs = true;
   }
   if (differs)
   {
      std::cout << "long double hypergeometric_1F1 differs from a 64-bit cpp_bin_float on this machine:\n";
      print_machine();
      std::cout << report.str();
      return 1;
   }
   std::cout << "long double and 64-bit cpp_bin_float agree\n";
   return 0;
}

#else

int main()
{
   std::cout << "long double is not the x87 80-bit type here, nothing to diagnose\n";
   return 0;
}

#endif
