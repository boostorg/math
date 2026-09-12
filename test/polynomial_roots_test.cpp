// Copyright Nick Thompson 2026.
// Distributed under the Boost Software License, Version 1.0.

#include <boost/math/tools/polynomial_roots.hpp>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

namespace {

double backward_error(std::vector<double> const& p, std::complex<double> z)
{
   std::complex<double> value = p.back();
   double scale = std::abs(p.back());
   for(std::size_t i = p.size() - 1; i-- > 0;)
   {
      value = value * z + p[i];
      scale = scale * std::abs(z) + std::abs(p[i]);
   }
   return std::abs(value) / (std::max)(scale, (std::numeric_limits<double>::min)());
}

void check_residuals(std::vector<double> const& p,
                     std::vector<std::complex<double> > const& roots,
                     double tolerance = 2e-8)
{
   assert(roots.size() + 1 == p.size());
   for(auto root : roots)
   {
      assert(std::isfinite(root.real()));
      assert(std::isfinite(root.imag()));
      assert(backward_error(p, root) <= tolerance);
   }
}

std::vector<double> polynomial_from_real_roots(std::vector<double> const& roots)
{
   std::vector<double> p(1, 1);
   for(double root : roots)
   {
      std::vector<double> next(p.size() + 1);
      for(std::size_t i = 0; i < p.size(); ++i)
      {
         next[i] -= root * p[i];
         next[i + 1] += p[i];
      }
      p.swap(next);
   }
   return p;
}

} // namespace

int main()
{
   using boost::math::tools::polynomial_roots;

   auto quadratic = polynomial_roots(std::vector<double>{1, 0, 1});
   check_residuals({1, 0, 1}, quadratic);
   assert(quadratic.size() == 2);

   std::vector<double> with_zero_roots{0, 0, -2, 1};
   auto zero_roots = polynomial_roots(with_zero_roots);
   check_residuals(with_zero_roots, zero_roots);
   assert(std::count(zero_roots.begin(), zero_roots.end(), std::complex<double>(0, 0)) == 2);

   std::mt19937_64 generator(0x5eed);
   std::uniform_real_distribution<double> distribution(-2, 2);
   for(unsigned degree = 2; degree <= 32; ++degree)
   {
      std::vector<double> known(degree);
      for(auto& root : known)
         root = distribution(generator);
      auto p = polynomial_from_real_roots(known);
      check_residuals(p, polynomial_roots(p));
   }

   // A sparse polynomial exercises the finite radial bounds needed when the
   // Newton polygon has a single edge.
   std::vector<double> roots_of_unity(257);
   roots_of_unity.front() = -1;
   roots_of_unity.back() = 1;
   check_residuals(roots_of_unity, polynomial_roots(roots_of_unity), 2e-8);

   bool threw = false;
   try { polynomial_roots(std::vector<double>{0, 0}); }
   catch(std::domain_error const&) { threw = true; }
   assert(threw);

   threw = false;
   try { polynomial_roots(std::vector<double>{1, std::numeric_limits<double>::infinity()}); }
   catch(std::domain_error const&) { threw = true; }
   assert(threw);
}
