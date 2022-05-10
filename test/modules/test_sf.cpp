// Copyright John Maddock 2022.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <iomanip>
#include <limits>
#include <tuple>
#include <complex>
#include <array>

import boost.math.special_functions;

int main()
{
   std::cout << boost::math::tgamma(2.5) << std::endl;
   std::cout << boost::math::lgamma(2.5) << std::endl;
   std::cout << boost::math::tgamma1pm1(2.5) << std::endl;
   std::cout << boost::math::digamma(2.5) << std::endl;
   std::cout << boost::math::trigamma(2.5) << std::endl;
   std::cout << boost::math::tgamma_ratio(2.5, 2.75) << std::endl;
   std::cout << boost::math::tgamma_delta_ratio(2.5, -0.25) << std::endl;
   std::cout << boost::math::gamma_p(2.5, 3.25) << std::endl;
   std::cout << boost::math::gamma_q(2.5, 3.25) << std::endl;
   std::cout << boost::math::tgamma(2.5, 3.25) << std::endl;
   std::cout << boost::math::tgamma_lower(2.5, 3.25) << std::endl;
   std::cout << boost::math::gamma_p_inv(2.5, 0.5) << std::endl;
   std::cout << boost::math::gamma_q_inv(2.5, 0.5) << std::endl;
   std::cout << boost::math::gamma_q_inva(2.5, 0.5) << std::endl;
   std::cout << boost::math::gamma_p_inva(2.5, 0.5) << std::endl;
   std::cout << boost::math::gamma_p_derivative(2.5, 2.75) << std::endl;
   std::cout << boost::math::factorial<double>(7) << std::endl;
   std::cout << boost::math::unchecked_factorial<double>(7) << std::endl;
   std::cout << boost::math::double_factorial<double>(7) << std::endl;
   std::cout << boost::math::rising_factorial(12.5, 7) << std::endl;
   std::cout << boost::math::falling_factorial(12.5, 7) << std::endl;
   std::cout << boost::math::erf(1.5) << std::endl;
   std::cout << boost::math::erfc(1.5) << std::endl;
   std::cout << boost::math::erf_inv(0.5) << std::endl;
   std::cout << boost::math::erfc_inv(0.5) << std::endl;
   std::cout << boost::math::zeta(0.5) << std::endl;
   std::cout << boost::math::sin_pi(0.5) << std::endl;
   std::cout << boost::math::cos_pi(0.5) << std::endl;
   std::cout << boost::math::log1p(0.5) << std::endl;
   std::cout << boost::math::expm1(0.5) << std::endl;
}
