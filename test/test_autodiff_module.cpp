//  (C) Copyright Matt Borland 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Module build check for forward-mode automatic differentiation. The fvar class
// template, the make_fvar factory, its arithmetic operators and the
// std::numeric_limits<fvar> specialization all reach a consumer through the
// boost.math module. Exercising make_fvar/derivative and that specialization
// confirms the differentiation component resolves when built against the module.

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/differentiation/autodiff.hpp>
#else
import boost.math;
#endif

#ifndef BOOST_MATH_BUILD_MODULE
#include <limits>
#endif
#include "math_unit_test.hpp"

template <typename T>
void test_spots(T, const char* type_name)
{
   std::cout << "Testing autodiff for type " << type_name << std::endl;

   using boost::math::differentiation::make_fvar;

   // f(x) = x * x at x = 3. value 9, f'(x) = 2x = 6, f''(x) = 2.
   {
      const auto x {make_fvar<T, 2>(static_cast<T>(3))};
      const auto y {x * x};
      CHECK_EQUAL(y.derivative(0), static_cast<T>(9));
      CHECK_EQUAL(y.derivative(1), static_cast<T>(6));
      CHECK_EQUAL(y.derivative(2), static_cast<T>(2));

      // The std::numeric_limits<fvar> specialization must cross the module boundary.
      static_assert(std::numeric_limits<decltype(x)>::is_specialized, "autodiff numeric_limits");
   }

   // Single independent variable f(x) = x at x = 10. value 10, f'(x) = 1, rest 0.
   {
      const auto x {make_fvar<T, 3>(static_cast<T>(10))};
      CHECK_EQUAL(x.derivative(0), static_cast<T>(10));
      CHECK_EQUAL(x.derivative(1), static_cast<T>(1));
      CHECK_EQUAL(x.derivative(2), static_cast<T>(0));
      CHECK_EQUAL(x.derivative(3), static_cast<T>(0));
   }

   // f(x) = 1 / x at x = 16 (a power of two, so every derivative is exact).
   // f'(x) = -1 / x^2, f''(x) = 2 / x^3, f'''(x) = -6 / x^4.
   {
      const T cx {16};
      const auto x {make_fvar<T, 3>(cx)};
      const auto quotient {static_cast<T>(1) / x};
      CHECK_ULP_CLOSE(static_cast<T>(1) / cx, quotient.derivative(0), 2);
      CHECK_ULP_CLOSE(static_cast<T>(-1) / (cx * cx), quotient.derivative(1), 2);
      CHECK_ULP_CLOSE(static_cast<T>(2) / (cx * cx * cx), quotient.derivative(2), 2);
      CHECK_ULP_CLOSE(static_cast<T>(-6) / (cx * cx * cx * cx), quotient.derivative(3), 2);
   }

   // Second independent variable f(x, y) = y at y = 100. df/dy = 1, others 0.
   {
      const auto y {make_fvar<T, 3, 4>(static_cast<T>(100))};
      CHECK_EQUAL(y.derivative(0, 0), static_cast<T>(100));
      CHECK_EQUAL(y.derivative(0, 1), static_cast<T>(1));
      CHECK_EQUAL(y.derivative(1, 0), static_cast<T>(0));
      CHECK_EQUAL(y.derivative(0, 2), static_cast<T>(0));
   }
}

int main()
{
   test_spots(0.0F, "float");
   test_spots(0.0, "double");

   return boost::math::test::report_errors();
}
