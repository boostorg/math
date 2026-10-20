//           Copyright Matthew Pulver 2018 - 2019.
// Distributed under the Boost Software License, Version 1.0.
//      (See accompanying file LICENSE_1_0.txt or copy at
//           https://www.boost.org/LICENSE_1_0.txt)

#include "test_autodiff.hpp"
#include <type_traits>
#include <vector>


template <class T>
void constructors() {
  constexpr std::size_t m = 3;
  constexpr std::size_t n = 4;
  // Verify value-initialized instance has all 0 entries.
  const autodiff_fvar<T, m> empty1 = autodiff_fvar<T, m>();
  for (auto i : boost::irange(m + 1)) {
    CHECK_EQUAL(empty1.derivative(i), static_cast<T>(0));
  }
  const auto empty2 = autodiff_fvar<T, m, n>();
  for (auto i : boost::irange(m + 1)) {
    for (auto j : boost::irange(n + 1)) {
      CHECK_EQUAL(empty2.derivative(i, j), static_cast<T>(0));
    }
  }
  // Single variable
  const T cx = 10;
  const auto x = make_fvar<T, m>(cx);
  for (auto i : boost::irange(m + 1)) {
    if (i == 0u) {
      CHECK_EQUAL(x.derivative(i), cx);
    } else if (i == 1) {
      CHECK_EQUAL(x.derivative(i), static_cast<T>(1));
    } else {
      CHECK_EQUAL(x.derivative(i), static_cast<T>(0));
    }
  }
  const autodiff_fvar<T, n> xn = x;
  for (auto i : boost::irange(n + 1)) {
    if (i == 0) {
      CHECK_EQUAL(xn.derivative(i), cx);
    } else if (i == 1) {
      CHECK_EQUAL(xn.derivative(i), static_cast<T>(1));
    } else {
      CHECK_EQUAL(xn.derivative(i), static_cast<T>(0));
    }
  }
  // Second independent variable
  const T cy = 100;
  const auto y = make_fvar<T, m, n>(cy);
  for (auto i : boost::irange(m + 1)) {
    for (auto j : boost::irange(n + 1)) {
      if (i == 0 && j == 0) {
        CHECK_EQUAL(y.derivative(i, j), cy);
      } else if (i == 0 && j == 1) {
        CHECK_EQUAL(y.derivative(i, j), static_cast<T>(1.0));
      } else {
        CHECK_EQUAL(y.derivative(i, j), static_cast<T>(0.0));
      }
    }
  }
}

template <class T>
void implicit_constructors() {
  constexpr std::size_t m = 3;
  const autodiff_fvar<T, m> x = 3;
  const autodiff_fvar<T, m> one = uncast_return(x);
  const autodiff_fvar<T, m> two_and_a_half = 2.5;
  CHECK_EQUAL(static_cast<T>(x), static_cast<T>(3.0));
  CHECK_EQUAL(static_cast<T>(one), static_cast<T>(1.0));
  CHECK_EQUAL(static_cast<T>(two_and_a_half), static_cast<T>(2.5));
  // Only types that T can be constructed from convert implicitly (#424: Eigen expressions must not).
  CHECK_TRUE((std::is_convertible<int, autodiff_fvar<T, m>>::value));
  CHECK_TRUE((std::is_convertible<double, autodiff_fvar<autodiff_fvar<T, 2>, m>>::value));
  CHECK_TRUE((!std::is_convertible<std::vector<T>, autodiff_fvar<T, m>>::value));
  CHECK_TRUE((!std::is_constructible<autodiff_fvar<T, m>, std::vector<T>>::value));
}

template <class T>
void assignment() {
  constexpr std::size_t m = 3;
  constexpr std::size_t n = 4;
  const T cx = 10;
  const T cy = 10;
  autodiff_fvar<T, m, n>
      empty; // Uninitialized variable<> may have non-zero values.
  // Single variable
  auto x = make_fvar<T, m>(cx);
  empty = static_cast<decltype(empty)>(
      x); // Test static_cast of single-variable to double-variable type.
  for (auto i : boost::irange(m + 1)) {
    for (auto j : boost::irange(n + 1)) {
      if (i == 0 && j == 0) {
        CHECK_EQUAL(empty.derivative(i, j), cx);
      } else if (i == 1 && j == 0) {
        CHECK_EQUAL(empty.derivative(i, j), static_cast<T>(1.0));
      } else {
        CHECK_EQUAL(empty.derivative(i, j), static_cast<T>(0.0));
      }
    }
  }
  auto y = make_fvar<T, m, n>(cy);
  empty = y; // default assignment operator
  for (auto i : boost::irange(m + 1)) {
    for (auto j : boost::irange(n + 1)) {
      if (i == 0 && j == 0) {
        CHECK_EQUAL(empty.derivative(i, j), cy);
      } else if (i == 0 && j == 1) {
        CHECK_EQUAL(empty.derivative(i, j), static_cast<T>(1.0));
      } else {
        CHECK_EQUAL(empty.derivative(i, j), static_cast<T>(0.0));
      }
    }
  }
  empty = cx; // set a constant
  for (auto i : boost::irange(m + 1)) {
    for (auto j : boost::irange(n + 1)) {
      if (i == 0 && j == 0) {
        CHECK_EQUAL(empty.derivative(i, j), cx);
      } else {
        CHECK_EQUAL(empty.derivative(i, j), static_cast<T>(0.0));
      }
    }
  }
}

template <class T>
void ostream() {
  constexpr std::size_t m = 3;
  const T cx = 10;
  const auto x = make_fvar<T, m>(cx);
  std::ostringstream ss;
  ss << "x = " << x;
  CHECK_EQUAL(ss.str(), std::string("x = depth(1)(10,1,0,0)"));
  ss.str(std::string());
  const auto scalar = make_fvar<T,0>(cx);
  ss << "scalar = " << scalar;
  CHECK_EQUAL(ss.str(), std::string("scalar = depth(1)(10)"));
}

template <class T>
void addition_assignment() {
  constexpr std::size_t m = 3;
  constexpr std::size_t n = 4;
  const T cx = 10;
  auto sum = autodiff_fvar<T, m, n>(); // zero-initialized
  // Single variable
  const auto x = make_fvar<T, m>(cx);
  sum += x;
  for (auto i : boost::irange(m + 1)) {
    for (auto j : boost::irange(n + 1)) {
      if (i == 0 && j == 0) {
        CHECK_EQUAL(sum.derivative(i, j), cx);
      } else if (i == 1 && j == 0) {
        CHECK_EQUAL(sum.derivative(i, j), T(1));
      } else {
        CHECK_EQUAL(sum.derivative(i, j), T(0));
      }
    }
  }
  // Arithmetic constant
  const T cy = 11;
  sum = 0;
  sum += cy;
  for (auto i : boost::irange(m + 1)) {
    for (auto j : boost::irange(n + 1)) {
      if (i == 0 && j == 0) {
        CHECK_EQUAL(sum.derivative(i, j), cy);
      } else {
        CHECK_EQUAL(sum.derivative(i, j), T(0));
      }
    }
  }
}

template <class T>
void subtraction_assignment() {
  constexpr std::size_t m = 3;
  constexpr std::size_t n = 4;
  const T cx = 10;
  auto sum = autodiff_fvar<T, m, n>(); // zero-initialized
  // Single variable
  const auto x = make_fvar<T, m>(cx);
  sum -= x;
  for (auto i : boost::irange(m + 1)) {
    for (auto j : boost::irange(n + 1)) {
      if (i == 0 && j == 0) {
        CHECK_EQUAL(sum.derivative(i, j), -cx);
      } else if (i == 1 && j == 0) {
        CHECK_EQUAL(sum.derivative(i, j), T(-1));
      } else {
        CHECK_EQUAL(sum.derivative(i, j), T(0));
      }
    }
  }
  // Arithmetic constant
  const T cy = 11;
  sum = 0;
  sum -= cy;
  for (auto i : boost::irange(m + 1)) {
    for (auto j : boost::irange(n + 1)) {
      if (i == 0 && j == 0) {
        CHECK_EQUAL(sum.derivative(i, j), -cy);
      } else {
        CHECK_EQUAL(sum.derivative(i, j), T(0));
      }
    }
  }
}

template <class T>
void multiplication_assignment() {
  // Try explicit bracing based on feedback. Doesn't add very much except 26
  // extra lines.
  constexpr std::size_t m = 3;
  constexpr std::size_t n = 4;
  const T cx = 10;
  auto product = autodiff_fvar<T, m, n>(1); // unit constant
  // Single variable
  auto x = make_fvar<T, m>(cx);
  product *= x;
  for (auto i : boost::irange(m + 1)) {
    for (auto j : boost::irange(n + 1)) {
      if (i == 0 && j == 0) {
        CHECK_EQUAL(product.derivative(i, j), cx);
      } else if (i == 1 && j == 0) {
        CHECK_EQUAL(product.derivative(i, j), T(1));
      } else {
        CHECK_EQUAL(product.derivative(i, j), T(0));
      }
    }
  }
  // Arithmetic constant
  const T cy = 11;
  product = 1;
  product *= cy;
  for (auto i : boost::irange(m + 1)) {
    for (auto j : boost::irange(n + 1)) {
      if (i == 0 && j == 0) {
        CHECK_EQUAL(product.derivative(i, j), cy);
      } else {
        CHECK_EQUAL(product.derivative(i, j), T(0));
      }
    }
  }
  // 0 * inf = nan
  x = make_fvar<T, m>(T(0.0));
  x *= std::numeric_limits<T>::infinity();
  // std::cout << "x = " << x << std::endl;
  for (auto i : boost::irange(m + 1)) {
    if (i == 0) {
      CHECK_TRUE(boost::math::isnan(static_cast<T>(x))); // Correct
      // BOOST_CHECK_EQUAL(x.derivative(i) == 0.0); // Wrong. See
      // multiply_assign_by_root_type().
    } else if (i == 1) {
      CHECK_TRUE(boost::math::isinf(x.derivative(i)));
    } else {
      CHECK_EQUAL(x.derivative(i), static_cast<T>(0.0));
    }
  }
}

template <class T>
void division_assignment() {
  constexpr std::size_t m = 3;
  constexpr std::size_t n = 4;
  const T cx = 16;
  auto quotient = autodiff_fvar<T, m, n>(1); // unit constant
  // Single variable
  const auto x = make_fvar<T, m>(cx);
  quotient /= x;
  CHECK_EQUAL(quotient.derivative(0, 0), 1 / cx);
  CHECK_EQUAL(quotient.derivative(1, 0), static_cast<T>(-1 / pow(cx, 2)));
  CHECK_EQUAL(quotient.derivative(2, 0), static_cast<T>(2 / pow(cx, 3)));
  CHECK_EQUAL(quotient.derivative(3, 0), static_cast<T>(-6 / pow(cx, 4)));
  for (auto i : boost::irange(m + 1)) {
    for (auto j : boost::irange(std::size_t(1), n + 1)) {
      CHECK_EQUAL(quotient.derivative(i, j), T(0));
    }
  }
  // Arithmetic constant
  const T cy = 32;
  quotient = 1;
  quotient /= cy;
  for (auto i : boost::irange(m + 1)) {
    for (auto j : boost::irange(n + 1)) {
      if (i == 0 && j == 0) {
        CHECK_EQUAL(quotient.derivative(i, j), 1 / cy);
      } else {
        CHECK_EQUAL(quotient.derivative(i, j), T(0));
      }
    }
  }
}

template <class T>
void unary_signs() {
  constexpr std::size_t m = 3;
  constexpr std::size_t n = 4;
  const T cx = 16;
  autodiff_fvar<T, m, n> lhs;
  // Single variable
  const auto x = make_fvar<T, m>(cx);
  lhs = static_cast<decltype(lhs)>(-x);
  for (auto i : boost::irange(m + 1)) {
    for (auto j : boost::irange(n + 1)) {
      if (i == 0 && j == 0) {
        CHECK_EQUAL(lhs.derivative(i, j), -cx);
      } else if (i == 1 && j == 0) {
        CHECK_EQUAL(lhs.derivative(i, j), T(-1));
      } else {
        CHECK_EQUAL(lhs.derivative(i, j), T(0));
      }
    }
  }
  lhs = static_cast<decltype(lhs)>(+x);
  for (auto i : boost::irange(m + 1)) {
    for (auto j : boost::irange(n + 1)) {
      if (i == 0 && j == 0) {
        CHECK_EQUAL(lhs.derivative(i, j), cx);
      } else if (i == 1 && j == 0) {
        CHECK_EQUAL(lhs.derivative(i, j), T(1));
      } else {
        CHECK_EQUAL(lhs.derivative(i, j), T(0));
      }
    }
  }
}

// TODO 3 tests for 3 operator+() definitions.
template <class T>
void cast_double() {
  const T ca(13);
  const T i(12);
  constexpr std::size_t m = 3;
  const auto x = make_fvar<T, m>(ca);
  CHECK_TRUE((i) < (x));
  CHECK_TRUE((i * x) == (i * ca));
}

template <class T>
void int_double_casting() {
  const T ca = 3;
  const auto x0 = make_fvar<T, 0>(ca);
  CHECK_EQUAL(static_cast<T>(x0), ca);
  const auto x1 = make_fvar<T, 1>(ca);
  CHECK_EQUAL(static_cast<T>(x1), ca);
  const auto x2 = make_fvar<T, 2>(ca);
  CHECK_EQUAL(static_cast<T>(x2), ca);
}

template <class T>
void scalar_addition() {
  const T ca = 3;
  const T cb = 4;
  const auto sum0 = autodiff_fvar<T, 0>(ca) + autodiff_fvar<T, 0>(cb);
  CHECK_EQUAL(ca + cb, static_cast<T>(sum0));
  const auto sum1 = autodiff_fvar<T, 0>(ca) + cb;
  CHECK_EQUAL(ca + cb, static_cast<T>(sum1));
  const auto sum2 = ca + autodiff_fvar<T, 0>(cb);
  CHECK_EQUAL(ca + cb, static_cast<T>(sum2));
}

template <class T>
void power8() {
  constexpr std::size_t n = 8u;
  const T ca = 3;
  auto x = make_fvar<T, n>(ca);
  // Test operator*=()
  x *= x;
  x *= x;
  x *= x;
  const T power_factorial = boost::math::factorial<T>(n);
  for (auto i : boost::irange(n + 1)) {
    CHECK_RELATIVE_ERROR(
        static_cast<T>(x.derivative(i)),
        static_cast<T>(power_factorial /
                       boost::math::factorial<T>(static_cast<unsigned>(n - i)) *
                       pow(ca, n - i)), (std::numeric_limits<T>::epsilon()) / 100);
  }
  x = make_fvar<T, n>(ca);
  // Test operator*()
  x = x * x * x * x * x * x * x * x;
  for (auto i : boost::irange(n + 1)) {
    CHECK_RELATIVE_ERROR(
        static_cast<T>(x.derivative(i)),
        static_cast<T>(power_factorial /
            boost::math::factorial<T>(static_cast<unsigned>(n - i)) *
            pow(ca, n - i)), (std::numeric_limits<T>::epsilon()) / 100);
  }
}

template <class T>
void dim1_multiplication() {
  constexpr std::size_t m = 2;
  constexpr std::size_t n = 3;
  const T cy = 4;
  auto y0 = make_fvar<T, m>(cy);
  auto y = make_fvar<T, n>(cy);
  y *= y0;
  CHECK_EQUAL(y.derivative(0), cy * cy);
  CHECK_EQUAL(y.derivative(1), 2 * cy);
  CHECK_EQUAL(y.derivative(2), T(2));
  CHECK_EQUAL(y.derivative(3), T(0));
  y = y * cy;
  CHECK_EQUAL(y.derivative(0), cy * cy * cy);
  CHECK_EQUAL(y.derivative(1), 2 * cy * cy);
  CHECK_EQUAL(y.derivative(2), static_cast<T>(2.0 * cy));
  CHECK_EQUAL(y.derivative(3), static_cast<T>(0.0));
}

template <class T>
void dim1and2_multiplication() {
  constexpr std::size_t m = 2;
  constexpr std::size_t n = 3;
  const T cx = 3;
  const T cy = 4;
  auto x = make_fvar<T, m>(cx);
  auto y = make_fvar<T, m, n>(cy);
  y *= x;
  CHECK_EQUAL(y.derivative(0, 0), cx * cy);
  CHECK_EQUAL(y.derivative(0, 1), cx);
  CHECK_EQUAL(y.derivative(1, 0), cy);
  CHECK_EQUAL(y.derivative(1, 1), T(1));
  for (auto i : boost::irange(std::size_t(1), m)) {
    for (auto j : boost::irange(std::size_t(1), n)) {
      if (i == 1 && j == 1) {
        CHECK_EQUAL(y.derivative(i, j), T(1));
      } else {
        CHECK_EQUAL(y.derivative(i, j), T(0));
      }
    }
  }
}

template <class T>
void dim2_addition() {
  constexpr std::size_t m = 2;
  constexpr std::size_t n = 3;
  const T cx = 3;
  const auto x = make_fvar<T, m>(cx);
  CHECK_EQUAL(x.derivative(0), cx);
  CHECK_EQUAL(x.derivative(1), T(1));
  CHECK_EQUAL(x.derivative(2), static_cast<T>(0.0));
  const T cy = 4;
  const auto y = make_fvar<T, m, n>(cy);
  CHECK_EQUAL(static_cast<T>(y.derivative(0)), cy);
  CHECK_EQUAL(static_cast<T>(y.derivative(1)), static_cast<T>(0.0)); // partial of y w.r.t. x.

  CHECK_EQUAL(y.derivative(0, 0), cy);
  CHECK_EQUAL(y.derivative(0, 1), T(1));
  CHECK_EQUAL(y.derivative(1, 0), static_cast<T>(0.0));
  CHECK_EQUAL(y.derivative(1, 1), static_cast<T>(0.0));
  const auto z = x + y;
  CHECK_EQUAL(z.derivative(0, 0), cx + cy);
  CHECK_EQUAL(z.derivative(0, 1), T(1));
  CHECK_EQUAL(z.derivative(1, 0), T(1));
  CHECK_EQUAL(z.derivative(1, 1), static_cast<T>(0.0));
  // The following 4 are unnecessarily more expensive than the previous 4.
  CHECK_EQUAL(z.derivative(0).derivative(0), cx + cy);
  CHECK_EQUAL(z.derivative(0).derivative(1), T(1));
  CHECK_EQUAL(z.derivative(1).derivative(0), T(1));
  CHECK_EQUAL(z.derivative(1).derivative(1), static_cast<T>(0.0));
}

template <class T>
void dim2_multiplication() {
  constexpr std::size_t m = 3;
  constexpr std::size_t n = 4;
  const T cx = 6;
  const auto x = make_fvar<T, m>(cx);
  const T cy = 5;
  const auto y = make_fvar<T, 0, n>(cy);
  const auto z = x * x * y * y * y;
  CHECK_EQUAL(z.derivative(0, 0), cx * cx * cy * cy * cy); // x^2 * y^3
  CHECK_EQUAL(z.derivative(0, 1), cx * cx * 3 * cy * cy);  // x^2 * 3y^2
  CHECK_EQUAL(z.derivative(0, 2), cx * cx * 6 * cy);       // x^2 * 6y
  CHECK_EQUAL(z.derivative(0, 3), cx * cx * 6);            // x^2 * 6
  CHECK_EQUAL(z.derivative(0, 4), T(0));                    // x^2 * 0
  CHECK_EQUAL(z.derivative(1, 0), 2 * cx * cy * cy * cy);  // 2x * y^3
  CHECK_EQUAL(z.derivative(1, 1), 2 * cx * 3 * cy * cy);   // 2x * 3y^2
  CHECK_EQUAL(z.derivative(1, 2), 2 * cx * 6 * cy);        // 2x * 6y
  CHECK_EQUAL(z.derivative(1, 3), 2 * cx * 6);             // 2x * 6
  CHECK_EQUAL(z.derivative(1, 4), T(0));                    // 2x * 0
  CHECK_EQUAL(z.derivative(2, 0), 2 * cy * cy * cy);       // 2 * y^3
  CHECK_EQUAL(z.derivative(2, 1), 2 * 3 * cy * cy);        // 2 * 3y^2
  CHECK_EQUAL(z.derivative(2, 2), 2 * 6 * cy);             // 2 * 6y
  CHECK_EQUAL(z.derivative(2, 3), static_cast<T>(2 * 6));                  // 2 * 6
  CHECK_EQUAL(z.derivative(2, 4), T(0));                    // 2 * 0
  CHECK_EQUAL(z.derivative(3, 0), T(0));                    // 0 * y^3
  CHECK_EQUAL(z.derivative(3, 1), T(0));                    // 0 * 3y^2
  CHECK_EQUAL(z.derivative(3, 2), T(0));                    // 0 * 6y
  CHECK_EQUAL(z.derivative(3, 3), T(0));                    // 0 * 6
  CHECK_EQUAL(z.derivative(3, 4), T(0));                    // 0 * 0
}

template <class T>
void dim2_multiplication_and_subtraction() {
  constexpr std::size_t m = 3;
  constexpr std::size_t n = 4;
  const T cx = 6;
  const auto x = make_fvar<T, m>(cx);
  const T cy = 5;
  const auto y = make_fvar<T, 0, n>(cy);
  const auto z = x * x - y * y;
  CHECK_EQUAL(z.derivative(0, 0), cx * cx - cy * cy);
  CHECK_EQUAL(z.derivative(0, 1), -2 * cy);
  CHECK_EQUAL(z.derivative(0, 2), T(-2));
  CHECK_EQUAL(z.derivative(0, 3), T(0));
  CHECK_EQUAL(z.derivative(0, 4), T(0));
  CHECK_EQUAL(z.derivative(1, 0), 2 * cx);
  CHECK_EQUAL(z.derivative(2, 0), T(2));
  for (auto i : boost::irange(std::size_t(1), m + 1)) {
    for (auto j : boost::irange(std::size_t(1), n + 1)) {
      CHECK_EQUAL(z.derivative(i, j), T(0));
    }
  }
}

template <class T>
void inverse() {
  constexpr std::size_t m = 3;
  const T cx = 4;
  const auto x = make_fvar<T, m>(cx);
  const auto xinv = x.inverse();
  CHECK_EQUAL(xinv.derivative(0), 1 / cx);
  CHECK_EQUAL(xinv.derivative(1), static_cast<T>(-1 / pow(cx, 2)));
  CHECK_EQUAL(xinv.derivative(2), static_cast<T>(2 / pow(cx, 3)));
  CHECK_EQUAL(xinv.derivative(3), static_cast<T>(-6 / pow(cx, 4)));
  const auto zero = make_fvar<T, m>(0);
  const auto inf = zero.inverse();
  for (auto i : boost::irange(m + 1)) {
    CHECK_EQUAL(inf.derivative(i),
                      (i % 2 == 1 ? -1 : 1) *
                          std::numeric_limits<T>::infinity());
  }
}

template <class T>
void division() {
  constexpr std::size_t m = 3;
  constexpr std::size_t n = 4;
  const T cx = 16;
  auto x = make_fvar<T, m>(cx);
  const T cy = 4;
  auto y = make_fvar<T, 1, n>(cy);
  auto z = x * x / (y * y);
  CHECK_EQUAL(z.derivative(0, 0), cx * cx / (cy * cy)); // x^2 * y^-2
  CHECK_EQUAL(z.derivative(0, 1), static_cast<T>(cx * cx * (-2) * pow(cy, -3)));
  CHECK_EQUAL(z.derivative(0, 2), static_cast<T>(cx * cx * (6) * pow(cy, -4)));
  CHECK_EQUAL(z.derivative(0, 3), static_cast<T>(cx * cx * (-24) * pow(cy, -5)));
  CHECK_EQUAL(z.derivative(0, 4), static_cast<T>(cx * cx * (120) * pow(cy, -6)));
  CHECK_EQUAL(z.derivative(1, 0), 2 * cx / (cy * cy));
  CHECK_EQUAL(z.derivative(1, 1), static_cast<T>(2 * cx * (-2) * pow(cy, -3)));
  CHECK_EQUAL(z.derivative(1, 2), static_cast<T>(2 * cx * (6) * pow(cy, -4)));
  CHECK_EQUAL(z.derivative(1, 3), static_cast<T>(2 * cx * (-24) * pow(cy, -5)));
  CHECK_EQUAL(z.derivative(1, 4), static_cast<T>(2 * cx * (120) * pow(cy, -6)));
  CHECK_EQUAL(z.derivative(2, 0), 2 / (cy * cy));
  CHECK_EQUAL(z.derivative(2, 1), static_cast<T>(2 * (-2) * pow(cy, -3)));
  CHECK_EQUAL(z.derivative(2, 2), static_cast<T>(2 * (6) * pow(cy, -4)));
  CHECK_EQUAL(z.derivative(2, 3), static_cast<T>(2 * (-24) * pow(cy, -5)));
  CHECK_EQUAL(z.derivative(2, 4), static_cast<T>(2 * (120) * pow(cy, -6)));
  for (auto j : boost::irange(n + 1)) {
    CHECK_EQUAL(z.derivative(3, j), static_cast<T>(0.0));
  }

  auto x1 = make_fvar<T, m>(cx);
  auto z1 = x1 / cy;
  CHECK_EQUAL(z1.derivative(0), cx / cy);
  CHECK_EQUAL(z1.derivative(1), 1 / cy);
  CHECK_EQUAL(z1.derivative(2), static_cast<T>(0.0));
  CHECK_EQUAL(z1.derivative(3), static_cast<T>(0.0));
  auto y2 = make_fvar<T, m, n>(cy);
  auto z2 = cx / y2;
  CHECK_EQUAL(z2.derivative(0, 0), cx / cy);
  CHECK_EQUAL(z2.derivative(0, 1), static_cast<T>(-cx / pow(cy, 2)));
  CHECK_EQUAL(z2.derivative(0, 2), static_cast<T>(2 * cx / pow(cy, 3)));
  CHECK_EQUAL(z2.derivative(0, 3), static_cast<T>(-6 * cx / pow(cy, 4)));
  CHECK_EQUAL(z2.derivative(0, 4), static_cast<T>(24 * cx / pow(cy, 5)));
  for (auto i : boost::irange(std::size_t(1), m + 1)) {
    for (auto j : boost::irange(n + 1)) {
      CHECK_EQUAL(z2.derivative(i, j), static_cast<T>(0.0));
    }
  }

  const auto z3 = y / x;
  CHECK_EQUAL(z3.derivative(0, 0), cy / cx);
  CHECK_EQUAL(z3.derivative(0, 1), 1 / cx);
  CHECK_EQUAL(z3.derivative(1, 0), static_cast<T>(-cy / pow(cx, 2)));
  CHECK_EQUAL(z3.derivative(1, 1), static_cast<T>(-1 / pow(cx, 2)));
  CHECK_EQUAL(z3.derivative(2, 0), static_cast<T>(2 * cy / pow(cx, 3)));
  CHECK_EQUAL(z3.derivative(2, 1), static_cast<T>(2 / pow(cx, 3)));
  CHECK_EQUAL(z3.derivative(3, 0), static_cast<T>(-6 * cy / pow(cx, 4)));
  CHECK_EQUAL(z3.derivative(3, 1), static_cast<T>(-6 / pow(cx, 4)));
  for (auto i : boost::irange(m + 1)) {
    for (auto j : boost::irange(std::size_t(2), n + 1)) {
      CHECK_EQUAL(z3.derivative(i, j), static_cast<T>(0.0));
    }
  }
}

template <class T>
void equality() {
  constexpr std::size_t m = 3;
  constexpr std::size_t n = 4;
  const T cx = 10;
  const T cy = 10;
  const auto x = make_fvar<T, m>(cx);
  const auto y = make_fvar<T, 0, n>(cy);
  CHECK_TRUE((x) == (y));
  CHECK_TRUE((x) == (cy));
  CHECK_TRUE((cx) == (y));
  CHECK_TRUE((cy) == (x));
  CHECK_TRUE((y) == (cx));
}

template <class T>
void inequality() {
  constexpr std::size_t m = 3;
  constexpr std::size_t n = 4;
  const T cx = 10;
  const T cy = 11;
  const auto x = make_fvar<T, m>(cx);
  const auto y = make_fvar<T, 0, n>(cy);
  CHECK_TRUE((x) != (y));
  CHECK_TRUE((x) != (cy));
  CHECK_TRUE((cx) != (y));
  CHECK_TRUE((cy) != (x));
  CHECK_TRUE((y) != (cx));
}

template <class T>
void less_than_or_equal_to() {
  constexpr std::size_t m = 3;
  constexpr std::size_t n = 4;
  const T cx = 10;
  const T cy = 11;
  const auto x = make_fvar<T, m>(cx);
  const auto y = make_fvar<T, 0, n>(cy);
  CHECK_TRUE((x) <= (y));
  CHECK_TRUE((x) <= (y - 1));
  CHECK_TRUE((x) < (y));
  CHECK_TRUE((x) <= (cy));
  CHECK_TRUE((x) <= (cy - 1));
  CHECK_TRUE((x) < (cy));
  CHECK_TRUE((cx) <= (y));
  CHECK_TRUE((cx) <= (y - 1));
  CHECK_TRUE((cx) < (y));
}

template <class T>
void greater_than_or_equal_to() {
  constexpr std::size_t m = 3;
  constexpr std::size_t n = 4;
  const T cx = 11;
  const T cy = 10;
  const auto x = make_fvar<T, m>(cx);
  const auto y = make_fvar<T, 0, n>(cy);
  CHECK_TRUE((x) >= (y));
  CHECK_TRUE((x) >= (y + 1));
  CHECK_TRUE((x) > (y));
  CHECK_TRUE((x) >= (cy));
  CHECK_TRUE((x) >= (cy + 1));
  CHECK_TRUE((x) > (cy));
  CHECK_TRUE((cx) >= (y));
  CHECK_TRUE((cx) >= (y + 1));
  CHECK_TRUE((cx) > (y));
}

template <class T>
void fabs_test() {
  using bmp::fabs;
  using detail::fabs;
  using std::fabs;
  constexpr std::size_t m = 3;
  const T cx = 11;
  const auto x = make_fvar<T, m>(cx);
  auto a = fabs(x);
  CHECK_EQUAL(a.derivative(0), fabs(cx));
  CHECK_EQUAL(a.derivative(1), T(1));
  CHECK_EQUAL(a.derivative(2), T(0));
  CHECK_EQUAL(a.derivative(3), T(0));
  a = fabs(-x);
  CHECK_EQUAL(a.derivative(0), fabs(cx));
  CHECK_EQUAL(a.derivative(1), T(1)); // fabs(-x) = fabs(x)
  CHECK_EQUAL(a.derivative(2), T(0));
  CHECK_EQUAL(a.derivative(3), T(0));
  const auto xneg = make_fvar<T, m>(-cx);
  a = fabs(xneg);
  CHECK_EQUAL(a.derivative(0), fabs(cx));
  CHECK_EQUAL(a.derivative(1), -T(1));
  CHECK_EQUAL(a.derivative(2), T(0));
  CHECK_EQUAL(a.derivative(3), T(0));
  const auto zero = make_fvar<T, m>(0);
  a = fabs(zero);
  for (auto i : boost::irange(m + 1)) {
    CHECK_EQUAL(a.derivative(i), T(0));
  }
}

template <class T>
void ceil_and_floor() {
  using bmp::ceil;
  using bmp::floor;
  using std::ceil;
  using std::floor;
  constexpr std::size_t m = 3;
  T tests[]{T(-1.5), T(0.0), T(1.5)};
  for (T &test : tests) {
    const auto x = make_fvar<T, m>(test);
    auto c = ceil(x);
    auto f = floor(x);
    CHECK_EQUAL(c.derivative(0), ceil(test));
    CHECK_EQUAL(f.derivative(0), floor(test));
    for (auto i : boost::irange(std::size_t(1), m + 1)) {
      CHECK_EQUAL(c.derivative(i), T(0));
      CHECK_EQUAL(f.derivative(i), T(0));
    }
  }
}

int main()
{
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { constructors<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { implicit_constructors<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { assignment<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { ostream<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { addition_assignment<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { subtraction_assignment<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { multiplication_assignment<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { division_assignment<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { unary_signs<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { cast_double<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { int_double_casting<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { scalar_addition<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { power8<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { dim1_multiplication<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { dim1and2_multiplication<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { dim2_addition<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { dim2_multiplication<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { dim2_multiplication_and_subtraction<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { inverse<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { division<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { equality<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { inequality<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { less_than_or_equal_to<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { greater_than_or_equal_to<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { fabs_test<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { ceil_and_floor<typename decltype(t)::type>(); });
   return boost::math::test::report_errors();
}
