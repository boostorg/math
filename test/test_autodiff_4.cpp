//           Copyright Matthew Pulver 2018 - 2019.
// Distributed under the Boost Software License, Version 1.0.
//      (See accompanying file LICENSE_1_0.txt or copy at
//           https://www.boost.org/LICENSE_1_0.txt)

#include "test_autodiff.hpp"

BOOST_AUTO_TEST_SUITE(test_autodiff_4)

BOOST_AUTO_TEST_CASE_TEMPLATE(lround_llround_lltrunc_truncl, T,
                              all_float_types) {
  using boost::math::llround;
  using boost::math::lltrunc;
  using boost::math::lround;
  using boost::multiprecision::llround;
  using boost::multiprecision::lltrunc;
  using boost::multiprecision::lround;
  using detail::llround;
  using detail::lltrunc;
  using detail::lround;
  using detail::truncl;
  using std::truncl;

  constexpr std::size_t m = 3;
  const auto &cx = static_cast<T>(3.25);
  auto x = make_fvar<T, m>(cx);
  auto yl = lround(x);
  BOOST_CHECK_EQUAL(yl, lround(cx));
  auto yll = llround(x);
  BOOST_CHECK_EQUAL(yll, llround(cx));
  BOOST_CHECK_EQUAL(lltrunc(cx), lltrunc(x));

#ifndef BOOST_NO_CXX17_IF_CONSTEXPR
  if constexpr (!bmp::is_number<T>::value &&
                !bmp::is_number_expression<T>::value) {
    BOOST_CHECK_EQUAL(truncl(x), truncl(cx));
  }
#endif
}

BOOST_AUTO_TEST_CASE_TEMPLATE(equality, T, all_float_types) {
  BOOST_MATH_STD_USING
  using boost::math::epsilon_difference;
  using boost::math::fpclassify;
  using boost::math::ulp;
  using std::fpclassify;

  constexpr std::size_t m = 3;
  // check zeros
  {
    auto x = make_fvar<T, m>(0.0);
    auto y = T(-0.0);
    BOOST_CHECK_EQUAL(x.derivative(0u), y);
  }
}

#if defined(BOOST_AUTODIFF_TESTING_INCLUDE_MULTIPRECISION)
BOOST_AUTO_TEST_CASE_TEMPLATE(multiprecision, T, multiprecision_float_types) {
  using boost::multiprecision::fabs;
  using detail::fabs;
  using std::fabs;

  const T eps = 3000 * std::numeric_limits<T>::epsilon();
  constexpr std::size_t Nw = 3;
  constexpr std::size_t Nx = 2;
  constexpr std::size_t Ny = 4;
  constexpr std::size_t Nz = 3;
  const auto w = make_fvar<T, Nw>(11);
  const auto x = make_fvar<T, 0, Nx>(12);
  const auto y = make_fvar<T, 0, 0, Ny>(13);
  const auto z = make_fvar<T, 0, 0, 0, Nz>(14);
  const auto v =
      mixed_partials_f(w, x, y, z); // auto = autodiff_fvar<T,Nw,Nx,Ny,Nz>
  // Calculated from Mathematica symbolic differentiation.
  const T answer = boost::lexical_cast<T>(
      "1976.3196007477977177798818752904187209081211892187"
      "5499076582535951111845769110560421820940516423255314");
  // BOOST_CHECK_CLOSE(v.derivative(Nw,Nx,Ny,Nz), answer, eps); // Doesn't work
  // for cpp_dec_float
  const T relative_error =
      static_cast<T>(fabs(v.derivative(Nw, Nx, Ny, Nz) / answer - 1));
  BOOST_CHECK_LT(relative_error, eps);
}
#endif

BOOST_AUTO_TEST_CASE_TEMPLATE(airy_hpp, T, all_float_types) {
  using boost::multiprecision::min;
  using std::min;

  using test_constants = test_constants_t<T>;
  static constexpr auto m = test_constants::order;

  test_detail::RandomSample<T> x_sampler(-100, 100);
  for (auto i : boost::irange(test_constants::n_samples)) {
    const auto &x = x_sampler.next();
    {
      auto autodiff_v = boost::math::airy_ai(make_fvar<T, m>(x));
      auto anchor_v = boost::math::airy_ai(x);
      BOOST_CHECK_CLOSE(autodiff_v.derivative(0u), anchor_v,
                        1e4 * test_constants::pct_epsilon());
    }

    {
      auto autodiff_v = boost::math::airy_ai_prime(make_fvar<T, m>(x));
      auto anchor_v = boost::math::airy_ai_prime(x);
      BOOST_CHECK_CLOSE(autodiff_v.derivative(0u), anchor_v,
                        1e4 * test_constants::pct_epsilon());
    }

    {
      auto x_ = (min)(x, T(26));
      auto autodiff_v = boost::math::airy_bi(make_fvar<T, m>(x_));
      auto anchor_v = boost::math::airy_bi(x_);
      BOOST_CHECK_CLOSE(autodiff_v.derivative(0u), anchor_v,
                        1e4 * test_constants::pct_epsilon());
    }

    {
      auto x_ = ((min))(x, T(26));
      auto autodiff_v = boost::math::airy_bi_prime(make_fvar<T, m>(x_));
      auto anchor_v = boost::math::airy_bi_prime(x_);
      BOOST_CHECK_CLOSE(autodiff_v.derivative(0u), anchor_v,
                        1e4 * test_constants::pct_epsilon());
    }

    if (i > 0) {
      {
        auto autodiff_v = boost::math::airy_ai_zero<autodiff_fvar<T, m>>(i);
        auto anchor_v = boost::math::airy_ai_zero<T>(i);
        BOOST_CHECK_CLOSE(autodiff_v.derivative(0u), anchor_v,
                          1e4 * test_constants::pct_epsilon());
      }

      {
        auto autodiff_v = boost::math::airy_bi_zero<autodiff_fvar<T, m>>(i);
        auto anchor_v = boost::math::airy_bi_zero<T>(i);
        BOOST_CHECK_CLOSE(autodiff_v.derivative(0u), anchor_v,
                          1e4 * test_constants::pct_epsilon());
      }
    }
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(acosh_hpp, T, all_float_types) {
  using boost::math::acosh;
  using test_constants = test_constants_t<T>;
  static constexpr auto m = test_constants::order;

  test_detail::RandomSample<T> x_sampler{1, 100};
  for (auto i : boost::irange(test_constants::n_samples)) {
    std::ignore = i;
    auto x = x_sampler.next();
    auto autodiff_v = acosh(make_fvar<T, m>(x));
    auto anchor_v = acosh(x);
    BOOST_CHECK_CLOSE(autodiff_v.derivative(0u), anchor_v,
                      1e3 * test_constants::pct_epsilon());
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(asinh_hpp, T, all_float_types) {
  using boost::math::asinh;
  using test_constants = test_constants_t<T>;
  static constexpr auto m = test_constants::order;

  test_detail::RandomSample<T> x_sampler{-100, 100};
  for (auto i : boost::irange(test_constants::n_samples)) {
    std::ignore = i;
    auto x = x_sampler.next();

    auto autodiff_v = asinh(make_fvar<T, m>(x));
    auto anchor_v = asinh(x);
    BOOST_CHECK_CLOSE(autodiff_v.derivative(0u), anchor_v,
                      1e3 * test_constants::pct_epsilon());
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(atanh_hpp, T, all_float_types) {
  using boost::math::nextafter;
  using std::nextafter;

  using boost::math::atanh;
  using test_constants = test_constants_t<T>;
  static constexpr auto m = test_constants::order;

  test_detail::RandomSample<T> x_sampler{-1, 1};
  for (auto i : boost::irange(test_constants::n_samples)) {
    std::ignore = i;
    auto x = nextafter(x_sampler.next(), T(0));

    auto autodiff_v = atanh(make_fvar<T, m>(x));
    auto anchor_v = atanh(x);
    BOOST_CHECK_CLOSE(autodiff_v.derivative(0u), anchor_v,
                      1e3 * test_constants::pct_epsilon());
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(atan_hpp, T, all_float_types) {
  using boost::math::float_prior;
  using boost::math::fpclassify;
  using boost::math::signbit;
  using boost::math::differentiation::detail::atan;
  using boost::multiprecision::atan;
  using boost::multiprecision::fabs;
  using boost::multiprecision::fpclassify;
  using boost::multiprecision::signbit;
  using detail::fabs;
  using std::atan;
  using std::fabs;

  using test_constants = test_constants_t<T>;
  static constexpr auto m = test_constants::order;

  test_detail::RandomSample<T> x_sampler{-1, 1};
  for (auto i : boost::irange(test_constants::n_samples)) {
    std::ignore = i;
    auto x = T(1);
    while (fpclassify(T(fabs(x) - 1)) == FP_ZERO) {
      x = x_sampler.next();
    }

    auto autodiff_v = atan(make_fvar<T, m>(x));
    auto anchor_v = atan(x);
    BOOST_CHECK_CLOSE(autodiff_v.derivative(0u), anchor_v,
                      1e3 * test_constants::pct_epsilon());
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(bernoulli_hpp, T, all_float_types) {

  using boost::multiprecision::min;
  using std::min;
  using test_constants = test_constants_t<T>;
  static constexpr auto m = test_constants::order;

  for (auto i : boost::irange(test_constants::n_samples)) {
    {
      auto autodiff_v = boost::math::bernoulli_b2n<autodiff_fvar<T, m>>(i);
      auto anchor_v = boost::math::bernoulli_b2n<T>(i);
      BOOST_CHECK_CLOSE(autodiff_v.derivative(0u), anchor_v,
                        50 * test_constants::pct_epsilon());
    }
    {
      auto i_ = (min)(19, i);
      auto autodiff_v = boost::math::tangent_t2n<autodiff_fvar<T, m>>(i_);
      auto anchor_v = boost::math::tangent_t2n<T>(i_);
      BOOST_CHECK_CLOSE(autodiff_v.derivative(0u), anchor_v,
                        50 * test_constants::pct_epsilon());
    }
  }
}

// TODO(kbhat): Something in here is very slow with boost::multiprecision
BOOST_AUTO_TEST_CASE_TEMPLATE(bessel_hpp, T, bin_float_types) {
  using boost::math::fpclassify;
  using boost::math::nextafter;
  using boost::math::signbit;
  using boost::math::tools::max;
  using boost::multiprecision::fabs;
  using boost::multiprecision::fpclassify;
  using boost::multiprecision::signbit;
  using detail::fabs;
  using std::fabs;
  using std::max;
  using std::nextafter;

  using test_constants = test_constants_t<T>;
  static constexpr auto m = test_constants::order;
  test_detail::RandomSample<T> v_sampler{-20, 20};
  test_detail::RandomSample<T> x_sampler{
      -boost::math::tools::log_max_value<T>() + 1,
      boost::math::tools::log_max_value<T>() - 1};
  for (auto i : boost::irange(test_constants::n_samples)) {
    auto v = v_sampler.next();
    auto x = x_sampler.next();
    v = ((fpclassify(v) == FP_ZERO) ? 1 : -1) *
        (max)(v, (nextafter)(T(0), ((std::numeric_limits<T>::max))()));
    if (signbit(x)) {
      v = static_cast<T>(boost::math::itrunc(v));
    }

    try {
      auto autodiff_v =
          boost::math::cyl_bessel_i(make_fvar<T, m>(v), make_fvar<T, m>(x));
      auto anchor_v = boost::math::cyl_bessel_i(v, x);
      BOOST_WARN_CLOSE(autodiff_v.derivative(0u), anchor_v,
                       1e4 * test_constants::pct_epsilon());
    } catch (const std::overflow_error &) {
      BOOST_CHECK_THROW(
          boost::math::cyl_bessel_i(make_fvar<T, m>(v), make_fvar<T, m>(x)),
          boost::wrapexcept<std::overflow_error>);
      BOOST_CHECK_THROW(boost::math::cyl_bessel_i(v, x),
                        boost::wrapexcept<std::overflow_error>);
    } catch (...) {
      std::cout << "Inputs: v: " << v << " x: " << x << std::endl;
      std::rethrow_exception(std::current_exception());
    }

    {
      auto x_j = fabs(x) + 1;
      try {
        auto autodiff_v =
            boost::math::cyl_bessel_j(make_fvar<T, m>(v), make_fvar<T, m>(x_j));
        auto anchor_v = boost::math::cyl_bessel_j(v, x_j);
        BOOST_WARN_CLOSE(autodiff_v.derivative(0u), anchor_v,
                         1e4 * test_constants::pct_epsilon());
      } catch (const std::overflow_error &) {
        BOOST_CHECK_THROW(
            boost::math::cyl_bessel_j(make_fvar<T, m>(v), make_fvar<T, m>(x_j)),
            boost::wrapexcept<std::overflow_error>);
        BOOST_CHECK_THROW(boost::math::cyl_bessel_j(v, x_j),
                          boost::wrapexcept<std::overflow_error>);
      } catch (...) {
        std::cout << "Inputs: v: " << v << " x: " << x_j << std::endl;
        std::rethrow_exception(std::current_exception());
      }
    }

    try {
      auto autodiff_v =
          boost::math::cyl_bessel_j_zero(make_fvar<T, m>(v), i + 1);
      auto anchor_v = boost::math::cyl_bessel_j_zero(v, i + 1);
      BOOST_WARN_CLOSE(autodiff_v.derivative(0u), anchor_v,
                       1e4 * test_constants::pct_epsilon());
    } catch (const std::overflow_error &) {
      BOOST_CHECK_THROW(
          boost::math::cyl_bessel_j_zero(make_fvar<T, m>(v), i + 1),
          boost::wrapexcept<std::overflow_error>);
      BOOST_CHECK_THROW(boost::math::cyl_bessel_j_zero(v, i + 1),
                        boost::wrapexcept<std::overflow_error>);
    } catch (...) {
      std::cout << "Inputs: v: " << v << " i: " << (i + 1) << std::endl;
      std::rethrow_exception(std::current_exception());
    }

    {
      auto x_k = fabs(x) + 1;
      try {
        auto autodiff_v =
            boost::math::cyl_bessel_k(make_fvar<T, m>(v), make_fvar<T, m>(x_k));
        auto anchor_v = boost::math::cyl_bessel_k(v, x_k);
        BOOST_WARN_CLOSE(autodiff_v.derivative(0u), anchor_v,
                         1e4 * test_constants::pct_epsilon());
      } catch (const std::overflow_error &) {
        BOOST_CHECK_THROW(
            boost::math::cyl_bessel_k(make_fvar<T, m>(v), make_fvar<T, m>(x_k)),
            boost::wrapexcept<std::overflow_error>);
        BOOST_CHECK_THROW(boost::math::cyl_bessel_k(v, x_k),
                          boost::wrapexcept<std::overflow_error>);
      } catch (...) {
        std::cout << "Inputs: v: " << v << " x: " << x_k << std::endl;
        std::rethrow_exception(std::current_exception());
      }
    }

    {
      auto x_neumann = fabs(x) + 1;
      try {
        auto autodiff_v = boost::math::cyl_neumann(make_fvar<T, m>(v),
                                                   make_fvar<T, m>(x_neumann));
        auto anchor_v = boost::math::cyl_neumann(v, x_neumann);
        BOOST_WARN_CLOSE(autodiff_v.derivative(0u), anchor_v,
                         1e4 * test_constants::pct_epsilon());
      } catch (const std::overflow_error &) {
        BOOST_CHECK_THROW(boost::math::cyl_neumann(make_fvar<T, m>(v),
                                                   make_fvar<T, m>(x_neumann)),
                          boost::wrapexcept<std::overflow_error>);
        BOOST_CHECK_THROW(boost::math::cyl_neumann(v, x_neumann),
                          boost::wrapexcept<std::overflow_error>);
      } catch (...) {
        std::cout << "Inputs: v: " << v << " x: " << x_neumann << std::endl;
        std::rethrow_exception(std::current_exception());
      }
    }

    {
      try {
        auto autodiff_v =
            boost::math::cyl_neumann_zero(make_fvar<T, m>(v), i + 1);
        auto anchor_v = boost::math::cyl_neumann_zero(v, i + 1);
        BOOST_WARN_CLOSE(autodiff_v.derivative(0u), anchor_v,
                         1e4 * test_constants::pct_epsilon());
      } catch (std::overflow_error &) {
        BOOST_CHECK_THROW(
            boost::math::cyl_neumann_zero((make_fvar<T, m>)(v), i + 1),
            boost::wrapexcept<std::overflow_error>);
        BOOST_CHECK_THROW(boost::math::cyl_neumann_zero(v, i + 1),
                          boost::wrapexcept<std::overflow_error>);
      } catch (...) {
        std::cout << "Inputs: v: " << v << " i: " << (i + 1) << std::endl;
        std::rethrow_exception(std::current_exception());
      }
    }

    {
      auto i_ = static_cast<unsigned>(i);
      try {
        auto autodiff_v = boost::math::sph_bessel<autodiff_fvar<T, m>>(
            i_, make_fvar<T, m>(fabs(v)));
        auto anchor_v = boost::math::sph_bessel<T>(i_, fabs(v));
        BOOST_WARN_CLOSE(autodiff_v.derivative(0u), anchor_v,
                         1e4 * test_constants::pct_epsilon());
      } catch (const std::overflow_error &) {
        BOOST_CHECK_THROW(
            (boost::math::sph_bessel<
                autodiff_fvar<T, m>>)(i_, (make_fvar<T, m>)(fabs(v))),
            boost::wrapexcept<std::overflow_error>);
        BOOST_CHECK_THROW(boost::math::sph_bessel<T>(i_, fabs(v)),
                          boost::wrapexcept<std::overflow_error>);
      } catch (...) {
        std::cout << "Inputs: i: " << i_ << " v: " << (fabs(v)) << std::endl;
        std::rethrow_exception(std::current_exception());
      }

      {
        auto v_ =
            (max)(T(fabs(v)),
                  nextafter(T(fabs(v)), 2 * (std::numeric_limits<T>::min)()));
        try {
          auto autodiff_v = boost::math::sph_neumann<autodiff_fvar<T, m>>(
              i_, make_fvar<T, m>(v_));
          auto anchor_v = boost::math::sph_neumann<T>(i_, v_);
          BOOST_WARN_CLOSE(autodiff_v.derivative(0u), anchor_v,
                           1e4 * test_constants::pct_epsilon());
        } catch (const std::overflow_error &) {
          BOOST_CHECK_THROW(((boost::math::sph_neumann<autodiff_fvar<T, m>>(
                                i_, make_fvar<T, m>(v_)))),
                            boost::wrapexcept<std::overflow_error>);
          BOOST_CHECK_THROW(((boost::math::sph_neumann<T>(i_, v_))),
                            boost::wrapexcept<std::overflow_error>);
        } catch (...) {
          std::cout << "Inputs: i: " << i_ << " v: " << v_ << std::endl;
          std::rethrow_exception(std::current_exception());
        }
      }

      try {
        auto autodiff_v = boost::math::cyl_bessel_i_prime(make_fvar<T, m>(v),
                                                          make_fvar<T, m>(x));
        auto anchor_v = boost::math::cyl_bessel_i_prime(v, x);
        BOOST_WARN_CLOSE(autodiff_v.derivative(0u), anchor_v,
                         1e4 * test_constants::pct_epsilon());
      } catch (const std::overflow_error &) {
        BOOST_CHECK_THROW(boost::math::cyl_bessel_i_prime((make_fvar<T, m>)(v),
                                                          (make_fvar<T, m>)(x)),
                          boost::wrapexcept<std::overflow_error>);
        BOOST_CHECK_THROW(boost::math::cyl_bessel_i_prime(v, x),
                          boost::wrapexcept<std::overflow_error>);
      } catch (...) {
        std::cout << "Inputs: v: " << v << " x: " << x << std::endl;
        std::rethrow_exception(std::current_exception());
      }

      {
        auto x_j = fabs(x) + 1;
        try {
          auto autodiff_v = boost::math::cyl_bessel_j_prime(
              make_fvar<T, m>(v), make_fvar<T, m>(x_j));
          auto anchor_v = boost::math::cyl_bessel_j_prime(v, x_j);
          BOOST_WARN_CLOSE(autodiff_v.derivative(0u), anchor_v,
                           1e4 * test_constants::pct_epsilon());
        } catch (const std::overflow_error &) {
          BOOST_CHECK_THROW(boost::math::cyl_bessel_j_prime(
                                (make_fvar<T, m>)(v), (make_fvar<T, m>)(x_j)),
                            boost::wrapexcept<std::overflow_error>);
          BOOST_CHECK_THROW(boost::math::cyl_bessel_j_prime(v, x_j),
                            boost::wrapexcept<std::overflow_error>);
        } catch (...) {
          std::cout << "Inputs: v: " << v << " x_k: " << x_j << std::endl;
          std::rethrow_exception(std::current_exception());
        }
      }

      {
        auto x_k = fabs(x) + 1;
        try {
          auto autodiff_v = boost::math::cyl_bessel_k_prime(
              make_fvar<T, m>(v), make_fvar<T, m>(x_k));
          auto anchor_v = boost::math::cyl_bessel_k_prime(v, x_k);
          BOOST_WARN_CLOSE(autodiff_v.derivative(0u), anchor_v,
                           1e4 * test_constants::pct_epsilon());
        } catch (const std::overflow_error &) {
          BOOST_CHECK_THROW(boost::math::cyl_bessel_k_prime(
                                (make_fvar<T, m>)(v), (make_fvar<T, m>)(x_k)),
                            boost::wrapexcept<std::overflow_error>);
          BOOST_CHECK_THROW(boost::math::cyl_bessel_k_prime(v, x_k),
                            boost::wrapexcept<std::overflow_error>);
        } catch (...) {
          std::cout << "Inputs: v: " << v << " x: " << x_k << std::endl;
          std::rethrow_exception(std::current_exception());
        }
      }

      {
        auto x_neumann = fabs(x) + 1;
        try {
          auto autodiff_v = boost::math::cyl_neumann_prime(
              make_fvar<T, m>(v), make_fvar<T, m>(x_neumann));
          auto anchor_v = boost::math::cyl_neumann_prime(v, x_neumann);
          BOOST_WARN_CLOSE(autodiff_v.derivative(0u), anchor_v,
                           1e4 * test_constants::pct_epsilon());
        } catch (const std::overflow_error &) {
          BOOST_CHECK_THROW(
              boost::math::cyl_neumann_prime((make_fvar<T, m>)(v),
                                             (make_fvar<T, m>)(x_neumann)),
              boost::wrapexcept<std::overflow_error>);
          BOOST_CHECK_THROW(boost::math::cyl_neumann_prime(v, x_neumann),
                            boost::wrapexcept<std::overflow_error>);
        } catch (...) {
          std::cout << "Inputs: v: " << v << " x: " << x_neumann << std::endl;
          std::rethrow_exception(std::current_exception());
        }
      }

      try {
        auto autodiff_v = boost::math::sph_bessel_prime<autodiff_fvar<T, m>>(
            i_, make_fvar<T, m>(fabs(v) + 1));
        auto anchor_v = boost::math::sph_bessel_prime<T>(i_, fabs(v) + 1);
        BOOST_WARN_CLOSE(autodiff_v.derivative(0u), anchor_v,
                         1e4 * test_constants::pct_epsilon());
      } catch (const std::overflow_error &) {
        BOOST_CHECK_THROW(
            (boost::math::sph_neumann<
                autodiff_fvar<T, m>>)(i_, (make_fvar<T, m>)(fabs(v) + 1)),
            boost::wrapexcept<std::overflow_error>);
        BOOST_CHECK_THROW((boost::math::sph_neumann<T>)(i_, fabs(v) + 1),
                          boost::wrapexcept<std::overflow_error>);
      } catch (...) {
        std::cout << "Inputs: i: " << i_ << " v: " << (fabs(v) + 1)
                  << std::endl;
        std::rethrow_exception(std::current_exception());
      }

      try {
        auto autodiff_v = boost::math::sph_neumann_prime<autodiff_fvar<T, m>>(
            i_, make_fvar<T, m>(fabs(v) + 1));
        auto anchor_v = boost::math::sph_neumann_prime<T>(i_, fabs(v) + 1);
        BOOST_WARN_CLOSE(autodiff_v.derivative(0u), anchor_v,
                         1e4 * test_constants::pct_epsilon());
      } catch (const std::overflow_error &) {
        BOOST_CHECK_THROW(
            (boost::math::sph_neumann_prime<
                autodiff_fvar<T, m>>)(i_, (make_fvar<T, m>)(fabs(v) + 1)),
            boost::wrapexcept<std::overflow_error>);
        BOOST_CHECK_THROW((boost::math::sph_neumann_prime<T>)(i_, fabs(v) + 1),
                          boost::wrapexcept<std::overflow_error>);
      } catch (...) {
        std::cout << "Inputs: i: " << i << " v: " << (fabs(v) + 1) << std::endl;
        std::rethrow_exception(std::current_exception());
      }
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
