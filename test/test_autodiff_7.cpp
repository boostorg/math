//           Copyright Matthew Pulver 2018 - 2019.
// Distributed under the Boost Software License, Version 1.0.
//      (See accompanying file LICENSE_1_0.txt or copy at
//           https://www.boost.org/LICENSE_1_0.txt)

#include "test_autodiff.hpp"

BOOST_AUTO_TEST_SUITE(test_autodiff_7)

BOOST_AUTO_TEST_CASE_TEMPLATE(expm1_hpp, T, all_float_types) {
  using boost::math::differentiation::detail::log;
  using boost::multiprecision::log;
  using std::log;
  using test_constants = test_constants_t<T>;
  static constexpr auto m = test_constants::order;
  test_detail::RandomSample<T> x_sampler{-log(T(2000)), log(T(2000))};
  for (auto i : boost::irange(test_constants::n_samples)) {
    std::ignore = i;
    auto x = x_sampler.next();
    BOOST_CHECK_CLOSE(boost::math::expm1(make_fvar<T, m>(x)).derivative(0u),
                      boost::math::expm1(x),
                      50 * test_constants::pct_epsilon());
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(factorials_hpp, T, all_float_types) {
  using test_constants = test_constants_t<T>;
  static constexpr auto m = test_constants::order;
  test_detail::RandomSample<T> x_sampler{0, 28};
  for (auto i :
       boost::irange(static_cast<unsigned>(test_constants::n_samples))) {
    {
      auto fact_i = boost::math::factorial<T>(i);
      auto autodiff_v = boost::math::factorial<autodiff_fvar<T, m>>(i);
      BOOST_CHECK_CLOSE(autodiff_v.derivative(0u), fact_i,
                        50 * test_constants::pct_epsilon());
    }

    {
      auto fact_i = boost::math::unchecked_factorial<T>(i);
      auto autodiff_v =
          boost::math::unchecked_factorial<autodiff_fvar<T, m>>(i);
      BOOST_CHECK_CLOSE(autodiff_v.derivative(0u), fact_i,
                        50 * test_constants::pct_epsilon());
    }

    {
      auto fact_i = boost::math::unchecked_factorial<T>(i);
      auto autodiff_v =
          boost::math::unchecked_factorial<autodiff_fvar<T, m>>(i);
      BOOST_CHECK_CLOSE(autodiff_v.derivative(0u), fact_i,
                        50 * test_constants::pct_epsilon());
    }

    {
      auto fact_i = boost::math::double_factorial<T>(i);
      auto autodiff_v = boost::math::double_factorial<autodiff_fvar<T, m>>(i);
      BOOST_CHECK_CLOSE(autodiff_v.derivative(0u), fact_i,
                        50 * test_constants::pct_epsilon());
    }

    auto x = x_sampler.next();
    {
      auto fact_i = boost::math::rising_factorial<T>(x, static_cast<int>(i));
      auto autodiff_v = make_fvar<T, m>(fact_i);
      BOOST_CHECK_CLOSE(autodiff_v.derivative(0u), fact_i,
                        50 * test_constants::pct_epsilon());
    }

    {
      auto fact_i =
          boost::math::falling_factorial<T>(x, test_constants::n_samples - i);
      auto autodiff_v = make_fvar<T, m>(fact_i);
      BOOST_CHECK_CLOSE(autodiff_v.derivative(0u), fact_i,
                        50 * test_constants::pct_epsilon());
    }
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(fpclassify_hpp, T, all_float_types) {
  using boost::math::fpclassify;
  using boost::math::isfinite;
  using boost::math::isinf;
  using boost::math::isnan;
  using boost::math::isnormal;
  using boost::multiprecision::fpclassify;
  using boost::multiprecision::isfinite;
  using boost::multiprecision::isinf;
  using boost::multiprecision::isnan;
  using boost::multiprecision::isnormal;

  using test_constants = test_constants_t<T>;
  static constexpr auto m = test_constants::order;
  test_detail::RandomSample<T> x_sampler{-1000, 1000};
  for (auto i : boost::irange(test_constants::n_samples)) {
    std::ignore = i;

    BOOST_CHECK_EQUAL(fpclassify(make_fvar<T, m>(0)), FP_ZERO);
    BOOST_CHECK_EQUAL(fpclassify(make_fvar<T, m>(10)), FP_NORMAL);
    BOOST_CHECK_EQUAL(
        fpclassify(make_fvar<T, m>(std::numeric_limits<T>::infinity())),
        FP_INFINITE);
    BOOST_CHECK_EQUAL(
        fpclassify(make_fvar<T, m>(std::numeric_limits<T>::quiet_NaN())),
        FP_NAN);
    if (std::numeric_limits<T>::has_denorm != std::denorm_absent) {
      BOOST_CHECK_EQUAL(
          fpclassify(make_fvar<T, m>(std::numeric_limits<T>::denorm_min())),
          FP_SUBNORMAL);
    }

    BOOST_CHECK(isfinite(make_fvar<T, m>(0)));
    BOOST_CHECK(isnormal(make_fvar<T, m>((std::numeric_limits<T>::min)())));
    BOOST_CHECK(
        !isnormal(make_fvar<T, m>(std::numeric_limits<T>::denorm_min())));
    BOOST_CHECK(isinf(make_fvar<T, m>(std::numeric_limits<T>::infinity())));
    BOOST_CHECK(isnan(make_fvar<T, m>(std::numeric_limits<T>::quiet_NaN())));
  }
}

// multiprecision types are breaking in here due to a lexical_cast error
BOOST_AUTO_TEST_CASE_TEMPLATE(gamma_hpp, T, bin_float_types) {
  using boost::math::nextafter;
  using boost::multiprecision::nextafter;

  using test_constants = test_constants_t<T>;
  static constexpr auto m = test_constants::order;
  test_detail::RandomSample<T> z_sampler{0, 34};
  test_detail::RandomSample<T> a_sampler{0, 34};
  for (auto i : boost::irange(test_constants::n_samples)) {
    std::ignore = i;
    auto z = z_sampler.next();
    {
      using boost::math::tgamma;
      using boost::multiprecision::tgamma;

      auto autodiff_v = tgamma(make_fvar<T, m>(z));
      auto anchor_v = tgamma(z);
      BOOST_CHECK_CLOSE(autodiff_v.derivative(0u), anchor_v,
                        2.5e5 * test_constants::pct_epsilon());
    }

    {
      using boost::math::tgamma1pm1;
      auto autodiff_v = tgamma1pm1(make_fvar<T, m>(z));
      auto anchor_v = tgamma1pm1(z);
      BOOST_CHECK_CLOSE(autodiff_v.derivative(0u), anchor_v,
                        2.5e5 * test_constants::pct_epsilon());
    }

    {
      int s1 = 0;
      int s2 = 0;
      using boost::math::lgamma;
      using boost::multiprecision::lgamma;
      BOOST_CHECK_CLOSE(
          lgamma(make_fvar<T, m>(z).derivative(0u), std::addressof(s1)),
          lgamma(z, std::addressof(s2)), 2.5e5 * test_constants::pct_epsilon());
      BOOST_CHECK(
          (std::addressof(s1) == nullptr && std::addressof(s2) == nullptr) ||
          (s1 == s2));
    }

    {
      using boost::math::tgamma_lower;
      auto a = nextafter(a_sampler.next(), ((std::numeric_limits<T>::max))());
      auto autodiff_v = tgamma_lower(make_fvar<T, m>(a), make_fvar<T, m>(z));
      auto anchor_v = tgamma_lower(a, z);
      BOOST_CHECK_CLOSE(autodiff_v.derivative(0u), anchor_v,
                        2.5e5 * test_constants::pct_epsilon());
    }

    {
      using boost::math::gamma_q;
      auto a = nextafter(a_sampler.next(), ((std::numeric_limits<T>::max))());
      auto autodiff_v = gamma_q(make_fvar<T, m>(a), make_fvar<T, m>(z));
      auto anchor_v = gamma_q(a, z);
      BOOST_CHECK_CLOSE(autodiff_v.derivative(0u), anchor_v,
                        2.5e5 * test_constants::pct_epsilon());
    }
    {
      using boost::math::gamma_p;
      auto a = nextafter(a_sampler.next(), ((std::numeric_limits<T>::max))());
      auto autodiff_v = gamma_p(make_fvar<T, m>(a), make_fvar<T, m>(z));
      auto anchor_v = gamma_p(a, z);
      BOOST_CHECK_CLOSE(autodiff_v.derivative(0u), anchor_v,
                        2.5e5 * test_constants::pct_epsilon());
    }

    auto z_normalized = z_sampler.normalize(z);
    {
      using boost::math::gamma_p_inv;
      auto a_normalized = a_sampler.normalize(a_sampler.next());
      auto autodiff_v = gamma_p_inv(make_fvar<T, m>(a_normalized),
                                    make_fvar<T, m>(z_normalized));
      auto anchor_v = gamma_p_inv(a_normalized, z_normalized);
      BOOST_CHECK_CLOSE(autodiff_v.derivative(0u), anchor_v,
                        2.5e5 * test_constants::pct_epsilon());
    }
    {
      using boost::math::gamma_q_inv;
      auto a_normalized = a_sampler.normalize(a_sampler.next());
      auto autodiff_v = gamma_q_inv(make_fvar<T, m>(a_normalized),
                                    make_fvar<T, m>(z_normalized));
      auto anchor_v = gamma_q_inv(a_normalized, z_normalized);
      BOOST_CHECK_CLOSE(autodiff_v.derivative(0u), anchor_v,
                        2.5e5 * test_constants::pct_epsilon());
    }
    {
      using boost::math::gamma_p_inva;
      auto a_normalized = a_sampler.normalize(a_sampler.next());
      auto autodiff_v = gamma_p_inva(make_fvar<T, m>(a_normalized),
                                     make_fvar<T, m>(z_normalized));
      auto anchor_v = gamma_p_inva(a_normalized, z_normalized);
      BOOST_CHECK_CLOSE(autodiff_v.derivative(0u), anchor_v,
                        2.5e5 * test_constants::pct_epsilon());
    }
    {
      using boost::math::gamma_q_inva;
      auto a_normalized = a_sampler.normalize(a_sampler.next());
      auto autodiff_v = gamma_q_inva(make_fvar<T, m>(a_normalized),
                                     make_fvar<T, m>(z_normalized));
      auto anchor_v = gamma_q_inva(a_normalized, z_normalized);
      BOOST_CHECK_CLOSE(autodiff_v.derivative(0u), anchor_v,
                        2.5e5 * test_constants::pct_epsilon());
    }
    {
      using boost::math::gamma_p_derivative;
      auto a_normalized = a_sampler.normalize(a_sampler.next());
      auto autodiff_v = gamma_p_derivative(make_fvar<T, m>(a_normalized),
                                           make_fvar<T, m>(z_normalized));
      auto anchor_v = gamma_p_derivative(a_normalized, z_normalized);
      BOOST_CHECK_CLOSE(autodiff_v.derivative(0u), anchor_v,
                        2.5e5 * test_constants::pct_epsilon());
    }
  }
}

// Requires pow(complex<autodiff_fvar<T,m>>, T)
/* BOOST_AUTO_TEST_CASE_TEMPLATE(hankel_hpp, T, all_float_types) {
  using test_constants = test_constants_t<T>;
  static constexpr auto m = test_constants::order;
  test_detail::RandomSample<T> v_sampler{-200, 200};
  test_detail::RandomSample<T> x_sampler{-200, 200};
  for (auto i : boost::irange(test_constants::n_samples)) {
    std::ignore = i;
    auto v = v_sampler.next();
    auto x = x_sampler.next();

    try {
      auto autodiff_v = boost::math::cyl_hankel_1(make_fvar<T,m>(v),
make_fvar<T,m>(x)); auto anchor_v = boost::math::cyl_hankel_1(v, x);
      BOOST_CHECK_CLOSE(autodiff_v.real(), anchor_v.real(),
                                   test_constants::pct_epsilon());
      BOOST_CHECK_CLOSE(autodiff_v.imag(), anchor_v.imag(),
                                   test_constants::pct_epsilon());
    } catch (const std::domain_error &) {
      BOOST_CHECK_THROW(boost::math::cyl_hankel_1(make_fvar<T,m>(v),
make_fvar<T,m>(x)), boost::wrapexcept<std::domain_error>);
BOOST_CHECK_THROW(boost::math::cyl_hankel_1(v, x),
boost::wrapexcept<std::domain_error>); } catch (const std::overflow_error &) {
      BOOST_CHECK_THROW(boost::math::cyl_hankel_1(make_fvar<T,m>(v),
make_fvar<T,m>(x)), boost::wrapexcept<std::overflow_error>);
BOOST_CHECK_THROW(boost::math::cyl_hankel_1(v, x),
boost::wrapexcept<std::overflow_error>); } catch (...) { std::cout <<
std::setprecision(20) << "Input: x: " << x<< " max: "<<
std::numeric_limits<T>::max() << std::endl;
std::rethrow_exception(std::exception_ptr(std::current_exception()));
    }

    try {
      auto autodiff_v = boost::math::cyl_hankel_2(make_fvar<T,m>(v),
make_fvar<T,m>(x)); auto anchor_v = boost::math::cyl_hankel_2(v, x);
      BOOST_CHECK_CLOSE(autodiff_v.real(), anchor_v.real(),
                                   test_constants::pct_epsilon());
      BOOST_CHECK_CLOSE(autodiff_v.imag(), anchor_v.imag(),
                                   test_constants::pct_epsilon());
    } catch (const std::domain_error &) {
      BOOST_CHECK_THROW(boost::math::cyl_hankel_2(make_fvar<T,m>(v),
make_fvar<T,m>(x)), boost::wrapexcept<std::domain_error>);
BOOST_CHECK_THROW(boost::math::cyl_hankel_2(v, x),
boost::wrapexcept<std::domain_error>); } catch (const std::overflow_error &) {
      BOOST_CHECK_THROW(boost::math::cyl_hankel_2(make_fvar<T,m>(v),
make_fvar<T,m>(x)), boost::wrapexcept<std::overflow_error>);
BOOST_CHECK_THROW(boost::math::cyl_hankel_2(v, x),
boost::wrapexcept<std::overflow_error>); } catch (...) { std::cout <<
std::setprecision(20) << "Input: x: " << x<< " max: "<<
std::numeric_limits<T>::max() << std::endl;
std::rethrow_exception(std::exception_ptr(std::current_exception()));
    }

    try {
      auto autodiff_v = boost::math::sph_hankel_1(make_fvar<T,m>(v),
make_fvar<T,m>(x)); auto anchor_v = boost::math::sph_hankel_1(v, x);
      BOOST_CHECK_CLOSE(autodiff_v.derivative(0u), anchor_v,
                                   test_constants::pct_epsilon());
    } catch (const std::domain_error &) {
      BOOST_CHECK_THROW(boost::math::sph_hankel_1(make_fvar<T,m>(v),
make_fvar<T,m>(x)), boost::wrapexcept<std::domain_error>);
BOOST_CHECK_THROW(boost::math::sph_hankel_1(v, x),
boost::wrapexcept<std::domain_error>); } catch (const std::overflow_error &) {
      BOOST_CHECK_THROW(boost::math::sph_hankel_1(make_fvar<T,m>(v),
make_fvar<T,m>(x)), boost::wrapexcept<std::overflow_error>);
BOOST_CHECK_THROW(boost::math::sph_hankel_1(v, x),
boost::wrapexcept<std::overflow_error>); } catch (...) { std::cout <<
std::setprecision(20) << "Input: x: " << x<< " max: "<<
std::numeric_limits<T>::max() << std::endl;
std::rethrow_exception(std::exception_ptr(std::current_exception()));
    }

    try {
      auto autodiff_v = boost::math::sph_hankel_2(make_fvar<T,m>(v),
make_fvar<T,m>(x)); auto anchor_v = boost::math::sph_hankel_2(v, x);
      BOOST_CHECK_CLOSE(autodiff_v.derivative(0u), anchor_v,
                                   test_constants::pct_epsilon());
    } catch (const std::domain_error &) {
      BOOST_CHECK_THROW(boost::math::sph_hankel_2(make_fvar<T,m>(v),
make_fvar<T,m>(x)), boost::wrapexcept<std::domain_error>);
BOOST_CHECK_THROW(boost::math::sph_hankel_2(v, x),
boost::wrapexcept<std::domain_error>); } catch (const std::overflow_error &) {
      BOOST_CHECK_THROW(boost::math::sph_hankel_2(make_fvar<T,m>(v),
make_fvar<T,m>(x)), boost::wrapexcept<std::overflow_error>);
BOOST_CHECK_THROW(boost::math::sph_hankel_2(v, x),
boost::wrapexcept<std::overflow_error>); } catch (...) { std::cout <<
std::setprecision(20) << "Input: x: " << x<< " max: "<<
std::numeric_limits<T>::max() << std::endl;
std::rethrow_exception(std::exception_ptr(std::current_exception()));
    }
  }
} */

BOOST_AUTO_TEST_SUITE_END()
