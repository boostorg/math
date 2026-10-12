//           Copyright Matthew Pulver 2018 - 2019.
// Distributed under the Boost Software License, Version 1.0.
//      (See accompanying file LICENSE_1_0.txt or copy at
//           https://www.boost.org/LICENSE_1_0.txt)

#include "test_autodiff.hpp"
#include <boost/math/tools/test_value.hpp>


template <class T>
void lround_llround_lltrunc_truncl() {
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
  CHECK_EQUAL(yl, lround(cx));
  auto yll = llround(x);
  CHECK_EQUAL(yll, llround(cx));
  CHECK_EQUAL(lltrunc(cx), lltrunc(x));

#ifndef BOOST_NO_CXX17_IF_CONSTEXPR
  if constexpr (!bmp::is_number<T>::value &&
                !bmp::is_number_expression<T>::value) {
    CHECK_EQUAL(truncl(x), truncl(cx));
  }
#endif
}

template <class T>
void equality() {
  BOOST_MATH_STD_USING
  using boost::math::epsilon_difference;
  using boost::math::fpclassify;
  using boost::math::ulp;
  using std::fpclassify;

  constexpr std::size_t m = 3;
  // check zeros
  {
    auto x = make_fvar<T, m>(T(0));
    auto y = T(-0.0);
    CHECK_EQUAL(x.derivative(0u), y);
  }
}

#if defined(BOOST_AUTODIFF_TESTING_INCLUDE_MULTIPRECISION)
template <class T>
void multiprecision() {
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
  const T answer = BOOST_MATH_TEST_VALUE(T, 1976.31960074779771777988187529041872090812118921875499076582535951111845769110560421820940516423255314);
  // BOOST_CHECK_CLOSE(v.derivative(Nw,Nx,Ny,Nz), answer, eps); // Doesn't work
  // for cpp_dec_float
  const T relative_error =
      static_cast<T>(fabs(v.derivative(Nw, Nx, Ny, Nz) / answer - 1));
  CHECK_TRUE((relative_error) < (eps));
}
#endif

template <class T>
void acosh_hpp() {
  using boost::math::acosh;
  using test_constants = test_constants_t<T>;
  static constexpr auto m = test_constants::order;

  test_detail::RandomSample<T> x_sampler{1, 100};
  for (auto i : boost::irange(test_constants::n_samples)) {
    std::ignore = i;
    auto x = x_sampler.next();
    auto autodiff_v = acosh(make_fvar<T, m>(x));
    auto anchor_v = acosh(x);
    CHECK_RELATIVE_ERROR(autodiff_v.derivative(0u), anchor_v, (1e3 * test_constants::pct_epsilon()) / 100);
  }
}

template <class T>
void asinh_hpp() {
  using boost::math::asinh;
  using test_constants = test_constants_t<T>;
  static constexpr auto m = test_constants::order;

  test_detail::RandomSample<T> x_sampler{-100, 100};
  for (auto i : boost::irange(test_constants::n_samples)) {
    std::ignore = i;
    auto x = x_sampler.next();

    auto autodiff_v = asinh(make_fvar<T, m>(x));
    auto anchor_v = asinh(x);
    CHECK_RELATIVE_ERROR(autodiff_v.derivative(0u), anchor_v, (1e3 * test_constants::pct_epsilon()) / 100);
  }
}

template <class T>
void atanh_hpp() {
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
    CHECK_RELATIVE_ERROR(autodiff_v.derivative(0u), anchor_v, (1e3 * test_constants::pct_epsilon()) / 100);
  }
}

template <class T>
void atan_hpp() {
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
    CHECK_RELATIVE_ERROR(autodiff_v.derivative(0u), anchor_v, (T(1e3) * test_constants::pct_epsilon()) / 100);
  }
}

template <class T>
void bernoulli_hpp() {

  using boost::multiprecision::min;
  using std::min;
  using test_constants = test_constants_t<T>;
  static constexpr auto m = test_constants::order;

  for (auto i : boost::irange(test_constants::n_samples)) {
    {
      auto autodiff_v = boost::math::bernoulli_b2n<autodiff_fvar<T, m>>(i);
      auto anchor_v = boost::math::bernoulli_b2n<T>(i);
      CHECK_RELATIVE_ERROR(autodiff_v.derivative(0u), anchor_v, (50 * test_constants::pct_epsilon()) / 100);
    }
    {
      auto i_ = (min)(19, i);
      auto autodiff_v = boost::math::tangent_t2n<autodiff_fvar<T, m>>(i_);
      auto anchor_v = boost::math::tangent_t2n<T>(i_);
      CHECK_RELATIVE_ERROR(autodiff_v.derivative(0u), anchor_v, (50 * test_constants::pct_epsilon()) / 100);
    }
  }
}

int main()
{
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { lround_llround_lltrunc_truncl<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { equality<typename decltype(t)::type>(); });
#if defined(BOOST_AUTODIFF_TESTING_INCLUDE_MULTIPRECISION)
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, multiprecision_float_types>>([](auto t) { multiprecision<typename decltype(t)::type>(); });
#endif
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { acosh_hpp<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { asinh_hpp<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { atanh_hpp<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { atan_hpp<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { bernoulli_hpp<typename decltype(t)::type>(); });
   return boost::math::test::report_errors();
}
