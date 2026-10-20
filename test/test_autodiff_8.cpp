//           Copyright Matthew Pulver 2018 - 2019.
// Distributed under the Boost Software License, Version 1.0.
//      (See accompanying file LICENSE_1_0.txt or copy at
//           https://www.boost.org/LICENSE_1_0.txt)

#include "test_autodiff.hpp"
#include <boost/math/special_functions.hpp>


template <class T>
void hermite_hpp() {
  using test_constants = test_constants_t<T>;
  static constexpr auto m = test_constants::order;
  test_detail::RandomSample<T> x_sampler{-1, 1};
  for (auto i : boost::irange(3u)) {
    auto x = x_sampler.next();
    auto autodiff_v = boost::math::hermite(i, make_fvar<T, m>(x));
    auto anchor_v = boost::math::hermite(i, x);
    CHECK_RELATIVE_ERROR(autodiff_v.derivative(0u), anchor_v, (100 * test_constants::pct_epsilon()) / 100);
  }
}

template <class T>
void heuman_lambda_hpp() {
  using test_constants = test_constants_t<T>;
  static constexpr auto m = test_constants::order;
  test_detail::RandomSample<T> x_sampler{-1, 1};
  test_detail::RandomSample<T> phi_sampler{-boost::math::constants::two_pi<T>(),
                                           boost::math::constants::two_pi<T>()};
  for (auto i : boost::irange(test_constants::n_samples)) {
    std::ignore = i;
    auto x = x_sampler.next();
    auto phi = phi_sampler.next();
    auto autodiff_v =
        boost::math::heuman_lambda(make_fvar<T, m>(x), make_fvar<T, m>(phi));
    auto anchor_v = boost::math::heuman_lambda(x, phi);
    CHECK_TRUE(isNearZero(autodiff_v.derivative(0u) - anchor_v));
  }
}

template <class T>
void hypot_hpp() {
  using test_constants = test_constants_t<T>;
  static constexpr auto m = test_constants::order;
  test_detail::RandomSample<T> x_sampler{-2000, 2000};
  test_detail::RandomSample<T> y_sampler{-2000, 2000};
  for (auto i : boost::irange(test_constants::n_samples)) {
    std::ignore = i;
    auto x = x_sampler.next();
    auto y = y_sampler.next();
    auto autodiff_v =
        boost::math::hypot(make_fvar<T, m>(x), make_fvar<T, m>(y));
    auto anchor_v = boost::math::hypot(x, y);
    CHECK_TRUE(isNearZero(autodiff_v.derivative(0u) - anchor_v));
  }
}

template <class T>
void jacobi_elliptic_hpp() {
  using test_constants = test_constants_t<T>;
  static constexpr auto m = test_constants::order;
  test_detail::RandomSample<T> k_sampler{0, 1};
  test_detail::RandomSample<T> theta_sampler{-100, 100};
  for (auto i : boost::irange(test_constants::n_samples)) {
    std::ignore = i;
    auto k = k_sampler.next();
    auto theta = theta_sampler.next();
    CHECK_TRUE(isNearZero(
        boost::math::jacobi_cd(make_fvar<T, m>(k), make_fvar<T, m>(theta))
            .derivative(0u) -
        boost::math::jacobi_cd(k, theta)));
    CHECK_TRUE(isNearZero(
        boost::math::jacobi_cn(make_fvar<T, m>(k), make_fvar<T, m>(theta))
            .derivative(0u) -
        boost::math::jacobi_cn(k, theta)));
    CHECK_TRUE(isNearZero(
        boost::math::jacobi_cs(make_fvar<T, m>(k), make_fvar<T, m>(theta))
            .derivative(0u) -
        boost::math::jacobi_cs(k, theta)));
    CHECK_TRUE(isNearZero(
        boost::math::jacobi_dc(make_fvar<T, m>(k), make_fvar<T, m>(theta))
            .derivative(0u) -
        boost::math::jacobi_dc(k, theta)));
    CHECK_TRUE(isNearZero(
        boost::math::jacobi_dn(make_fvar<T, m>(k), make_fvar<T, m>(theta))
            .derivative(0u) -
        boost::math::jacobi_dn(k, theta)));
    CHECK_TRUE(isNearZero(
        boost::math::jacobi_ds(make_fvar<T, m>(k), make_fvar<T, m>(theta))
            .derivative(0u) -
        boost::math::jacobi_ds(k, theta)));
    CHECK_TRUE(isNearZero(
        boost::math::jacobi_nc(make_fvar<T, m>(k), make_fvar<T, m>(theta))
            .derivative(0u) -
        boost::math::jacobi_nc(k, theta)));
    CHECK_TRUE(isNearZero(
        boost::math::jacobi_nd(make_fvar<T, m>(k), make_fvar<T, m>(theta))
            .derivative(0u) -
        boost::math::jacobi_nd(k, theta)));
    CHECK_TRUE(isNearZero(
        boost::math::jacobi_ns(make_fvar<T, m>(k), make_fvar<T, m>(theta))
            .derivative(0u) -
        boost::math::jacobi_ns(k, theta)));
    CHECK_TRUE(isNearZero(
        boost::math::jacobi_sc(make_fvar<T, m>(k), make_fvar<T, m>(theta))
            .derivative(0u) -
        boost::math::jacobi_sc(k, theta)));
    CHECK_TRUE(isNearZero(
        boost::math::jacobi_sd(make_fvar<T, m>(k), make_fvar<T, m>(theta))
            .derivative(0u) -
        boost::math::jacobi_sd(k, theta)));
    CHECK_TRUE(isNearZero(
        boost::math::jacobi_sn(make_fvar<T, m>(k), make_fvar<T, m>(theta))
            .derivative(0u) -
        boost::math::jacobi_sn(k, theta)));
  }
}

template <class T>
void jacobi_zeta_hpp() {
  using test_constants = test_constants_t<T>;
  static constexpr auto m = test_constants::order;
  test_detail::RandomSample<T> x_sampler{-1, 1};
  test_detail::RandomSample<T> phi_sampler{-boost::math::constants::two_pi<T>(),
                                           boost::math::constants::two_pi<T>()};
  for (auto i : boost::irange(test_constants::n_samples)) {
    std::ignore = i;
    auto x = x_sampler.next();
    auto phi = phi_sampler.next();
    CHECK_TRUE(isNearZero(
        boost::math::jacobi_zeta(make_fvar<T, m>(x), make_fvar<T, m>(phi))
            .derivative(0u) -
        boost::math::jacobi_zeta(x, phi)));
  }
}

template <class T>
void laguerre_hpp() {
  using test_constants = test_constants_t<T>;
  static constexpr auto m = test_constants::order;
  int n = 3;
  int r = 2;
  for (T x : {T(0), T(0.25), T(0.5), T(0.75), T(1)}) {
    {
      auto autodiff_v = boost::math::laguerre(n, make_fvar<T, m>(x));
      auto anchor_v = boost::math::laguerre(n, x);

      CHECK_RELATIVE_ERROR(autodiff_v.derivative(0u), anchor_v, (1000 * test_constants::pct_epsilon()) / 100);
    }

    {
      auto autodiff_v = boost::math::laguerre(n, r, make_fvar<T, m>(x));
      auto anchor_v = boost::math::laguerre(n, r, x);

      CHECK_RELATIVE_ERROR(autodiff_v.derivative(0u), anchor_v, (500 * test_constants::pct_epsilon()) / 100);
    }
  }
}

template <class T>
void lambert_w_hpp() {
  using boost::math::nextafter;
  using boost::math::tools::max;
  using boost::multiprecision::fabs;
  using boost::multiprecision::min;
  using detail::fabs;
  using std::fabs;
  using std::max;
  using std::min;
  using std::nextafter;

  using promoted_t = promote<T, double>;
  using test_constants = test_constants_t<T>;
  static constexpr auto m = test_constants::order;
  test_detail::RandomSample<T> x_sampler{static_cast<T>(-1 / std::exp(-1)),
                                         ((std::numeric_limits<T>::max))()};
  for (auto i : boost::irange(test_constants::n_samples)) {
    std::ignore = i;
    auto x = x_sampler.next();
    {
      auto x_ = (min<T>)(static_cast<T>(((max<promoted_t>))(
                             -exp(promoted_t(-1)), promoted_t(x))),
                         ((std::numeric_limits<T>::max))());
      {
        auto autodiff_v = boost::math::lambert_w0(make_fvar<T, m>(x_));
        auto anchor_v = boost::math::lambert_w0(x_);
        CHECK_TRUE(isNearZero(autodiff_v.derivative(0u) - anchor_v));
      }
      {
        auto autodiff_v = boost::math::lambert_w0_prime(make_fvar<T, m>(x_));
        auto anchor_v = boost::math::lambert_w0_prime(x_);
        CHECK_TRUE(isNearZero(autodiff_v.derivative(0u) - anchor_v));
      }
    }

    {
      auto x_ = nextafter(
          static_cast<T>(nextafter(
              ((max))(static_cast<promoted_t>(-exp(-1)), -fabs(promoted_t(x))),
              ((std::numeric_limits<promoted_t>::max))())),
          ((std::numeric_limits<T>::max))());
      x_ = (max)(static_cast<T>(-0.36), x_);
      CHECK_TRUE(isNearZero(
          boost::math::lambert_wm1(make_fvar<T, m>(x_)).derivative(0u) -
          boost::math::lambert_wm1(x_)));
      CHECK_TRUE(isNearZero(
          boost::math::lambert_wm1_prime(make_fvar<T, m>(x_)).derivative(0u) -
          boost::math::lambert_wm1_prime(x_)));
    }
  }
}

template <class T>
void log1p_hpp() {
  using boost::math::log1p;
  using boost::multiprecision::log1p;

  using test_constants = test_constants_t<T>;
  static constexpr auto m = test_constants::order;
  test_detail::RandomSample<T> x_sampler{-1, 2000};
  for (auto i : boost::irange(test_constants::n_samples)) {
    std::ignore = i;
    auto x = x_sampler.next();
    CHECK_TRUE(
        isNearZero(log1p(make_fvar<T, m>(x)).derivative(0u) - log1p(x)));
  }
}

template <class T>
void next_hpp() {
  using boost::math::float_advance;
  using boost::math::float_distance;
  using boost::math::float_next;
  using boost::math::float_prior;
  using boost::math::nextafter;
  using boost::multiprecision::nextafter;

  using test_constants = test_constants_t<T>;
  static constexpr auto m = test_constants::order;
  for (auto i : boost::irange(test_constants::n_samples)) {
    const auto j = static_cast<T>(i);
    auto fvar_j = make_fvar<T, m>(j);
    CHECK_TRUE(isNearZero(float_next(fvar_j).derivative(0u) - float_next(j)));
    CHECK_TRUE(
        isNearZero(float_prior(fvar_j).derivative(0u) - float_prior(j)));

    CHECK_TRUE(isNearZero(
        nextafter(fvar_j, make_fvar<T, m>(static_cast<T>(1))).derivative(0u) -
        nextafter(j, static_cast<T>(1))));
    CHECK_TRUE(
        isNearZero(nextafter(fvar_j, make_fvar<T, m>(static_cast<T>(i + 2))) -
                   nextafter(j, static_cast<T>(i + 2))));
    CHECK_TRUE(
        isNearZero(nextafter(fvar_j, make_fvar<T, m>(static_cast<T>(i + 1)))
                       .derivative(0u) -
                   nextafter(j, static_cast<T>(i + 2))));
    CHECK_TRUE(isNearZero(
        nextafter(fvar_j, make_fvar<T, m>(static_cast<T>(-1))).derivative(0u) -
        nextafter(j, static_cast<T>(-1))));
    CHECK_TRUE(isNearZero(
        nextafter(fvar_j, make_fvar<T, m>(static_cast<T>(-1 * (i + 2))))
            .derivative(0u) -
        nextafter(j, -1 * static_cast<T>(i + 2))));

    CHECK_TRUE(isNearZero(
        nextafter(fvar_j, make_fvar<T, m>(static_cast<T>(-1 * (i + 1))))
            .derivative(0u) -
        nextafter(j, -1 * static_cast<T>(i + 2))));

    CHECK_TRUE(isNearZero(nextafter(fvar_j, fvar_j) - fvar_j));

    CHECK_TRUE(isNearZero(float_advance(fvar_j, 1).derivative(0u) -
                           float_advance(j, 1)));
    CHECK_TRUE(isNearZero(float_advance(fvar_j, i + 2).derivative(0u) -
                           float_advance(j, i + 2)));
    CHECK_TRUE(isNearZero(
        float_advance(fvar_j, i + 1).derivative(0u) -
        float_advance(float_advance(fvar_j, i + 2), -1).derivative(0u)));

    CHECK_TRUE(isNearZero(float_advance(fvar_j, -1).derivative(0u) -
                           float_advance(j, -1)));

    CHECK_TRUE(isNearZero(
        float_advance(fvar_j, -i - 1).derivative(0u) -
        float_advance(float_advance(fvar_j, -i - 2), 1).derivative(0u)));

    CHECK_TRUE(isNearZero(float_advance(fvar_j, 0) - fvar_j));

    CHECK_TRUE(isNearZero(float_distance(fvar_j, j).derivative(0u) -
                           static_cast<T>(0)));
    CHECK_TRUE(
        isNearZero(float_distance(float_next(fvar_j), fvar_j).derivative(0u) -
                   ((make_fvar<T, m>(-1))).derivative(0u)));
    CHECK_TRUE(
        isNearZero(float_distance(float_prior(fvar_j), fvar_j).derivative(0u) -
                   (make_fvar<T, m>(1)).derivative(0u)));
  }
}

template <class T>
void owens_t_hpp() {
  BOOST_MATH_STD_USING;
  using test_constants = test_constants_t<T>;
  static constexpr auto m = test_constants::order;
  test_detail::RandomSample<T> h_sampler{-2000, 2000};
  test_detail::RandomSample<T> a_sampler{-2000, 2000};
  for (auto i : boost::irange(test_constants::n_samples)) {
    std::ignore = i;
    auto h = h_sampler.next();
    auto a = a_sampler.next();
    auto autodiff_v =
        boost::math::owens_t(make_fvar<T, m>(h), make_fvar<T, m>(a));
    auto anchor_v = boost::math::owens_t(h, a);
    CHECK_TRUE(isNearZero(autodiff_v.derivative(0u) - anchor_v));
  }
}


template <class T>
void pow_hpp() {
  using test_constants = test_constants_t<T>;
  static constexpr auto m = test_constants::order;
  for (auto i : boost::irange(10)) {
    CHECK_RELATIVE_ERROR(
        boost::math::pow<0>(make_fvar<T, m>(static_cast<T>(i))).derivative(0u),
        boost::math::pow<0>(static_cast<T>(i)), (test_constants::pct_epsilon()) / 100);
    CHECK_RELATIVE_ERROR(
        boost::math::pow<1>(make_fvar<T, m>(static_cast<T>(i))).derivative(0u),
        boost::math::pow<1>(static_cast<T>(i)), (test_constants::pct_epsilon()) / 100);
    CHECK_RELATIVE_ERROR(
        boost::math::pow<2>(make_fvar<T, m>(static_cast<T>(i))).derivative(0u),
        boost::math::pow<2>(static_cast<T>(i)), (test_constants::pct_epsilon()) / 100);
    CHECK_RELATIVE_ERROR(
        boost::math::pow<3>(make_fvar<T, m>(static_cast<T>(i))).derivative(0u),
        boost::math::pow<3>(static_cast<T>(i)), (test_constants::pct_epsilon()) / 100);
    CHECK_RELATIVE_ERROR(
        boost::math::pow<4>(make_fvar<T, m>(static_cast<T>(i))).derivative(0u),
        boost::math::pow<4>(static_cast<T>(i)), (test_constants::pct_epsilon()) / 100);
    CHECK_RELATIVE_ERROR(
        boost::math::pow<5>(make_fvar<T, m>(static_cast<T>(i))).derivative(0u),
        boost::math::pow<5>(static_cast<T>(i)), (test_constants::pct_epsilon()) / 100);
    CHECK_RELATIVE_ERROR(
        boost::math::pow<6>(make_fvar<T, m>(static_cast<T>(i))).derivative(0u),
        boost::math::pow<6>(static_cast<T>(i)), (test_constants::pct_epsilon()) / 100);
    CHECK_RELATIVE_ERROR(
        boost::math::pow<7>(make_fvar<T, m>(static_cast<T>(i))).derivative(0u),
        boost::math::pow<7>(static_cast<T>(i)), (test_constants::pct_epsilon()) / 100);
    CHECK_RELATIVE_ERROR(
        boost::math::pow<8>(make_fvar<T, m>(static_cast<T>(i))).derivative(0u),
        boost::math::pow<8>(static_cast<T>(i)), (test_constants::pct_epsilon()) / 100);
    CHECK_RELATIVE_ERROR(
        boost::math::pow<9>(make_fvar<T, m>(static_cast<T>(i))).derivative(0u),
        boost::math::pow<9>(static_cast<T>(i)), (test_constants::pct_epsilon()) / 100);
  }
}

template <class T>
void polygamma_hpp() {
  using test_constants = test_constants_t<T>;
  static constexpr auto m = test_constants::order;
  test_detail::RandomSample<T> x_sampler{0, 2000};
  for (auto i : boost::irange(test_constants::n_samples)) {
    auto x = x_sampler.next();
    try {
      auto autodiff_v = boost::math::polygamma(i, make_fvar<T, m>(x));
      auto anchor_v = boost::math::polygamma(i, x);
      CHECK_TRUE(isNearZero(autodiff_v.derivative(0u) - anchor_v));
    } catch (const std::overflow_error &) {
      std::cout << "Overflow Error thrown with inputs i: " << i << " x: " << x
                << std::endl;
      CHECK_THROW(boost::math::polygamma(i, make_fvar<T, m>(x)),
                        boost::wrapexcept<std::overflow_error>);
      CHECK_THROW(boost::math::polygamma(i, x),
                        boost::wrapexcept<std::overflow_error>);
    }
  }
}

template <class T>
void powm1_hpp() {
  using boost::math::tools::max;
  using boost::multiprecision::max;
  using std::max;

  using boost::multiprecision::log;
  using detail::log;
  using std::log;

  using boost::multiprecision::min;
  using std::min;

  using boost::multiprecision::sqrt;
  using detail::sqrt;
  using std::sqrt;

  using boost::math::nextafter;
  using boost::multiprecision::nextafter;

  using test_constants = test_constants_t<T>;
  static constexpr auto m = test_constants::order;
  test_detail::RandomSample<T> x_sampler{0, 2000};
  for (auto i : boost::irange(test_constants::n_samples)) {
    std::ignore = i;
    auto x = ((max))(x_sampler.next(),
                     boost::math::nextafter(static_cast<T>(0),
                                            ((std::numeric_limits<T>::max))()));

    auto y =
        ((min))(x_sampler.next(),
                log(sqrt(((std::numeric_limits<T>::max))()) + 1) / log(x + 1));
    auto autodiff_v =
        boost::math::powm1(make_fvar<T, m>(x), make_fvar<T, m>(y));
    auto anchor_v = boost::math::powm1(x, y);
    CHECK_TRUE(isNearZero(autodiff_v.derivative(0u) - anchor_v));
  }
}

#if __clang_major__ > 5 || __GNUC__ > 5 || defined(_MSC_VER)

template <class T>
void sin_pi_hpp() {
  using test_constants = test_constants_t<T>;
  static constexpr auto m = test_constants::order;
  test_detail::RandomSample<T> x_sampler{-2000, 2000};
  for (auto i : boost::irange(test_constants::n_samples)) {
    std::ignore = i;
    auto x = x_sampler.next();
    CHECK_TRUE(
        isNearZero(boost::math::sin_pi(make_fvar<T, m>(x)).derivative(0u) -
                   boost::math::sin_pi(x)));
  }
}

template <class T>
void sinhc_hpp() {
  using test_constants = test_constants_t<T>;
  static constexpr auto m = test_constants::order;
  test_detail::RandomSample<T> x_sampler{ -80, 80 };
  for (auto i : boost::irange(test_constants::n_samples)) {
    std::ignore = i;
    auto x = x_sampler.next();
    if (x != 0) {
      auto autodiff_v = boost::math::sinhc_pi(make_fvar<T, m>(x));
      auto anchor_v = boost::math::sinhc_pi(x);
      CHECK_RELATIVE_ERROR(autodiff_v.derivative(0u), anchor_v, (50 * test_constants::pct_epsilon()) / 100);
    }
  }
}

template <class T>
void spherical_harmonic_hpp() {
  using test_constants = test_constants_t<T>;
  static constexpr auto m = test_constants::order;
  test_detail::RandomSample<T> theta_sampler{0,
                                             boost::math::constants::pi<T>()};
  test_detail::RandomSample<T> phi_sampler{0,
                                           boost::math::constants::two_pi<T>()};
  test_detail::RandomSample<int> r_sampler{0, test_constants::n_samples};
  for (auto n : boost::irange(
           1u, static_cast<unsigned>(test_constants::n_samples) + 1u)) {
    auto theta = theta_sampler.next();
    auto phi = phi_sampler.next();
    auto r = (std::min)(static_cast<int>(n) - 1, r_sampler.next());
    {
      auto autodiff_v = boost::math::spherical_harmonic(
          n, r, make_fvar<T, m>(theta), make_fvar<T, m>(phi));
      auto anchor_v = boost::math::spherical_harmonic(n, r, theta, phi);
      CHECK_TRUE(
          isNearZero(autodiff_v.real().derivative(0u) - anchor_v.real()));
      CHECK_TRUE(
          isNearZero(autodiff_v.imag().derivative(0u) - anchor_v.imag()));
    }

    {
      auto autodiff_v = boost::math::spherical_harmonic_r(
          n, r, make_fvar<T, m>(theta), make_fvar<T, m>(phi));
      auto anchor_v = boost::math::spherical_harmonic_r(n, r, theta, phi);
      CHECK_TRUE(isNearZero(autodiff_v.derivative(0u) - anchor_v));
    }

    {
      auto autodiff_v = boost::math::spherical_harmonic_i(
          n, r, make_fvar<T, m>(theta), make_fvar<T, m>(phi));
      auto anchor_v = boost::math::spherical_harmonic_i(n, r, theta, phi);
      CHECK_TRUE(isNearZero(autodiff_v.derivative(0u) - anchor_v));
    }
  }
}

template <class T>
void sqrt1pm1_hpp() {
  using test_constants = test_constants_t<T>;
  static constexpr auto m = test_constants::order;
  test_detail::RandomSample<T> x_sampler{-1, 2000};
  for (auto i : boost::irange(test_constants::n_samples)) {
    std::ignore = i;
    auto x = x_sampler.next();
    CHECK_TRUE(
        isNearZero(boost::math::sqrt1pm1(make_fvar<T, m>(x)).derivative(0u) -
                   boost::math::sqrt1pm1(x)));
  }
}

template <class T>
void trigamma_hpp() {
  using test_constants = test_constants_t<T>;
  static constexpr auto m = test_constants::order;
  test_detail::RandomSample<T> x_sampler{0, 2000};
  for (auto i : boost::irange(test_constants::n_samples)) {
    std::ignore = i;
    auto x = x_sampler.next();
    CHECK_TRUE(
        isNearZero(boost::math::trigamma(make_fvar<T, m>(x)).derivative(0u) -
                   boost::math::trigamma(x)));
  }
}

#endif // Compiler guard

int main()
{
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { hermite_hpp<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { heuman_lambda_hpp<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { hypot_hpp<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { jacobi_elliptic_hpp<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { jacobi_zeta_hpp<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { laguerre_hpp<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { lambert_w_hpp<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { log1p_hpp<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { next_hpp<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, bin_float_types>>([](auto t) { owens_t_hpp<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { pow_hpp<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { polygamma_hpp<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { powm1_hpp<typename decltype(t)::type>(); });
#if __clang_major__ > 5 || __GNUC__ > 5 || defined(_MSC_VER)
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { sin_pi_hpp<typename decltype(t)::type>(); });
#endif
#if __clang_major__ > 5 || __GNUC__ > 5 || defined(_MSC_VER)
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { sinhc_hpp<typename decltype(t)::type>(); });
#endif
#if __clang_major__ > 5 || __GNUC__ > 5 || defined(_MSC_VER)
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { spherical_harmonic_hpp<typename decltype(t)::type>(); });
#endif
#if __clang_major__ > 5 || __GNUC__ > 5 || defined(_MSC_VER)
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { sqrt1pm1_hpp<typename decltype(t)::type>(); });
#endif
#if __clang_major__ > 5 || __GNUC__ > 5 || defined(_MSC_VER)
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { trigamma_hpp<typename decltype(t)::type>(); });
#endif
   return boost::math::test::report_errors();
}
