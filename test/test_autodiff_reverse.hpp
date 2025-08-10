#ifndef TEST_AUTODIFF_REVERSE_HPP
#define TEST_AUTODIFF_REVERSE_HPP

#ifndef BOOST_TEST_MODULE
#define BOOST_TEST_MODULE test_autodiff
#endif

#include <boost/test/included/unit_test.hpp>

#include <algorithm>
#include <boost/math/differentiation/autodiff_reverse.hpp>
#include <boost/math/tools/config.hpp>
#include <boost/mp11/function.hpp>
#include <boost/mp11/integral.hpp>
#include <boost/mp11/list.hpp>
#include <boost/mp11/utility.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/test/included/unit_test.hpp>
#include <cfenv>
#include <cstdlib>
#include <random>

namespace mp11         = boost::mp11;
namespace bmp          = boost::multiprecision;
namespace rdiff_detail = boost::math::differentiation::reverse_mode::detail;
namespace rdiff        = boost::math::differentiation::reverse_mode;

#if defined(BOOST_USE_VALGRIND) || defined(BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS)
using bin_float_types = mp11::mp_list<float>;
#elif defined(__STDCPP_FLOAT32_T__) && defined(__STDCPP_FLOAT64_T__)
using bin_float_types = mp11::mp_list<std::float32_t, std::float64_t>;
#else
using bin_float_types = mp11::mp_list<float, double, long double>;
#endif

#if !defined(BOOST_VERSION) || BOOST_VERSION < 107000 || defined(BOOST_USE_VALGRIND) \
    || defined(BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS) || defined(BOOST_NO_STRESS_TEST) \
    || defined(BOOST_MATH_STANDALONE)
using multiprecision_float_types = mp11::mp_list<>;
#else
#define BOOST_AUTODIFF_TESTING_INCLUDE_MULTIPRECISION
using multiprecision_float_types = mp11::mp_list<bmp::cpp_bin_float_50>;
#endif

using all_float_types
    = bin_float_types; //mp11::mp_append<bin_float_types, multiprecision_float_types>;

using namespace boost::math::differentiation;
#endif // TEST_AUTODIFF_REVERSE_HPP
template<typename T>
using is_multiprecision_t = mp11::mp_or<bmp::is_number<T>, bmp::is_number_expression<T>>;

template<bool IfValue, typename ThenType, typename ElseType>
using if_c = mp11::mp_eval_if_c<IfValue, ThenType, mp11::mp_identity_t, ElseType>;

template<typename IfType, typename ThenType, typename ElseType>
using if_t = if_c<IfType::value, ThenType, ElseType>;
/**
 * struct to emit pseudo-random values from a given interval.
 * Endpoints are closed or open depending on whether or not they're infinite).
 */
/**
 * Simple struct to hold constants that are used in each test
 * since BOOST_AUTO_TEST_CASE_TEMPLATE doesn't support fixtures.
 */
template<typename T, std::size_t OrderValue>
struct test_constants_t
{
    static constexpr auto n_samples
        = if_t<mp11::mp_or<bmp::is_number<T>, bmp::is_number_expression<T>>,
               mp11::mp_int<10>,
               mp11::mp_int<25>>::value;
    static constexpr auto order = OrderValue;
    static constexpr T    pct_epsilon() noexcept
    {
        return (is_multiprecision_t<T>::value ? 2 : 1) * std::numeric_limits<T>::epsilon() * 100;
    }
};

template<typename T, std::size_t Order = 5>
using test_constants = test_constants_t<T, Order>;

template<typename T>
struct RandomSample
{
    using numeric_limits_t = std::numeric_limits<T>;
    using is_integer_t     = mp11::mp_bool<std::numeric_limits<T>::is_integer>;

    using distribution_param_t
        = if_t<is_multiprecision_t<T>,
               if_t<is_integer_t, if_c<numeric_limits_t::is_signed, int64_t, uint64_t>, long double>,
               T>;
    static_assert((std::numeric_limits<T>::is_integer
                   && std::numeric_limits<distribution_param_t>::is_integer)
                      || (!std::numeric_limits<T>::is_integer
                          && !std::numeric_limits<distribution_param_t>::is_integer),
                  "T and distribution_param_t must either both be integral or "
                  "both be not integral");

    using dist_t = if_t<is_integer_t,
                        std::uniform_int_distribution<distribution_param_t>,
                        std::uniform_real_distribution<distribution_param_t>>;

    struct get_integral_endpoint
    {
        template<typename V>
        constexpr distribution_param_t operator()(V finish) const noexcept
        {
            return static_cast<distribution_param_t>(finish);
        }
    };

    struct get_real_endpoint
    {
        template<typename V>
        constexpr distribution_param_t operator()(V finish) const noexcept
        {
            return std::nextafter(static_cast<distribution_param_t>(finish),
                                  (std::numeric_limits<distribution_param_t>::max)());
        }
    };

    using get_endpoint_t = if_t<is_integer_t, get_integral_endpoint, get_real_endpoint>;

    template<typename U, typename V>
    RandomSample(U start, V finish)
        : rng_(std::random_device{}())
        , dist_(static_cast<distribution_param_t>(start), get_endpoint_t{}(finish))
    {}

    T            next() noexcept { return static_cast<T>(dist_(rng_)); }
    T            normalize(const T& x) noexcept { return x / ((dist_.max)() - (dist_.min)()); }

    std::mt19937 rng_;
    dist_t       dist_;
};
static_assert(
    std::is_same<RandomSample<float>::dist_t, std::uniform_real_distribution<float>>::value, "");
static_assert(
    std::is_same<RandomSample<int64_t>::dist_t, std::uniform_int_distribution<int64_t>>::value, "");
static_assert(std::is_same<RandomSample<bmp::uint512_t>::dist_t,
                           std::uniform_int_distribution<uint64_t>>::value,
              "");
static_assert(std::is_same<RandomSample<bmp::cpp_bin_float_50>::dist_t,
                           std::uniform_real_distribution<long double>>::value,
              "");

template<typename T>
constexpr T boost_close_tol(T scale_factor = 1e5)
{
    static_assert(std::is_floating_point<T>::value, "T must be floating point");
    return std::numeric_limits<T>::epsilon() * scale_factor;
}
