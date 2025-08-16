#ifndef REVERSE_MODE_AUTODIFF_ERF_OVERLOADS_HPP
#define REVERSE_MODE_AUTODIFF_ERF_OVERLOADS_HPP

#include <boost/math/constants/constants.hpp>
#include <boost/math/differentiation/detail/reverse_mode_autodiff_expression_template_base.hpp>
#include <boost/math/differentiation/detail/reverse_mode_autodiff_stl_overloads.hpp>
#include <boost/math/differentiation/detail/reverse_mode_autodiff_utilities.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
namespace boost {
namespace math {
namespace differentiation {
namespace reverse_mode {

template<typename T, size_t order, typename ARG>
struct erf_expr;

template<typename T, size_t order, typename ARG>
struct erfc_expr;

template<typename T, size_t order, typename ARG>
struct erf_inv_expr;

template<typename T, size_t order, typename ARG>
struct erfc_inv_expr;

template<typename T, size_t order, typename ARG>
erf_expr<T, order, ARG> erf(const expression<T, order, ARG> &arg)
{
    return erf_expr<T, order, ARG>(arg, 0.0);
}

template<typename T, size_t order, typename ARG>
erfc_expr<T, order, ARG> erfc(const expression<T, order, ARG> &arg)
{
    return erfc_expr<T, order, ARG>(arg, 0.0);
}

template<typename T, size_t order, typename ARG>
erf_inv_expr<T, order, ARG> erf_inv(const expression<T, order, ARG> &arg)
{
    return erf_inv_expr<T, order, ARG>(arg, 0.0);
}

template<typename T, size_t order, typename ARG>
erfc_inv_expr<T, order, ARG> erfc_inv(const expression<T, order, ARG> &arg)
{
    return erfc_inv_expr<T, order, ARG>(arg, 0.0);
}

template<typename T, size_t order, typename ARG>
struct erf_expr : public abstract_unary_expression<T, order, ARG, erf_expr<T, order, ARG>>
{
    /** @brief erf(x)
    *
    * d/dx erf(x) = 2*exp(x^2)/sqrt(pi)
    *
    * */
    using arg_type   = ARG;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;

    explicit erf_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, erf_expr<T, order, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        return detail::if_functional_dispatch<(order > 1)>(
            [this](auto &&x) { return reverse_mode::erf(std::forward<decltype(x)>(x)); },
            [this](auto &&x) { return boost::math::erf(std::forward<decltype(x)>(x)); },
            this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        using multiprecision::sqrt;
        using std::exp;
        using std::sqrt;
        return 2 * exp(-argv * argv) / sqrt(constants::pi<T>());
    }
};

template<typename T, size_t order, typename ARG>
struct erfc_expr : public abstract_unary_expression<T, order, ARG, erfc_expr<T, order, ARG>>
{
    /** @brief erfc(x)
    *
    * d/dx erf(x) = -2*exp(x^2)/sqrt(pi)
    *
    * */
    using arg_type   = ARG;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;

    explicit erfc_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, erfc_expr<T, order, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        return detail::if_functional_dispatch<(order > 1)>(
            [this](auto &&x) { return reverse_mode::erfc(std::forward<decltype(x)>(x)); },
            [this](auto &&x) { return boost::math::erfc(std::forward<decltype(x)>(x)); },
            this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        using multiprecision::sqrt;
        using std::exp;
        using std::sqrt;
        return -2 * exp(-argv * argv) / sqrt(constants::pi<T>());
    }
};

template<typename T, size_t order, typename ARG>
struct erf_inv_expr : public abstract_unary_expression<T, order, ARG, erf_inv_expr<T, order, ARG>>
{
    /** @brief erf(x)
    *
    * d/dx erf(x) = 2*exp(x^2)/sqrt(pi)
    *
    * */
    using arg_type   = ARG;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;

    explicit erf_inv_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, erf_inv_expr<T, order, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        return detail::if_functional_dispatch<(order > 1)>(
            [this](auto &&x) { return reverse_mode::erf_inv(std::forward<decltype(x)>(x)); },
            [this](auto &&x) { return boost::math::erf_inv(std::forward<decltype(x)>(x)); },
            this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        using multiprecision::exp;
        using multiprecision::pow;
        using multiprecision::sqrt;
        using std::exp;
        using std::pow;
        using std::sqrt;
        return detail::if_functional_dispatch<(order > 1)>(
            [](auto &&x) {
                return 0.5 * sqrt(constants::pi<T>())
                       * reverse_mode::exp(reverse_mode::pow(reverse_mode::erf_inv(x), 2.0));
            },
            [](auto &&x) {
                return 0.5 * sqrt(constants::pi<T>()) * exp(pow(boost::math::erf_inv(x), 2));
            },
            argv);
    }
};

template<typename T, size_t order, typename ARG>
struct erfc_inv_expr : public abstract_unary_expression<T, order, ARG, erfc_inv_expr<T, order, ARG>>
{
    /** @brief erfc(x)
    *
    * d/dx erf(x) = -2*exp(x^2)/sqrt(pi)
    *
    * */
    using arg_type   = ARG;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;

    explicit erfc_inv_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, erfc_inv_expr<T, order, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        return detail::if_functional_dispatch<(order > 1)>(
            [this](auto &&x) { return reverse_mode::erfc_inv(std::forward<decltype(x)>(x)); },
            [this](auto &&x) { return boost::math::erfc_inv(std::forward<decltype(x)>(x)); },
            this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        using multiprecision::exp;
        using multiprecision::pow;
        using multiprecision::sqrt;
        using std::exp;
        using std::pow;
        using std::sqrt;
        return detail::if_functional_dispatch<(order > 1)>(
            [](auto &&x) {
                return -0.5 * sqrt(constants::pi<T>())
                       * reverse_mode::exp(reverse_mode::pow(reverse_mode::erfc_inv(x), 2.0));
            },
            [](auto &&x) {
                return -0.5 * sqrt(constants::pi<T>()) * exp(pow(boost::math::erfc_inv(x), 2));
            },
            argv);
    }
};

} // namespace reverse_mode
} // namespace differentiation
} // namespace math
} // namespace boost

#endif // REVERSE_MODE_AUTODIFF_ERF_OVERLOADS_HPP
