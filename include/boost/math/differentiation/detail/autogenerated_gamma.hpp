#ifndef REVERSE_MODE_AUTODIFF_GENERATED_H
#define REVERSE_MODE_AUTODIFF_GENERATED_H

#include <boost/math/differentiation/detail/reverse_mode_autodiff_basic_operator_overloads.hpp>
#include <boost/math/differentiation/detail/reverse_mode_autodiff_expression_template_base.hpp>
#include <boost/math/differentiation/detail/reverse_mode_autodiff_stl_overloads.hpp>
#include <boost/math/differentiation/detail/reverse_mode_autodiff_utilities.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <string>

namespace boost {
namespace math {
namespace differentiation {
namespace reverse_mode {

// Forward declarations for all expression structs
template<typename T, size_t order, typename ARG>
struct tgamma_expr;

template<typename T, size_t order, typename ARG>
struct lgamma_expr;

template<typename T, size_t order, typename ARG>
struct digamma_expr;

template<typename T, size_t order, typename ARG>
struct trigamma_expr;

template<typename T, size_t order, typename ARG>
struct polygamma_expr;

template<typename T, size_t order, typename LHS, typename RHS>
struct tgamma_ratio_expr;

template<typename T, size_t order, typename ARG>
struct tgamma_ratio_left_float_expr;

template<typename T, size_t order, typename ARG>
struct tgamma_ratio_right_float_expr;

template<typename T, size_t order, typename LHS, typename RHS>
struct tgamma_delta_ratio_expr;

template<typename T, size_t order, typename ARG>
struct tgamma_delta_ratio_left_float_expr;

template<typename T, size_t order, typename ARG>
struct tgamma_delta_ratio_right_float_expr;

template<typename T, size_t order, typename LHS, typename RHS>
struct gamma_p_expr;

template<typename T, size_t order, typename ARG>
struct gamma_p_left_float_expr;

template<typename T, size_t order, typename ARG>
struct gamma_p_right_float_expr;


/*******************************************************************************************/

template<typename T, size_t order, typename ARG>
tgamma_expr<T, order, ARG> tgamma(const expression<T, order, ARG> &arg)
{
    return tgamma_expr<T, order, ARG>(arg, 0.0);
}
template<typename T, size_t order, typename ARG>
lgamma_expr<T, order, ARG> lgamma(const expression<T, order, ARG> &arg)
{
    return lgamma_expr<T, order, ARG>(arg, 0.0);
}
template<typename T, size_t order, typename ARG>
digamma_expr<T, order, ARG> digamma(const expression<T, order, ARG> &arg)
{
    return digamma_expr<T, order, ARG>(arg, 0.0);
}
template<typename T, size_t order, typename ARG>
trigamma_expr<T, order, ARG> trigamma(const expression<T, order, ARG> &arg)
{
    return trigamma_expr<T, order, ARG>(arg, 0.0);
}
template<typename T, size_t order, typename ARG>
polygamma_expr<T, order, ARG> polygamma(const expression<T, order, ARG> &arg)
{
    return polygamma_expr<T, order, ARG>(arg, 0.0);
}
template<typename T, size_t order, typename LHS, typename RHS>
tgamma_ratio_expr<T, order, LHS, RHS> tgamma_ratio(const expression<T, order, LHS> &left,
                                                     const expression<T, order, RHS> &right)
{
    return tgamma_ratio_expr<T, order, LHS, RHS>(left, right);
}
template<typename T, size_t order, typename ARG>
tgamma_ratio_left_float_expr<T, order, ARG> tgamma_ratio(const T                         &left,
                                                          const expression<T, order, ARG> &arg)
{
    return tgamma_ratio_left_float_expr<T, order, ARG>(arg, left);
}

template<typename T, size_t order, typename ARG>
tgamma_ratio_right_float_expr<T, order, ARG> tgamma_ratio(const expression<T, order, ARG> &arg,
                                                           const T                         &right)
{
    return tgamma_ratio_right_float_expr<T, order, ARG>(arg, right);
}
template<typename T, size_t order, typename LHS, typename RHS>
tgamma_delta_ratio_expr<T, order, LHS, RHS> tgamma_delta_ratio(const expression<T, order, LHS> &left,
                                                     const expression<T, order, RHS> &right)
{
    return tgamma_delta_ratio_expr<T, order, LHS, RHS>(left, right);
}
template<typename T, size_t order, typename ARG>
tgamma_delta_ratio_left_float_expr<T, order, ARG> tgamma_delta_ratio(const T                         &left,
                                                          const expression<T, order, ARG> &arg)
{
    return tgamma_delta_ratio_left_float_expr<T, order, ARG>(arg, left);
}

template<typename T, size_t order, typename ARG>
tgamma_delta_ratio_right_float_expr<T, order, ARG> tgamma_delta_ratio(const expression<T, order, ARG> &arg,
                                                           const T                         &right)
{
    return tgamma_delta_ratio_right_float_expr<T, order, ARG>(arg, right);
}
template<typename T, size_t order, typename LHS, typename RHS>
gamma_p_expr<T, order, LHS, RHS> gamma_p(const expression<T, order, LHS> &left,
                                                     const expression<T, order, RHS> &right)
{
    return gamma_p_expr<T, order, LHS, RHS>(left, right);
}
template<typename T, size_t order, typename ARG>
gamma_p_left_float_expr<T, order, ARG> gamma_p(const T                         &left,
                                                          const expression<T, order, ARG> &arg)
{
    return gamma_p_left_float_expr<T, order, ARG>(arg, left);
}

template<typename T, size_t order, typename ARG>
gamma_p_right_float_expr<T, order, ARG> gamma_p(const expression<T, order, ARG> &arg,
                                                           const T                         &right)
{
    return gamma_p_right_float_expr<T, order, ARG>(arg, right);
}

/*******************************************************************************************/

template<typename T, size_t order, typename ARG>
struct tgamma_expr : public abstract_unary_expression<T, order, ARG, tgamma_expr<T, order, ARG>>
{
    /** @brief tgamma(x)
     * d/dx tgamma(x) = gamma(x)*digamma(x) 
      * */
    using arg_type   = ARG;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;

    explicit tgamma_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, tgamma_expr<T, order, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        return detail::if_functional_dispatch<(order > 1)>(
            [this](auto &&x) { return reverse_mode::tgamma(std::forward<decltype(x)>(x)); },
            [this](auto &&x) { return boost::math::tgamma(std::forward<decltype(x)>(x)); },
            this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        return detail::if_functional_dispatch<(
            order > 1)>([](auto &&x) { return reverse_mode::tgamma(x) * reverse_mode::digamma(x); },
                        [](auto &&x) { return boost::math::tgamma(x) * boost::math::digamma(x); },
                        argv);
    }
};

template<typename T, size_t order, typename ARG>
struct lgamma_expr : public abstract_unary_expression<T, order, ARG, lgamma_expr<T, order, ARG>>
{
    /** @brief lgamma(x)
     * d/dx lgamma(x) = digamma(x) 
      * */
    using arg_type   = ARG;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;

    explicit lgamma_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, lgamma_expr<T, order, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        return detail::if_functional_dispatch<(order > 1)>(
            [this](auto &&x) { return reverse_mode::lgamma(std::forward<decltype(x)>(x)); },
            [this](auto &&x) { return boost::math::lgamma(std::forward<decltype(x)>(x)); },
            this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        return detail::if_functional_dispatch<(
            order > 1)>([](auto &&x) { return reverse_mode::digamma(x); },
                        [](auto &&x) { return boost::math::digamma(x); },
                        argv);
    }
};

template<typename T, size_t order, typename ARG>
struct digamma_expr : public abstract_unary_expression<T, order, ARG, digamma_expr<T, order, ARG>>
{
    /** @brief digamma(x)
     * d/dx digamma(x) = polygamma(1,x) 
      * */
    using arg_type   = ARG;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;

    explicit digamma_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, digamma_expr<T, order, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        return detail::if_functional_dispatch<(order > 1)>(
            [this](auto &&x) { return reverse_mode::digamma(std::forward<decltype(x)>(x)); },
            [this](auto &&x) { return boost::math::digamma(std::forward<decltype(x)>(x)); },
            this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        return detail::if_functional_dispatch<(
            order > 1)>([](auto &&x) { return reverse_mode::polygamma(1, x); },
                        [](auto &&x) { return boost::math::polygamma(1, x); },
                        argv);
    }
};

template<typename T, size_t order, typename ARG>
struct trigamma_expr : public abstract_unary_expression<T, order, ARG, trigamma_expr<T, order, ARG>>
{
    /** @brief trigamma(x)
     * d/dx trigamma(x) = polygamma(2,x) 
      * */
    using arg_type   = ARG;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;

    explicit trigamma_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, trigamma_expr<T, order, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        return detail::if_functional_dispatch<(order > 1)>(
            [this](auto &&x) { return reverse_mode::trigamma(std::forward<decltype(x)>(x)); },
            [this](auto &&x) { return boost::math::trigamma(std::forward<decltype(x)>(x)); },
            this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        return detail::if_functional_dispatch<(
            order > 1)>([](auto &&x) { return reverse_mode::polygamma(2, x); },
                        [](auto &&x) { return boost::math::polygamma(2, x); },
                        argv);
    }
};

template<typename T, size_t order, typename ARG>
struct polygamma_expr : public abstract_unary_expression<T, order, ARG, polygamma_expr<T, order, ARG>>
{
    /** @brief polygamma(x)
     * d/dx polygamma(x) = polygamma(n+1,x) 
      * */
    using arg_type   = ARG;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;

    explicit polygamma_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, polygamma_expr<T, order, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        return detail::if_functional_dispatch<(order > 1)>(
            [this](auto &&x) { return reverse_mode::polygamma(std::forward<decltype(x)>(x)); },
            [this](auto &&x) { return boost::math::polygamma(std::forward<decltype(x)>(x)); },
            this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        return detail::if_functional_dispatch<(
            order > 1)>([constant](auto &&x) { return reverse_mode::polygamma(static_cast<int>(constant) + 1, x); },
                        [constant](auto &&x) { return boost::math::polygamma(static_cast<int>(constant) + 1, x); },
                        argv);
    }
};

template<typename T, size_t order, typename LHS, typename RHS>
struct tgamma_ratio_expr
    : public abstract_binary_expression<T, order, LHS, RHS, tgamma_ratio_expr<T, order, LHS, RHS>>
{
    using lhs_type   = LHS;
    using rhs_type   = RHS;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;
    // Explicitly define constructor to forward to base class
    explicit tgamma_ratio_expr(const expression<T, order, LHS> &left_hand_expr,
                               const expression<T, order, RHS> &right_hand_expr)
        : abstract_binary_expression<T, order, LHS, RHS, tgamma_ratio_expr<T, order, LHS, RHS>>(
              left_hand_expr, right_hand_expr)
    {}

    inner_t evaluate() const
    {
        return detail::if_functional_dispatch<(order > 1)>(
            [this](auto &&x, auto &&y) {
                return reverse_mode::tgamma_ratio(std::forward<decltype(x)>(x),
                                                  std::forward<decltype(y)>(y));
            },
            [this](auto &&x, auto &&y) {
                return boost::math::tgamma_ratio(std::forward<decltype(x)>(x),
                                                 std::forward<decltype(y)>(y));
            },
            this->lhs.evaluate(),
            this->rhs.evaluate());
    };
    static const inner_t left_derivative(const inner_t &l, const inner_t &r, const inner_t &v)
    {
        return detail::if_functional_dispatch<(order > 1)>(
            [](auto &&x, auto &&y) {
                return reverse_mode::tgamma_ratio(x, y) * reverse_mode::digamma(x);
            },
            [](auto &&x, auto &&y) {
                return boost::math::tgamma_ratio(x, y) * boost::math::digamma(x);
            },
            l,
            r);
    };
    static const inner_t right_derivative(const inner_t &l, const inner_t &r, const inner_t &v)
    {
        return detail::if_functional_dispatch<(order > 1)>(
            [](auto &&x, auto &&y) {
                return -reverse_mode::tgamma_ratio(x, y) * reverse_mode::digamma(y);
            },
            [](auto &&x, auto &&y) {
                return -boost::math::tgamma_ratio(x, y) * boost::math::digamma(y);
            },
            l,
            r);
    };
};

template<typename T, size_t order, typename ARG>
struct tgamma_ratio_left_float_expr
    : public abstract_unary_expression<T, order, ARG, tgamma_ratio_left_float_expr<T, order, ARG>>
{
    using arg_type   = ARG;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;

    explicit tgamma_ratio_left_float_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, tgamma_ratio_left_float_expr<T, order, ARG>>(
              arg_expr, v){};

    inner_t evaluate() const
    {
        return detail::if_functional_dispatch<(order > 1)>(
            [this](auto &&x) {
                return reverse_mode::tgamma_ratio(this->constant, std::forward<decltype(x)>(x));
            },
            [this](auto &&x) {
                return boost::math::tgamma_ratio(this->constant, std::forward<decltype(x)>(x));
            },
            this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        return detail::if_functional_dispatch<(order > 1)>(
            [constant](auto &&x) {
                return -reverse_mode::tgamma_ratio(x, constant) * reverse_mode::digamma(constant);
            },
            [constant](auto &&x) {
                return -boost::math::tgamma_ratio(x, constant) * boost::math::digamma(constant);
            },
            argv);
    }
};

template<typename T, size_t order, typename ARG>
struct tgamma_ratio_right_float_expr
    : public abstract_unary_expression<T, order, ARG, tgamma_ratio_right_float_expr<T, order, ARG>>
{
    using arg_type   = ARG;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;

    explicit tgamma_ratio_right_float_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, tgamma_ratio_right_float_expr<T, order, ARG>>(
              arg_expr, v){};

    inner_t evaluate() const
    {
        return detail::if_functional_dispatch<(order > 1)>(
            [this](auto &&x) {
                return reverse_mode::tgamma_ratio(std::forward<decltype(x)>(x), this->constant);
            },
            [this](auto &&x) {
                return boost::math::tgamma_ratio(std::forward<decltype(x)>(x), this->constant);
            },
            this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        return detail::if_functional_dispatch<(order > 1)>(
            [constant](auto &&x) {
                return reverse_mode::tgamma_ratio(x, constant) * reverse_mode::digamma(x);
            },
            [constant](auto &&x) {
                return boost::math::tgamma_ratio(x, constant) * boost::math::digamma(x);
            },
            argv);
    }
};

template<typename T, size_t order, typename LHS, typename RHS>
struct tgamma_delta_ratio_expr
    : public abstract_binary_expression<T, order, LHS, RHS, tgamma_delta_ratio_expr<T, order, LHS, RHS>>
{
    using lhs_type   = LHS;
    using rhs_type   = RHS;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;
    // Explicitly define constructor to forward to base class
    explicit tgamma_delta_ratio_expr(const expression<T, order, LHS> &left_hand_expr,
                               const expression<T, order, RHS> &right_hand_expr)
        : abstract_binary_expression<T, order, LHS, RHS, tgamma_delta_ratio_expr<T, order, LHS, RHS>>(
              left_hand_expr, right_hand_expr)
    {}

    inner_t evaluate() const
    {
        return detail::if_functional_dispatch<(order > 1)>(
            [this](auto &&x, auto &&y) {
                return reverse_mode::tgamma_delta_ratio(std::forward<decltype(x)>(x),
                                                  std::forward<decltype(y)>(y));
            },
            [this](auto &&x, auto &&y) {
                return boost::math::tgamma_delta_ratio(std::forward<decltype(x)>(x),
                                                 std::forward<decltype(y)>(y));
            },
            this->lhs.evaluate(),
            this->rhs.evaluate());
    };
    static const inner_t left_derivative(const inner_t &l, const inner_t &r, const inner_t &v)
    {
        return detail::if_functional_dispatch<(order > 1)>(
            [](auto &&x, auto &&y) {
                return reverse_mode::tgamma_delta_ratio(x, y) * (reverse_mode::digamma(x) - reverse_mode::digamma(x + y));
            },
            [](auto &&x, auto &&y) {
                return boost::math::tgamma_delta_ratio(x, y) * (boost::math::digamma(x) - boost::math::digamma(x + y));
            },
            l,
            r);
    };
    static const inner_t right_derivative(const inner_t &l, const inner_t &r, const inner_t &v)
    {
        return detail::if_functional_dispatch<(order > 1)>(
            [](auto &&x, auto &&y) {
                return -reverse_mode::tgamma_delta_ratio(x, y) * (reverse_mode::digamma(x + y));
            },
            [](auto &&x, auto &&y) {
                return -boost::math::tgamma_delta_ratio(x, y) * (boost::math::digamma(x + y));
            },
            l,
            r);
    };
};

template<typename T, size_t order, typename ARG>
struct tgamma_delta_ratio_left_float_expr
    : public abstract_unary_expression<T, order, ARG, tgamma_delta_ratio_left_float_expr<T, order, ARG>>
{
    using arg_type   = ARG;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;

    explicit tgamma_delta_ratio_left_float_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, tgamma_delta_ratio_left_float_expr<T, order, ARG>>(
              arg_expr, v){};

    inner_t evaluate() const
    {
        return detail::if_functional_dispatch<(order > 1)>(
            [this](auto &&x) {
                return reverse_mode::tgamma_delta_ratio(this->constant, std::forward<decltype(x)>(x));
            },
            [this](auto &&x) {
                return boost::math::tgamma_delta_ratio(this->constant, std::forward<decltype(x)>(x));
            },
            this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        return detail::if_functional_dispatch<(order > 1)>(
            [constant](auto &&x) {
                return -reverse_mode::tgamma_delta_ratio(x, constant) * (reverse_mode::digamma(x + constant));
            },
            [constant](auto &&x) {
                return -boost::math::tgamma_delta_ratio(x, constant) * (boost::math::digamma(x + constant));
            },
            argv);
    }
};

template<typename T, size_t order, typename ARG>
struct tgamma_delta_ratio_right_float_expr
    : public abstract_unary_expression<T, order, ARG, tgamma_delta_ratio_right_float_expr<T, order, ARG>>
{
    using arg_type   = ARG;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;

    explicit tgamma_delta_ratio_right_float_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, tgamma_delta_ratio_right_float_expr<T, order, ARG>>(
              arg_expr, v){};

    inner_t evaluate() const
    {
        return detail::if_functional_dispatch<(order > 1)>(
            [this](auto &&x) {
                return reverse_mode::tgamma_delta_ratio(std::forward<decltype(x)>(x), this->constant);
            },
            [this](auto &&x) {
                return boost::math::tgamma_delta_ratio(std::forward<decltype(x)>(x), this->constant);
            },
            this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        return detail::if_functional_dispatch<(order > 1)>(
            [constant](auto &&x) {
                return reverse_mode::tgamma_delta_ratio(x, constant) * (reverse_mode::digamma(x) - reverse_mode::digamma(x + constant));
            },
            [constant](auto &&x) {
                return boost::math::tgamma_delta_ratio(x, constant) * (boost::math::digamma(x) - boost::math::digamma(x + constant));
            },
            argv);
    }
};

template<typename T, size_t order, typename LHS, typename RHS>
struct gamma_p_expr
    : public abstract_binary_expression<T, order, LHS, RHS, gamma_p_expr<T, order, LHS, RHS>>
{
    using lhs_type   = LHS;
    using rhs_type   = RHS;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;
    // Explicitly define constructor to forward to base class
    explicit gamma_p_expr(const expression<T, order, LHS> &left_hand_expr,
                               const expression<T, order, RHS> &right_hand_expr)
        : abstract_binary_expression<T, order, LHS, RHS, gamma_p_expr<T, order, LHS, RHS>>(
              left_hand_expr, right_hand_expr)
    {}

    inner_t evaluate() const
    {
        return detail::if_functional_dispatch<(order > 1)>(
            [this](auto &&x, auto &&y) {
                return reverse_mode::gamma_p(std::forward<decltype(x)>(x),
                                                  std::forward<decltype(y)>(y));
            },
            [this](auto &&x, auto &&y) {
                return boost::math::gamma_p(std::forward<decltype(x)>(x),
                                                 std::forward<decltype(y)>(y));
            },
            this->lhs.evaluate(),
            this->rhs.evaluate());
    };
    static const inner_t left_derivative(const inner_t &l, const inner_t &r, const inner_t &v)
    {
        return detail::if_functional_dispatch<(order > 1)>(
            [](auto &&x, auto &&y) {
                return reverse_mode::tgamma_delta_ratio(x, y) * (reverse_mode::digamma(x) - reverse_mode::digamma(x + y));
            },
            [](auto &&x, auto &&y) {
                return boost::math::tgamma_delta_ratio(x, y) * (boost::math::digamma(x) - boost::math::digamma(x + y));
            },
            l,
            r);
    };
    static const inner_t right_derivative(const inner_t &l, const inner_t &r, const inner_t &v)
    {
        return detail::if_functional_dispatch<(order > 1)>(
            [](auto &&x, auto &&y) {
                return -reverse_mode::tgamma_delta_ratio(x, y) * (reverse_mode::digamma(x + y));
            },
            [](auto &&x, auto &&y) {
                return -boost::math::tgamma_delta_ratio(x, y) * (boost::math::digamma(x + y));
            },
            l,
            r);
    };
};

template<typename T, size_t order, typename ARG>
struct gamma_p_left_float_expr
    : public abstract_unary_expression<T, order, ARG, gamma_p_left_float_expr<T, order, ARG>>
{
    using arg_type   = ARG;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;

    explicit gamma_p_left_float_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, gamma_p_left_float_expr<T, order, ARG>>(
              arg_expr, v){};

    inner_t evaluate() const
    {
        return detail::if_functional_dispatch<(order > 1)>(
            [this](auto &&x) {
                return reverse_mode::gamma_p(this->constant, std::forward<decltype(x)>(x));
            },
            [this](auto &&x) {
                return boost::math::gamma_p(this->constant, std::forward<decltype(x)>(x));
            },
            this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        return detail::if_functional_dispatch<(order > 1)>(
            [constant](auto &&x) {
                return -reverse_mode::tgamma_delta_ratio(x, constant) * (reverse_mode::digamma(x + constant));
            },
            [constant](auto &&x) {
                return -boost::math::tgamma_delta_ratio(x, constant) * (boost::math::digamma(x + constant));
            },
            argv);
    }
};

template<typename T, size_t order, typename ARG>
struct gamma_p_right_float_expr
    : public abstract_unary_expression<T, order, ARG, gamma_p_right_float_expr<T, order, ARG>>
{
    using arg_type   = ARG;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;

    explicit gamma_p_right_float_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, gamma_p_right_float_expr<T, order, ARG>>(
              arg_expr, v){};

    inner_t evaluate() const
    {
        return detail::if_functional_dispatch<(order > 1)>(
            [this](auto &&x) {
                return reverse_mode::gamma_p(std::forward<decltype(x)>(x), this->constant);
            },
            [this](auto &&x) {
                return boost::math::gamma_p(std::forward<decltype(x)>(x), this->constant);
            },
            this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        return detail::if_functional_dispatch<(order > 1)>(
            [constant](auto &&x) {
                return reverse_mode::tgamma_delta_ratio(x, constant) * (reverse_mode::digamma(x) - reverse_mode::digamma(x + constant));
            },
            [constant](auto &&x) {
                return boost::math::tgamma_delta_ratio(x, constant) * (boost::math::digamma(x) - boost::math::digamma(x + constant));
            },
            argv);
    }
};

} // namespace reverse_mode
} // namespace differentiation
} // namespace math
} // namespace boost

#endif

