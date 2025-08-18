//           Copyright Maksym Zhelyenzyakov 2025-2026.
// Distributed under the Boost Software License, Version 1.0.
//      (See accompanying file LICENSE_1_0.txt or copy at
//           https://www.boost.org/LICENSE_1_0.txt)
#ifndef REVERSE_MODE_AUTODIFF_STL_OVERLOADS
#define REVERSE_MODE_AUTODIFF_STL_OVERLOADS
/* stl support : expressions */
#include <boost/math/differentiation/detail/reverse_mode_autodiff_basic_operator_overloads.hpp>
#include <boost/math/differentiation/detail/reverse_mode_autodiff_expression_template_base.hpp>
#include <boost/math/special_functions/round.hpp>
#include <boost/math/special_functions/trunc.hpp>
#include <cmath>
#include <complex>
namespace boost {
namespace math {
namespace differentiation {
namespace reverse_mode {
template<typename T, size_t order, typename ARG>
struct fabs_expr : public abstract_unary_expression<T, order, ARG, fabs_expr<T, order, ARG>>
{
    /** @brief
    * |x|
    * d/dx |x| = 1 if x > 0
    *          -1 if x <= 0
    *
    * the choice is arbitrary and for optimization it is most likely
    * more correct to chose this convention over d/dx = 0  at x = 0
    * to avoid vanishing gradients
    * */
    using inner_t    = rvar_t<T, order - 1>;

    explicit fabs_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, fabs_expr<T, order, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return fabs(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        return argv > 0.0 ? inner_t{1.0} : inner_t{-1.0};
    }
};

template<typename T, size_t order, typename ARG>
struct ceil_expr : public abstract_unary_expression<T, order, ARG, ceil_expr<T, order, ARG>>
{
    /** @brief ceil(1.11) = 2.0
    *
    * d/dx ceil(x) = 0.0 for all x
    *
    * we avoid problematic points at x = 1,2,3...
    * as with optimization its most likely intented
    * this function's derivative is 0.0;
    * */
    using inner_t    = rvar_t<T, order - 1>;

    explicit ceil_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, ceil_expr<T, order, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return ceil(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        return inner_t{0.0};
    }
};

template<typename T, size_t order, typename ARG>
struct floor_expr : public abstract_unary_expression<T, order, ARG, floor_expr<T, order, ARG>>
{
    /** @brief floor(1.11) = 1.0, floor(-1.11) = 2
    *
    * d/dx floor(x) = 0.0 for all x
    *
    * we avoid problematic points at x = 1,2,3...
    * as with optimization its most likely intented
    * this function's derivative is 0.0;
    * */
    using inner_t    = rvar_t<T, order - 1>;

    explicit floor_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, floor_expr<T, order, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return floor(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        return inner_t{0.0};
    }
};

template<typename T, size_t order, typename ARG>
struct trunc_expr : public abstract_unary_expression<T, order, ARG, trunc_expr<T, order, ARG>>
{
    /** @brief trunc(1.11) = 1.0, trunc(-1.11) = -1.0
    *
    * d/dx trunc(x) = 0.0 for all x
    *
    * we avoid problematic points at x = 1,2,3...
    * as with optimization its most likely intented
    * this function's derivative is 0.0;
    * */
    using inner_t    = rvar_t<T, order - 1>;

    explicit trunc_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, trunc_expr<T, order, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return trunc(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        return inner_t{0.0};
    }
};

template<typename T, size_t order, typename ARG>
struct exp_expr : public abstract_unary_expression<T, order, ARG, exp_expr<T, order, ARG>>
{
    /** @brief exp(x)
    *
    * d/dx exp(x) = exp(x)
    *
    * */
    using inner_t    = rvar_t<T, order - 1>;

    explicit exp_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, exp_expr<T, order, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return exp(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        BOOST_MATH_STD_USING
        return exp(argv);
    }
};

template<typename T, size_t order, typename LHS, typename RHS>
struct pow_expr
    : public abstract_binary_expression<T, order, LHS, RHS, pow_expr<T, order, LHS, RHS>>
{
    /** @brief pow(x,y)
     *  d/dx pow(x,y) = y pow (x, y-1)
     *  d/dy pow(x,y) = pow(x,y) log(x)
    * */
    using inner_t    = rvar_t<T, order - 1>;
    // Explicitly define constructor to forward to base class
    explicit pow_expr(const expression<T, order, LHS> &left_hand_expr,
                      const expression<T, order, RHS> &right_hand_expr)
        : abstract_binary_expression<T, order, LHS, RHS, pow_expr<T, order, LHS, RHS>>(
              left_hand_expr, right_hand_expr)
    {}

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return pow(this->lhs.evaluate(), this->rhs.evaluate());
    };
    static const inner_t left_derivative(const inner_t &l, const inner_t &r, const inner_t &v)
    {
        BOOST_MATH_STD_USING
        return r * pow(l, r - static_cast<T>(1.0));
    };
    static const inner_t right_derivative(const inner_t &l, const inner_t &r, const inner_t &v)
    {
        BOOST_MATH_STD_USING
        return pow(l, r) * log(l);
    };
};

template<typename T, size_t order, typename ARG>
struct expr_pow_float_expr
    : public abstract_unary_expression<T, order, ARG, expr_pow_float_expr<T, order, ARG>>
{
    /** @brief pow(rvar,float)
      */
    using inner_t    = rvar_t<T, order - 1>;

    explicit expr_pow_float_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, expr_pow_float_expr<T, order, ARG>>(arg_expr,
                                                                                       v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return pow(this->arg.evaluate(), this->constant);
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        BOOST_MATH_STD_USING
        return inner_t{constant} * pow(argv, inner_t{constant - 1});
    }
};

template<typename T, size_t order, typename ARG>
struct float_pow_expr_expr
    : public abstract_unary_expression<T, order, ARG, float_pow_expr_expr<T, order, ARG>>
{
    /** @brief pow(float, rvar)
      * */
    using inner_t    = rvar_t<T, order - 1>;

    explicit float_pow_expr_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, float_pow_expr_expr<T, order, ARG>>(arg_expr,
                                                                                       v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return pow(this->constant, this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        BOOST_MATH_STD_USING
        return pow(constant, argv) * log(constant);
    }
};

template<typename T, size_t order, typename ARG>
struct sqrt_expr : public abstract_unary_expression<T, order, ARG, sqrt_expr<T, order, ARG>>
{
    /** @brief  sqrt(x)
     *  d/dx sqrt(x) = 1/(2 sqrt(x))
    * */
    using inner_t    = rvar_t<T, order - 1>;

    explicit sqrt_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, sqrt_expr<T, order, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return sqrt(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        BOOST_MATH_STD_USING
        return static_cast<T>(1.0) / (static_cast<T>(2.0) * sqrt(argv));
    }
};

template<typename T, size_t order, typename ARG>
struct log_expr : public abstract_unary_expression<T, order, ARG, log_expr<T, order, ARG>>
{
    /** @brief log(x)
     *  d/dx log(x) = 1/x
      * */
    using inner_t    = rvar_t<T, order - 1>;

    explicit log_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, log_expr<T, order, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return log(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        return static_cast<T>(1.0) / argv;
    }
};

template<typename T, size_t order, typename ARG>
struct cos_expr : public abstract_unary_expression<T, order, ARG, cos_expr<T, order, ARG>>
{
    /** @brief cos(x)
     *  d/dx cos(x) = -sin(x)
      * */
    using inner_t    = rvar_t<T, order - 1>;

    explicit cos_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, cos_expr<T, order, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return cos(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        BOOST_MATH_STD_USING
        return -sin(argv);
    }
};

template<typename T, size_t order, typename ARG>
struct sin_expr : public abstract_unary_expression<T, order, ARG, sin_expr<T, order, ARG>>
{
    /** @brief sin(x)
     *  d/dx sin(x) = cos(x)
      * */
    using arg_type   = ARG;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;

    explicit sin_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, sin_expr<T, order, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return sin(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        BOOST_MATH_STD_USING
        return cos(argv);
    }
};

template<typename T, size_t order, typename ARG>
struct tan_expr : public abstract_unary_expression<T, order, ARG, tan_expr<T, order, ARG>>
{
    /** @brief tan(x)
     *  d/dx tan(x) = 1/cos^2(x)
      * */
    using inner_t    = rvar_t<T, order - 1>;

    explicit tan_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, tan_expr<T, order, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return tan(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        BOOST_MATH_STD_USING
        return static_cast<T>(1.0) / (cos(argv) * cos(argv));
    }
};

template<typename T, size_t order, typename ARG>
struct acos_expr : public abstract_unary_expression<T, order, ARG, acos_expr<T, order, ARG>>
{
    /** @brief acos(x)
     *  d/dx acos(x) = -1/sqrt(1-x^2)
      * */
    using inner_t    = rvar_t<T, order - 1>;

    explicit acos_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, acos_expr<T, order, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return acos(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        BOOST_MATH_STD_USING
        return static_cast<T>(-1.0) / sqrt(static_cast<T>(1.0) - argv * argv);
    }
};

template<typename T, size_t order, typename ARG>
struct asin_expr : public abstract_unary_expression<T, order, ARG, asin_expr<T, order, ARG>>
{
    /** @brief asin(x)
     *  d/dx asin =  1/sqrt(1-x^2)
      * */
    using arg_type   = ARG;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;

    explicit asin_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, asin_expr<T, order, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return asin(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        BOOST_MATH_STD_USING
        return static_cast<T>(1.0) / sqrt(static_cast<T>(1.0) - argv * argv);
    }
};

template<typename T, size_t order, typename ARG>
struct atan_expr : public abstract_unary_expression<T, order, ARG, atan_expr<T, order, ARG>>
{
    /** @brief atan(x)
     *  d/dx atan(x) = 1/x^2+1
      * */
    using inner_t    = rvar_t<T, order - 1>;

    explicit atan_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, atan_expr<T, order, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return atan(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        BOOST_MATH_STD_USING
        return static_cast<T>(1.0) / (static_cast<T>(1.0) + argv * argv);
    }
};
template<typename T, size_t order, typename LHS, typename RHS>
struct atan2_expr
    : public abstract_binary_expression<T, order, LHS, RHS, atan2_expr<T, order, LHS, RHS>>
{
    /** @brief atan2(x,y)
    * */
    using inner_t    = rvar_t<T, order - 1>;
    // Explicitly define constructor to forward to base class
    explicit atan2_expr(const expression<T, order, LHS> &left_hand_expr,
                        const expression<T, order, RHS> &right_hand_expr)
        : abstract_binary_expression<T, order, LHS, RHS, atan2_expr<T, order, LHS, RHS>>(
              left_hand_expr, right_hand_expr)
    {}

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return atan2(this->lhs.evaluate(), this->rhs.evaluate());
    };
    static const inner_t left_derivative(const inner_t &l, const inner_t &r, const inner_t &v)
    {
        return r / (l * l + r * r);
    };
    static const inner_t right_derivative(const inner_t &l, const inner_t &r, const inner_t &v)
    {
        return -l / (l * l + r * r);
    };
};

template<typename T, size_t order, typename ARG>
struct atan2_left_float_expr
    : public abstract_unary_expression<T, order, ARG, atan2_left_float_expr<T, order, ARG>>
{
    /** @brief atan2(float,rvar) 
      * */
    using inner_t    = rvar_t<T, order - 1>;

    explicit atan2_left_float_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, atan2_left_float_expr<T, order, ARG>>(arg_expr,
                                                                                         v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return atan2(this->constant, this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        return -constant / (constant * constant + argv * argv);
    }
};

template<typename T, size_t order, typename ARG>
struct atan2_right_float_expr
    : public abstract_unary_expression<T, order, ARG, atan2_right_float_expr<T, order, ARG>>
{
    /** @brief atan2(rvar,float) 
      * */
    using inner_t    = rvar_t<T, order - 1>;

    explicit atan2_right_float_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, atan2_right_float_expr<T, order, ARG>>(arg_expr,
                                                                                          v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return atan2(this->arg.evaluate(), this->constant);
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        return constant / (constant * constant + argv * argv);
    }
};

template<typename T, size_t order, typename ARG>
struct round_expr : public abstract_unary_expression<T, order, ARG, round_expr<T, order, ARG>>
{
    /** @brief round(x)
     *  d/dx round = 0
      * */
    using inner_t    = rvar_t<T, order - 1>;

    explicit round_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, round_expr<T, order, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return round(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        return inner_t{0.0};
    }
};

template<typename T, size_t order, typename ARG>
struct sinh_expr : public abstract_unary_expression<T, order, ARG, sinh_expr<T, order, ARG>>
{
    /** @brief sinh(x)
     *  d/dx sinh(x) = cosh
      * */
    using inner_t    = rvar_t<T, order - 1>;

    explicit sinh_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, sinh_expr<T, order, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return sinh(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        BOOST_MATH_STD_USING
        return cosh(argv);
    }
};

template<typename T, size_t order, typename ARG>
struct cosh_expr : public abstract_unary_expression<T, order, ARG, cosh_expr<T, order, ARG>>
{
    /** @brief cosh(x)
     *  d/dx cosh(x) = sinh
      * */
    using inner_t    = rvar_t<T, order - 1>;

    explicit cosh_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, cosh_expr<T, order, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return cosh(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        BOOST_MATH_STD_USING
        return sinh(argv);
    }
};
template<typename T, size_t order, typename ARG>
struct tanh_expr : public abstract_unary_expression<T, order, ARG, tanh_expr<T, order, ARG>>
{
    /** @brief tanh(x)
     *  d/dx tanh(x) = 1/cosh^2
      * */
    using inner_t    = rvar_t<T, order - 1>;

    explicit tanh_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, tanh_expr<T, order, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return tanh(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        BOOST_MATH_STD_USING
        return static_cast<T>(1.0) / (cosh(argv) * cosh(argv));
    }
};

template<typename T, size_t order, typename ARG>
struct log10_expr : public abstract_unary_expression<T, order, ARG, log10_expr<T, order, ARG>>
{
    /** @brief log10(x)
     *  d/dx log10(x) = 1/(x * log(10))
      * */
    using inner_t    = rvar_t<T, order - 1>;

    explicit log10_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, log10_expr<T, order, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return log10(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        BOOST_MATH_STD_USING
        return static_cast<T>(1.0) / (argv * log(static_cast<T>(10.0)));
    }
};

template<typename T, size_t order, typename ARG>
struct acosh_expr : public abstract_unary_expression<T, order, ARG, acosh_expr<T, order, ARG>>
{
    /** @brief acosh(x)
     *  d/dx acosh(x) = 1/(sqrt(x-1)sqrt(x+1)
      * */
    using inner_t    = rvar_t<T, order - 1>;

    explicit acosh_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, acosh_expr<T, order, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return acosh(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        BOOST_MATH_STD_USING
        return static_cast<T>(1.0)
               / (sqrt(argv - static_cast<T>(1.0)) * sqrt(argv + static_cast<T>(1.0)));
    }
};

template<typename T, size_t order, typename ARG>
struct asinh_expr : public abstract_unary_expression<T, order, ARG, asinh_expr<T, order, ARG>>
{
    /** @brief asinh(x)
     *  d/dx asinh(x) = 1/(sqrt(1+x^2))
      * */
    using inner_t    = rvar_t<T, order - 1>;

    explicit asinh_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, asinh_expr<T, order, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return asinh(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        BOOST_MATH_STD_USING
        return static_cast<T>(1.0) / (sqrt(static_cast<T>(1.0) + argv * argv));
    }
};

template<typename T, size_t order, typename ARG>
struct atanh_expr : public abstract_unary_expression<T, order, ARG, atanh_expr<T, order, ARG>>
{
    /** @brief atanh(x)
     *  d/dx atanh(x) = 1/(1-x^2)
      * */
    using inner_t    = rvar_t<T, order - 1>;

    explicit atanh_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, atanh_expr<T, order, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return atanh(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        BOOST_MATH_STD_USING
        return static_cast<T>(1.0) / (static_cast<T>(1.0) - argv * argv);
    }
};
template<typename T, size_t order, typename LHS, typename RHS>
struct fmod_expr
    : public abstract_binary_expression<T, order, LHS, RHS, fmod_expr<T, order, LHS, RHS>>
{
    /** @brief 
    * */
    using inner_t    = rvar_t<T, order - 1>;
    // Explicitly define constructor to forward to base class
    explicit fmod_expr(const expression<T, order, LHS> &left_hand_expr,
                       const expression<T, order, RHS> &right_hand_expr)
        : abstract_binary_expression<T, order, LHS, RHS, fmod_expr<T, order, LHS, RHS>>(
              left_hand_expr, right_hand_expr)
    {}

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return fmod(this->lhs.evaluate(), this->rhs.evaluate());
    };
    static const inner_t left_derivative(const inner_t &l, const inner_t &r, const inner_t &v)
    {
        return inner_t{1.0};
    };
    static const inner_t right_derivative(const inner_t &l, const inner_t &r, const inner_t &v)
    {
        BOOST_MATH_STD_USING
        return static_cast<T>(-1.0) * trunc(l / r);
    };
};

template<typename T, size_t order, typename ARG>
struct fmod_left_float_expr
    : public abstract_unary_expression<T, order, ARG, fmod_left_float_expr<T, order, ARG>>
{
    /** @brief 
      * */
    using inner_t    = rvar_t<T, order - 1>;

    explicit fmod_left_float_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, fmod_left_float_expr<T, order, ARG>>(arg_expr,
                                                                                        v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return fmod(this->constant, this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        return static_cast<T>(-1.0) * trunc(constant / argv);
    }
};

template<typename T, size_t order, typename ARG>
struct fmod_right_float_expr
    : public abstract_unary_expression<T, order, ARG, fmod_right_float_expr<T, order, ARG>>
{
    /** @brief
      * */
    using inner_t    = rvar_t<T, order - 1>;

    explicit fmod_right_float_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, fmod_right_float_expr<T, order, ARG>>(arg_expr,
                                                                                         v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return fmod(this->arg.evaluate(), this->constant);
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        return inner_t{1.0};
    }
};
/**************************************************************************************************/
template<typename T, size_t order, typename ARG>
fabs_expr<T, order, ARG> fabs(const expression<T, order, ARG> &arg)
{
    return fabs_expr<T, order, ARG>(arg, static_cast<T>(0.0));
}
template<typename T, size_t order, typename ARG>
auto abs(const expression<T, order, ARG> &arg)
{
    return fabs(arg);
}
template<typename T, size_t order, typename ARG>
ceil_expr<T, order, ARG> ceil(const expression<T, order, ARG> &arg)
{
    return ceil_expr<T, order, ARG>(arg, static_cast<T>(0.0));
}
template<typename T, size_t order, typename ARG>
floor_expr<T, order, ARG> floor(const expression<T, order, ARG> &arg)
{
    return floor_expr<T, order, ARG>(arg, static_cast<T>(0.0));
}

template<typename T, size_t order, typename ARG>
exp_expr<T, order, ARG> exp(const expression<T, order, ARG> &arg)
{
    return exp_expr<T, order, ARG>(arg, static_cast<T>(0.0));
}

template<typename T, size_t order, typename LHS, typename RHS>
pow_expr<T, order, LHS, RHS> pow(const expression<T, order, LHS> &lhs,
                                 const expression<T, order, RHS> &rhs)
{
    return pow_expr<T, order, LHS, RHS>(lhs, rhs);
}

template<typename U,
         typename T,
         size_t order,
         typename ARG,
         typename = typename std::enable_if<!detail::is_expression<U>::value>::type>
expr_pow_float_expr<T, order, ARG> pow(const expression<T, order, ARG> &arg, const U &v)
{
    return expr_pow_float_expr<T, order, ARG>(arg, static_cast<T>(v));
};

template<typename T, size_t order, typename ARG>
float_pow_expr_expr<T, order, ARG> pow(const T &v, const expression<T, order, ARG> &arg)
{
    return float_pow_expr_expr<T, order, ARG>(arg, v);
};

template<typename T, size_t order, typename ARG>
log_expr<T, order, ARG> log(const expression<T, order, ARG> &arg)
{
    return log_expr<T, order, ARG>(arg, static_cast<T>(0.0));
};

template<typename T, size_t order, typename ARG>
sqrt_expr<T, order, ARG> sqrt(const expression<T, order, ARG> &arg)
{
    return sqrt_expr<T, order, ARG>(arg, static_cast<T>(0.0));
};

template<typename T, size_t order, typename ARG>
auto frexp(const expression<T, order, ARG> &arg, int *i)
{
    BOOST_MATH_STD_USING
    frexp(arg.evaluate(), i);
    return arg / pow(static_cast<T>(2.0), *i);
}

template<typename T, size_t order, typename ARG>
auto ldexp(const expression<T, order, ARG> &arg, const int &i)
{
    BOOST_MATH_STD_USING
    return arg * pow(static_cast<T>(2.0), i);
}

template<typename T, size_t order, typename ARG>
cos_expr<T, order, ARG> cos(const expression<T, order, ARG> &arg)
{
    return cos_expr<T, order, ARG>(arg, static_cast<T>(0.0));
};

template<typename T, size_t order, typename ARG>
sin_expr<T, order, ARG> sin(const expression<T, order, ARG> &arg)
{
    return sin_expr<T, order, ARG>(arg, static_cast<T>(0.0));
};

template<typename T, size_t order, typename ARG>
tan_expr<T, order, ARG> tan(const expression<T, order, ARG> &arg)
{
    return tan_expr<T, order, ARG>(arg, static_cast<T>(0.0));
};

template<typename T, size_t order, typename ARG>
acos_expr<T, order, ARG> acos(const expression<T, order, ARG> &arg)
{
    return acos_expr<T, order, ARG>(arg, static_cast<T>(0.0));
};

template<typename T, size_t order, typename ARG>
asin_expr<T, order, ARG> asin(const expression<T, order, ARG> &arg)
{
    return asin_expr<T, order, ARG>(arg, static_cast<T>(0.0));
};

template<typename T, size_t order, typename ARG>
atan_expr<T, order, ARG> atan(const expression<T, order, ARG> &arg)
{
    return atan_expr<T, order, ARG>(arg, static_cast<T>(0.0));
};

template<typename T, size_t order, typename LHS, typename RHS>
atan2_expr<T, order, LHS, RHS> atan2(const expression<T, order, LHS> &lhs,
                                     const expression<T, order, RHS> &rhs)
{
    return atan2_expr<T, order, LHS, RHS>(lhs, rhs);
}

template<typename T, size_t order, typename ARG>
atan2_right_float_expr<T, order, ARG> atan2(const expression<T, order, ARG> &arg, const T &v)
{
    return atan2_right_float_expr<T, order, ARG>(arg, v);
};

template<typename T, size_t order, typename ARG>
atan2_left_float_expr<T, order, ARG> atan2(const T &v, const expression<T, order, ARG> &arg)
{
    return atan2_left_float_expr<T, order, ARG>(arg, v);
};

template<typename T, size_t order, typename ARG>
trunc_expr<T, order, ARG> trunc(const expression<T, order, ARG> &arg)
{
    return trunc_expr<T, order, ARG>(arg, static_cast<T>(0.0));
}

template<typename T, size_t order, typename LHS, typename RHS>
auto fmod(const expression<T, order, LHS> &lhs, const expression<T, order, RHS> &rhs)
{
    return fmod_expr<T, order, LHS, RHS>(lhs, rhs);
}

template<typename T, size_t order, typename ARG>
auto fmod(const expression<T, order, ARG> &lhs, const T rhs)
{
    return fmod_right_float_expr<T, order, ARG>(lhs, rhs);
}

template<typename T, size_t order, typename ARG>
auto fmod(const T lhs, const expression<T, order, ARG> &rhs)
{
    return fmod_left_float_expr<T, order, ARG>(rhs, lhs);
}

template<typename T, size_t order, typename ARG>
round_expr<T, order, ARG> round(const expression<T, order, ARG> &arg)
{
    return round_expr<T, order, ARG>(arg, static_cast<T>(0.0));
}

template<typename T, size_t order, typename ARG>
int iround(const expression<T, order, ARG> &arg)
{
    rvar<T, order> tmp = arg.evaluate();
    return iround(tmp.item());
}
template<typename T, size_t order, typename ARG>
long lround(const expression<T, order, ARG> &arg)
{
    BOOST_MATH_STD_USING
    rvar<T, order> tmp = arg.evaluate();
    return lround(tmp.item());
}

template<typename T, size_t order, typename ARG>
long long llround(const expression<T, order, ARG> &arg)
{
    rvar<T, order> tmp = arg.evaluate();
    return llround(tmp.item());
}

template<typename T, size_t order, typename ARG>
int itrunc(const expression<T, order, ARG> &arg)
{
    rvar<T, order> tmp = arg.evaluate();
    return itrunc(tmp.item());
}

template<typename T, size_t order, typename ARG>
long ltrunc(const expression<T, order, ARG> &arg)
{
    rvar<T, order> tmp = arg.evaluate();
    return ltrunc(tmp.item());
}

template<typename T, size_t order, typename ARG>
long long lltrunc(const expression<T, order, ARG> &arg)
{
    rvar<T, order> tmp = arg.evaluate();
    return lltrunc(tmp.item());
}

template<typename T, size_t order, typename ARG>
sinh_expr<T, order, ARG> sinh(const expression<T, order, ARG> &arg)
{
    return sinh_expr<T, order, ARG>(arg, static_cast<T>(0.0));
}

template<typename T, size_t order, typename ARG>
cosh_expr<T, order, ARG> cosh(const expression<T, order, ARG> &arg)
{
    return cosh_expr<T, order, ARG>(arg, static_cast<T>(0.0));
}

template<typename T, size_t order, typename ARG>
tanh_expr<T, order, ARG> tanh(const expression<T, order, ARG> &arg)
{
    return tanh_expr<T, order, ARG>(arg, static_cast<T>(0.0));
}

template<typename T, size_t order, typename ARG>
log10_expr<T, order, ARG> log10(const expression<T, order, ARG> &arg)
{
    return log10_expr<T, order, ARG>(arg, static_cast<T>(0.0));
}
template<typename T, size_t order, typename ARG>
asinh_expr<T, order, ARG> asinh(const expression<T, order, ARG> &arg)
{
    return asinh_expr<T, order, ARG>(arg, static_cast<T>(0.0));
}
template<typename T, size_t order, typename ARG>
acosh_expr<T, order, ARG> acosh(const expression<T, order, ARG> &arg)
{
    return acosh_expr<T, order, ARG>(arg, static_cast<T>(0.0));
}
template<typename T, size_t order, typename ARG>
atanh_expr<T, order, ARG> atanh(const expression<T, order, ARG> &arg)
{
    return atanh_expr<T, order, ARG>(arg, static_cast<T>(0.0));
}
} // namespace reverse_mode
} // namespace differentiation
} // namespace math
} // namespace boost
#endif
