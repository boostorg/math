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
template<typename RealType, size_t DerivativeOrder, typename ARG>
struct fabs_expr : public abstract_unary_expression<RealType, DerivativeOrder, ARG, fabs_expr<RealType,DerivativeOrder, ARG>>
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
    using inner_t = rvar_t<RealType, DerivativeOrder - 1>;

    explicit fabs_expr(const expression<RealType, DerivativeOrder, ARG> &arg_expr, const RealType &v)
        : abstract_unary_expression<RealType,
                                    DerivativeOrder,
                                    ARG,
                                    fabs_expr<RealType, DerivativeOrder, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return fabs(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv,
                                    const inner_t & /*v*/,
                                    const RealType & /*constant*/)
    {
        return argv > 0.0 ? inner_t{1.0} : inner_t{-1.0};
    }
};

template<typename RealType, size_t DerivativeOrder, typename ARG>
struct ceil_expr : public abstract_unary_expression<RealType, DerivativeOrder, ARG, ceil_expr<RealType,DerivativeOrder, ARG>>
{
    /** @brief ceil(1.11) = 2.0
    *
    * d/dx ceil(x) = 0.0 for all x
    *
    * we avoid problematic points at x = 1,2,3...
    * as with optimization its most likely intented
    * this function's derivative is 0.0;
    * */
    using inner_t    = rvar_t<RealType, DerivativeOrder - 1>;

    explicit ceil_expr(const expression<RealType, DerivativeOrder, ARG> &arg_expr, const RealType &v)
        : abstract_unary_expression<RealType,
                                    DerivativeOrder,
                                    ARG,
                                    ceil_expr<RealType, DerivativeOrder, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return ceil(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t & /*argv*/,
                                    const inner_t & /*v*/,
                                    const RealType & /*constant*/)
    {
        return inner_t{0.0};
    }
};

template<typename RealType, size_t DerivativeOrder, typename ARG>
struct floor_expr : public abstract_unary_expression<RealType, DerivativeOrder, ARG, floor_expr<RealType,DerivativeOrder, ARG>>
{
    /** @brief floor(1.11) = 1.0, floor(-1.11) = 2
    *
    * d/dx floor(x) = 0.0 for all x
    *
    * we avoid problematic points at x = 1,2,3...
    * as with optimization its most likely intented
    * this function's derivative is 0.0;
    * */
    using inner_t    = rvar_t<RealType, DerivativeOrder - 1>;

    explicit floor_expr(const expression<RealType, DerivativeOrder, ARG> &arg_expr,
                        const RealType                                   &v)
        : abstract_unary_expression<RealType,
                                    DerivativeOrder,
                                    ARG,
                                    floor_expr<RealType, DerivativeOrder, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return floor(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t & /*argv*/,
                                    const inner_t & /*v*/,
                                    const RealType & /*constant*/)
    {
        return inner_t{0.0};
    }
};

template<typename RealType, size_t DerivativeOrder, typename ARG>
struct trunc_expr : public abstract_unary_expression<RealType, DerivativeOrder, ARG, trunc_expr<RealType,DerivativeOrder, ARG>>
{
    /** @brief trunc(1.11) = 1.0, trunc(-1.11) = -1.0
    *
    * d/dx trunc(x) = 0.0 for all x
    *
    * we avoid problematic points at x = 1,2,3...
    * as with optimization its most likely intented
    * this function's derivative is 0.0;
    * */
    using inner_t    = rvar_t<RealType, DerivativeOrder - 1>;

    explicit trunc_expr(const expression<RealType, DerivativeOrder, ARG> &arg_expr,
                        const RealType                                   &v)
        : abstract_unary_expression<RealType,
                                    DerivativeOrder,
                                    ARG,
                                    trunc_expr<RealType, DerivativeOrder, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return trunc(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t & /*argv*/,
                                    const inner_t & /*v*/,
                                    const RealType & /*constant*/)
    {
        return inner_t{0.0};
    }
};

template<typename RealType, size_t DerivativeOrder, typename ARG>
struct exp_expr : public abstract_unary_expression<RealType, DerivativeOrder, ARG, exp_expr<RealType,DerivativeOrder, ARG>>
{
    /** @brief exp(x)
    *
    * d/dx exp(x) = exp(x)
    *
    * */
    using inner_t    = rvar_t<RealType, DerivativeOrder - 1>;

    explicit exp_expr(const expression<RealType, DerivativeOrder, ARG> &arg_expr, const RealType &v)
        : abstract_unary_expression<RealType,
                                    DerivativeOrder,
                                    ARG,
                                    exp_expr<RealType, DerivativeOrder, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return exp(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv,
                                    const inner_t & /*v*/,
                                    const RealType & /*constant*/)
    {
        BOOST_MATH_STD_USING
        return exp(argv);
    }
};

template<typename RealType, size_t DerivativeOrder, typename LHS, typename RHS>
struct pow_expr
    : public abstract_binary_expression<RealType, DerivativeOrder, LHS, RHS, pow_expr<RealType, DerivativeOrder, LHS, RHS>>
{
    /** @brief pow(x,y)
     *  d/dx pow(x,y) = y pow (x, y-1)
     *  d/dy pow(x,y) = pow(x,y) log(x)
    * */
    using inner_t    = rvar_t<RealType, DerivativeOrder - 1>;
    // Explicitly define constructor to forward to base class
    explicit pow_expr(const expression<RealType, DerivativeOrder, LHS> &left_hand_expr,
                      const expression<RealType, DerivativeOrder, RHS> &right_hand_expr)
        : abstract_binary_expression<RealType, DerivativeOrder, LHS, RHS, pow_expr<RealType, DerivativeOrder, LHS, RHS>>(
              left_hand_expr, right_hand_expr)
    {}

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return pow(this->lhs.evaluate(), this->rhs.evaluate());
    };
    static const inner_t left_derivative(const inner_t &l, const inner_t &r, const inner_t & /*v*/)
    {
        BOOST_MATH_STD_USING
        return r * pow(l, r - static_cast<RealType>(1.0));
    };
    static const inner_t right_derivative(const inner_t &l, const inner_t &r, const inner_t & /*v*/)
    {
        BOOST_MATH_STD_USING
        return pow(l, r) * log(l);
    };
};

template<typename RealType, size_t DerivativeOrder, typename ARG>
struct expr_pow_float_expr
    : public abstract_unary_expression<RealType, DerivativeOrder, ARG, expr_pow_float_expr<RealType,DerivativeOrder, ARG>>
{
    /** @brief pow(rvar,float)
      */
    using inner_t    = rvar_t<RealType, DerivativeOrder - 1>;

    explicit expr_pow_float_expr(const expression<RealType, DerivativeOrder, ARG> &arg_expr,
                                 const RealType                                   &v)
        : abstract_unary_expression<RealType,
                                    DerivativeOrder,
                                    ARG,
                                    expr_pow_float_expr<RealType, DerivativeOrder, ARG>>(arg_expr,
                                                                                         v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return pow(this->arg.evaluate(), this->constant);
    }
    static const inner_t derivative(const inner_t &argv, const inner_t & /*v*/, const RealType &constant)
    {
        BOOST_MATH_STD_USING
        return inner_t{constant} * pow(argv, inner_t{constant - 1});
    }
};

template<typename RealType, size_t DerivativeOrder, typename ARG>
struct float_pow_expr_expr
    : public abstract_unary_expression<RealType, DerivativeOrder, ARG, float_pow_expr_expr<RealType,DerivativeOrder, ARG>>
{
    /** @brief pow(float, rvar)
      * */
    using inner_t    = rvar_t<RealType, DerivativeOrder - 1>;

    explicit float_pow_expr_expr(const expression<RealType, DerivativeOrder, ARG> &arg_expr,
                                 const RealType                                   &v)
        : abstract_unary_expression<RealType,
                                    DerivativeOrder,
                                    ARG,
                                    float_pow_expr_expr<RealType, DerivativeOrder, ARG>>(arg_expr,
                                                                                         v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return pow(this->constant, this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t & /*v*/, const RealType &constant)
    {
        BOOST_MATH_STD_USING
        return pow(constant, argv) * log(constant);
    }
};

template<typename RealType, size_t DerivativeOrder, typename ARG>
struct sqrt_expr : public abstract_unary_expression<RealType, DerivativeOrder, ARG, sqrt_expr<RealType,DerivativeOrder, ARG>>
{
    /** @brief  sqrt(x)
     *  d/dx sqrt(x) = 1/(2 sqrt(x))
    * */
    using inner_t    = rvar_t<RealType, DerivativeOrder - 1>;

    explicit sqrt_expr(const expression<RealType, DerivativeOrder, ARG> &arg_expr, const RealType &v)
        : abstract_unary_expression<RealType,
                                    DerivativeOrder,
                                    ARG,
                                    sqrt_expr<RealType, DerivativeOrder, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return sqrt(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv,
                                    const inner_t & /*v*/,
                                    const RealType & /*constant*/)
    {
        BOOST_MATH_STD_USING
        return static_cast<RealType>(1.0) / (static_cast<RealType>(2.0) * sqrt(argv));
    }
};

template<typename RealType, size_t DerivativeOrder, typename ARG>
struct log_expr : public abstract_unary_expression<RealType, DerivativeOrder, ARG, log_expr<RealType,DerivativeOrder, ARG>>
{
    /** @brief log(x)
     *  d/dx log(x) = 1/x
      * */
    using inner_t    = rvar_t<RealType, DerivativeOrder - 1>;

    explicit log_expr(const expression<RealType, DerivativeOrder, ARG> &arg_expr, const RealType &v)
        : abstract_unary_expression<RealType,
                                    DerivativeOrder,
                                    ARG,
                                    log_expr<RealType, DerivativeOrder, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return log(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv,
                                    const inner_t & /*v*/,
                                    const RealType & /*constant*/)
    {
        return static_cast<RealType>(1.0) / argv;
    }
};

template<typename RealType, size_t DerivativeOrder, typename ARG>
struct cos_expr : public abstract_unary_expression<RealType, DerivativeOrder, ARG, cos_expr<RealType,DerivativeOrder, ARG>>
{
    /** @brief cos(x)
     *  d/dx cos(x) = -sin(x)
      * */
    using inner_t    = rvar_t<RealType, DerivativeOrder - 1>;

    explicit cos_expr(const expression<RealType, DerivativeOrder, ARG> &arg_expr, const RealType &v)
        : abstract_unary_expression<RealType,
                                    DerivativeOrder,
                                    ARG,
                                    cos_expr<RealType, DerivativeOrder, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return cos(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv,
                                    const inner_t & /*v*/,
                                    const RealType & /*constant*/)
    {
        BOOST_MATH_STD_USING
        return -sin(argv);
    }
};

template<typename RealType, size_t DerivativeOrder, typename ARG>
struct sin_expr : public abstract_unary_expression<RealType, DerivativeOrder, ARG, sin_expr<RealType,DerivativeOrder, ARG>>
{
    /** @brief sin(x)
     *  d/dx sin(x) = cos(x)
      * */
    using arg_type   = ARG;
    using value_type = RealType;
    using inner_t    = rvar_t<RealType, DerivativeOrder - 1>;

    explicit sin_expr(const expression<RealType, DerivativeOrder, ARG> &arg_expr, const RealType &v)
        : abstract_unary_expression<RealType,
                                    DerivativeOrder,
                                    ARG,
                                    sin_expr<RealType, DerivativeOrder, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return sin(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv,
                                    const inner_t & /*v*/,
                                    const RealType & /*constant*/)
    {
        BOOST_MATH_STD_USING
        return cos(argv);
    }
};

template<typename RealType, size_t DerivativeOrder, typename ARG>
struct tan_expr : public abstract_unary_expression<RealType, DerivativeOrder, ARG, tan_expr<RealType,DerivativeOrder, ARG>>
{
    /** @brief tan(x)
     *  d/dx tan(x) = 1/cos^2(x)
      * */
    using inner_t    = rvar_t<RealType, DerivativeOrder - 1>;

    explicit tan_expr(const expression<RealType, DerivativeOrder, ARG> &arg_expr, const RealType &v)
        : abstract_unary_expression<RealType, DerivativeOrder, ARG, tan_expr<RealType,DerivativeOrder, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return tan(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv,
                                    const inner_t & /*v*/,
                                    const RealType & /*constant*/)
    {
        BOOST_MATH_STD_USING
        return static_cast<RealType>(1.0) / (cos(argv) * cos(argv));
    }
};

template<typename RealType, size_t DerivativeOrder, typename ARG>
struct acos_expr : public abstract_unary_expression<RealType, DerivativeOrder, ARG, acos_expr<RealType,DerivativeOrder, ARG>>
{
    /** @brief acos(x)
     *  d/dx acos(x) = -1/sqrt(1-x^2)
      * */
    using inner_t    = rvar_t<RealType, DerivativeOrder - 1>;

    explicit acos_expr(const expression<RealType, DerivativeOrder, ARG> &arg_expr, const RealType &v)
        : abstract_unary_expression<RealType, DerivativeOrder, ARG, acos_expr<RealType,DerivativeOrder, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return acos(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv,
                                    const inner_t & /*v*/,
                                    const RealType & /*constant*/)
    {
        BOOST_MATH_STD_USING
        return static_cast<RealType>(-1.0) / sqrt(static_cast<RealType>(1.0) - argv * argv);
    }
};

template<typename RealType, size_t DerivativeOrder, typename ARG>
struct asin_expr : public abstract_unary_expression<RealType, DerivativeOrder, ARG, asin_expr<RealType,DerivativeOrder, ARG>>
{
    /** @brief asin(x)
     *  d/dx asin =  1/sqrt(1-x^2)
      * */
    using arg_type   = ARG;
    using value_type = RealType;
    using inner_t    = rvar_t<RealType, DerivativeOrder - 1>;

    explicit asin_expr(const expression<RealType, DerivativeOrder, ARG> &arg_expr, const RealType &v)
        : abstract_unary_expression<RealType, DerivativeOrder, ARG, asin_expr<RealType,DerivativeOrder, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return asin(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv,
                                    const inner_t & /*v*/,
                                    const RealType & /*constant*/)
    {
        BOOST_MATH_STD_USING
        return static_cast<RealType>(1.0) / sqrt(static_cast<RealType>(1.0) - argv * argv);
    }
};

template<typename RealType, size_t DerivativeOrder, typename ARG>
struct atan_expr : public abstract_unary_expression<RealType, DerivativeOrder, ARG, atan_expr<RealType,DerivativeOrder, ARG>>
{
    /** @brief atan(x)
     *  d/dx atan(x) = 1/x^2+1
      * */
    using inner_t    = rvar_t<RealType, DerivativeOrder - 1>;

    explicit atan_expr(const expression<RealType, DerivativeOrder, ARG> &arg_expr, const RealType &v)
        : abstract_unary_expression<RealType, DerivativeOrder, ARG, atan_expr<RealType,DerivativeOrder, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return atan(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv,
                                    const inner_t & /*v*/,
                                    const RealType & /*constant*/)
    {
        BOOST_MATH_STD_USING
        return static_cast<RealType>(1.0) / (static_cast<RealType>(1.0) + argv * argv);
    }
};
template<typename RealType, size_t DerivativeOrder, typename LHS, typename RHS>
struct atan2_expr
    : public abstract_binary_expression<RealType, DerivativeOrder, LHS, RHS, atan2_expr<RealType, DerivativeOrder, LHS, RHS>>
{
    /** @brief atan2(x,y)
    * */
    using inner_t    = rvar_t<RealType, DerivativeOrder - 1>;
    // Explicitly define constructor to forward to base class
    explicit atan2_expr(const expression<RealType, DerivativeOrder, LHS> &left_hand_expr,
                        const expression<RealType, DerivativeOrder, RHS> &right_hand_expr)
        : abstract_binary_expression<RealType, DerivativeOrder, LHS, RHS, atan2_expr<RealType, DerivativeOrder, LHS, RHS>>(
              left_hand_expr, right_hand_expr)
    {}

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return atan2(this->lhs.evaluate(), this->rhs.evaluate());
    };
    static const inner_t left_derivative(const inner_t &l, const inner_t &r, const inner_t & /*v*/)
    {
        return r / (l * l + r * r);
    };
    static const inner_t right_derivative(const inner_t &l, const inner_t &r, const inner_t & /*v*/)
    {
        return -l / (l * l + r * r);
    };
};

template<typename RealType, size_t DerivativeOrder, typename ARG>
struct atan2_left_float_expr
    : public abstract_unary_expression<RealType, DerivativeOrder, ARG, atan2_left_float_expr<RealType,DerivativeOrder, ARG>>
{
    /** @brief atan2(float,rvar) 
      * */
    using inner_t    = rvar_t<RealType, DerivativeOrder - 1>;

    explicit atan2_left_float_expr(const expression<RealType, DerivativeOrder, ARG> &arg_expr, const RealType &v)
        : abstract_unary_expression<RealType, DerivativeOrder, ARG, atan2_left_float_expr<RealType,DerivativeOrder, ARG>>(arg_expr,
                                                                                         v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return atan2(this->constant, this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t & /*v*/, const RealType &constant)
    {
        return -constant / (constant * constant + argv * argv);
    }
};

template<typename RealType, size_t DerivativeOrder, typename ARG>
struct atan2_right_float_expr
    : public abstract_unary_expression<RealType, DerivativeOrder, ARG, atan2_right_float_expr<RealType,DerivativeOrder, ARG>>
{
    /** @brief atan2(rvar,float) 
      * */
    using inner_t    = rvar_t<RealType, DerivativeOrder - 1>;

    explicit atan2_right_float_expr(const expression<RealType, DerivativeOrder, ARG> &arg_expr, const RealType &v)
        : abstract_unary_expression<RealType, DerivativeOrder, ARG, atan2_right_float_expr<RealType,DerivativeOrder, ARG>>(arg_expr,
                                                                                          v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return atan2(this->arg.evaluate(), this->constant);
    }
    static const inner_t derivative(const inner_t &argv, const inner_t & /*v*/, const RealType &constant)
    {
        return constant / (constant * constant + argv * argv);
    }
};

template<typename RealType, size_t DerivativeOrder, typename ARG>
struct round_expr : public abstract_unary_expression<RealType, DerivativeOrder, ARG, round_expr<RealType,DerivativeOrder, ARG>>
{
    /** @brief round(x)
     *  d/dx round = 0
      * */
    using inner_t    = rvar_t<RealType, DerivativeOrder - 1>;

    explicit round_expr(const expression<RealType, DerivativeOrder, ARG> &arg_expr, const RealType &v)
        : abstract_unary_expression<RealType, DerivativeOrder, ARG, round_expr<RealType,DerivativeOrder, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return round(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t & /*argv*/,
                                    const inner_t & /*v*/,
                                    const RealType & /*constant*/)
    {
        return inner_t{0.0};
    }
};

template<typename RealType, size_t DerivativeOrder, typename ARG>
struct sinh_expr : public abstract_unary_expression<RealType, DerivativeOrder, ARG, sinh_expr<RealType,DerivativeOrder, ARG>>
{
    /** @brief sinh(x)
     *  d/dx sinh(x) = cosh
      * */
    using inner_t    = rvar_t<RealType, DerivativeOrder - 1>;

    explicit sinh_expr(const expression<RealType, DerivativeOrder, ARG> &arg_expr, const RealType &v)
        : abstract_unary_expression<RealType, DerivativeOrder, ARG, sinh_expr<RealType,DerivativeOrder, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return sinh(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv,
                                    const inner_t & /*v*/,
                                    const RealType & /*constant*/)
    {
        BOOST_MATH_STD_USING
        return cosh(argv);
    }
};

template<typename RealType, size_t DerivativeOrder, typename ARG>
struct cosh_expr : public abstract_unary_expression<RealType, DerivativeOrder, ARG, cosh_expr<RealType,DerivativeOrder, ARG>>
{
    /** @brief cosh(x)
     *  d/dx cosh(x) = sinh
      * */
    using inner_t    = rvar_t<RealType, DerivativeOrder - 1>;

    explicit cosh_expr(const expression<RealType, DerivativeOrder, ARG> &arg_expr, const RealType &v)
        : abstract_unary_expression<RealType, DerivativeOrder, ARG, cosh_expr<RealType,DerivativeOrder, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return cosh(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv,
                                    const inner_t & /*v*/,
                                    const RealType & /*constant*/)
    {
        BOOST_MATH_STD_USING
        return sinh(argv);
    }
};
template<typename RealType, size_t DerivativeOrder, typename ARG>
struct tanh_expr : public abstract_unary_expression<RealType, DerivativeOrder, ARG, tanh_expr<RealType,DerivativeOrder, ARG>>
{
    /** @brief tanh(x)
     *  d/dx tanh(x) = 1/cosh^2
      * */
    using inner_t    = rvar_t<RealType, DerivativeOrder - 1>;

    explicit tanh_expr(const expression<RealType, DerivativeOrder, ARG> &arg_expr, const RealType &v)
        : abstract_unary_expression<RealType, DerivativeOrder, ARG, tanh_expr<RealType,DerivativeOrder, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return tanh(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv,
                                    const inner_t & /*v*/,
                                    const RealType & /*constant*/)
    {
        BOOST_MATH_STD_USING
        return static_cast<RealType>(1.0) / (cosh(argv) * cosh(argv));
    }
};

template<typename RealType, size_t DerivativeOrder, typename ARG>
struct log10_expr : public abstract_unary_expression<RealType, DerivativeOrder, ARG, log10_expr<RealType,DerivativeOrder, ARG>>
{
    /** @brief log10(x)
     *  d/dx log10(x) = 1/(x * log(10))
      * */
    using inner_t    = rvar_t<RealType, DerivativeOrder - 1>;

    explicit log10_expr(const expression<RealType, DerivativeOrder, ARG> &arg_expr, const RealType &v)
        : abstract_unary_expression<RealType, DerivativeOrder, ARG, log10_expr<RealType,DerivativeOrder, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return log10(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv,
                                    const inner_t & /*v*/,
                                    const RealType & /*constant*/)
    {
        BOOST_MATH_STD_USING
        return static_cast<RealType>(1.0) / (argv * log(static_cast<RealType>(10.0)));
    }
};

template<typename RealType, size_t DerivativeOrder, typename ARG>
struct acosh_expr : public abstract_unary_expression<RealType, DerivativeOrder, ARG, acosh_expr<RealType,DerivativeOrder, ARG>>
{
    /** @brief acosh(x)
     *  d/dx acosh(x) = 1/(sqrt(x-1)sqrt(x+1)
      * */
    using inner_t    = rvar_t<RealType, DerivativeOrder - 1>;

    explicit acosh_expr(const expression<RealType, DerivativeOrder, ARG> &arg_expr, const RealType &v)
        : abstract_unary_expression<RealType, DerivativeOrder, ARG, acosh_expr<RealType,DerivativeOrder, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return acosh(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv,
                                    const inner_t & /*v*/,
                                    const RealType & /*constant*/)
    {
        BOOST_MATH_STD_USING
        return static_cast<RealType>(1.0)
               / (sqrt(argv - static_cast<RealType>(1.0)) * sqrt(argv + static_cast<RealType>(1.0)));
    }
};

template<typename RealType, size_t DerivativeOrder, typename ARG>
struct asinh_expr : public abstract_unary_expression<RealType, DerivativeOrder, ARG, asinh_expr<RealType,DerivativeOrder, ARG>>
{
    /** @brief asinh(x)
     *  d/dx asinh(x) = 1/(sqrt(1+x^2))
      * */
    using inner_t    = rvar_t<RealType, DerivativeOrder - 1>;

    explicit asinh_expr(const expression<RealType, DerivativeOrder, ARG> &arg_expr, const RealType &v)
        : abstract_unary_expression<RealType, DerivativeOrder, ARG, asinh_expr<RealType,DerivativeOrder, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return asinh(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv,
                                    const inner_t & /*v*/,
                                    const RealType & /*constant*/)
    {
        BOOST_MATH_STD_USING
        return static_cast<RealType>(1.0) / (sqrt(static_cast<RealType>(1.0) + argv * argv));
    }
};

template<typename RealType, size_t DerivativeOrder, typename ARG>
struct atanh_expr : public abstract_unary_expression<RealType, DerivativeOrder, ARG, atanh_expr<RealType,DerivativeOrder, ARG>>
{
    /** @brief atanh(x)
     *  d/dx atanh(x) = 1/(1-x^2)
      * */
    using inner_t    = rvar_t<RealType, DerivativeOrder - 1>;

    explicit atanh_expr(const expression<RealType, DerivativeOrder, ARG> &arg_expr, const RealType &v)
        : abstract_unary_expression<RealType, DerivativeOrder, ARG, atanh_expr<RealType,DerivativeOrder, ARG>>(arg_expr, v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return atanh(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv,
                                    const inner_t & /*v*/,
                                    const RealType & /*constant*/)
    {
        BOOST_MATH_STD_USING
        return static_cast<RealType>(1.0) / (static_cast<RealType>(1.0) - argv * argv);
    }
};
template<typename RealType, size_t DerivativeOrder, typename LHS, typename RHS>
struct fmod_expr
    : public abstract_binary_expression<RealType, DerivativeOrder, LHS, RHS, fmod_expr<RealType, DerivativeOrder, LHS, RHS>>
{
    /** @brief 
    * */
    using inner_t    = rvar_t<RealType, DerivativeOrder - 1>;
    // Explicitly define constructor to forward to base class
    explicit fmod_expr(const expression<RealType, DerivativeOrder, LHS> &left_hand_expr,
                       const expression<RealType, DerivativeOrder, RHS> &right_hand_expr)
        : abstract_binary_expression<RealType, DerivativeOrder, LHS, RHS, fmod_expr<RealType, DerivativeOrder, LHS, RHS>>(
              left_hand_expr, right_hand_expr)
    {}

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return fmod(this->lhs.evaluate(), this->rhs.evaluate());
    };
    static const inner_t left_derivative(const inner_t & /*l*/,
                                         const inner_t & /*r*/,
                                         const inner_t & /*v*/)
    {
        return inner_t{1.0};
    };
    static const inner_t right_derivative(const inner_t &l, const inner_t &r, const inner_t & /*v*/)
    {
        BOOST_MATH_STD_USING
        return static_cast<RealType>(-1.0) * trunc(l / r);
    };
};

template<typename RealType, size_t DerivativeOrder, typename ARG>
struct fmod_left_float_expr
    : public abstract_unary_expression<RealType, DerivativeOrder, ARG, fmod_left_float_expr<RealType,DerivativeOrder, ARG>>
{
    /** @brief 
      * */
    using inner_t    = rvar_t<RealType, DerivativeOrder - 1>;

    explicit fmod_left_float_expr(const expression<RealType, DerivativeOrder, ARG> &arg_expr, const RealType &v)
        : abstract_unary_expression<RealType, DerivativeOrder, ARG, fmod_left_float_expr<RealType,DerivativeOrder, ARG>>(arg_expr,
                                                                                        v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return fmod(this->constant, this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t & /*v*/, const RealType &constant)
    {
        return static_cast<RealType>(-1.0) * trunc(constant / argv);
    }
};

template<typename RealType, size_t DerivativeOrder, typename ARG>
struct fmod_right_float_expr
    : public abstract_unary_expression<RealType, DerivativeOrder, ARG, fmod_right_float_expr<RealType,DerivativeOrder, ARG>>
{
    /** @brief
      * */
    using inner_t    = rvar_t<RealType, DerivativeOrder - 1>;

    explicit fmod_right_float_expr(const expression<RealType, DerivativeOrder, ARG> &arg_expr,
                                   const RealType                                   &v)
        : abstract_unary_expression<RealType,
                                    DerivativeOrder,
                                    ARG,
                                    fmod_right_float_expr<RealType, DerivativeOrder, ARG>>(arg_expr,
                                                                                           v){};

    inner_t evaluate() const
    {
        BOOST_MATH_STD_USING
        return fmod(this->arg.evaluate(), this->constant);
    }
    static const inner_t derivative(const inner_t & /*argv*/,
                                    const inner_t & /*v*/,
                                    const RealType & /*constant*/)
    {
        return inner_t{1.0};
    }
};
/**************************************************************************************************/
template<typename RealType, size_t DerivativeOrder, typename ARG>
fabs_expr<RealType, DerivativeOrder, ARG> fabs(const expression<RealType,DerivativeOrder, ARG> &arg)
{
    return fabs_expr<RealType, DerivativeOrder, ARG>(arg, static_cast<RealType>(0.0));
}
template<typename RealType, size_t DerivativeOrder, typename ARG>
auto abs(const expression<RealType, DerivativeOrder, ARG> &arg)
{
    return fabs(arg);
}
template<typename RealType, size_t DerivativeOrder, typename ARG>
ceil_expr<RealType, DerivativeOrder, ARG> ceil(const expression<RealType,DerivativeOrder, ARG> &arg)
{
    return ceil_expr<RealType, DerivativeOrder, ARG>(arg, static_cast<RealType>(0.0));
}
template<typename RealType, size_t DerivativeOrder, typename ARG>
floor_expr<RealType, DerivativeOrder, ARG> floor(const expression<RealType,DerivativeOrder, ARG> &arg)
{
    return floor_expr<RealType, DerivativeOrder, ARG>(arg, static_cast<RealType>(0.0));
}

template<typename RealType, size_t DerivativeOrder, typename ARG>
exp_expr<RealType, DerivativeOrder, ARG> exp(const expression<RealType,DerivativeOrder, ARG> &arg)
{
    return exp_expr<RealType, DerivativeOrder, ARG>(arg, static_cast<RealType>(0.0));
}

template<typename RealType, size_t DerivativeOrder, typename LHS, typename RHS>
pow_expr<RealType, DerivativeOrder, LHS, RHS> pow(const expression<RealType, DerivativeOrder, LHS> &lhs,
                                 const expression<RealType, DerivativeOrder, RHS> &rhs)
{
    return pow_expr<RealType, DerivativeOrder, LHS, RHS>(lhs, rhs);
}

template<typename RealType2,
         typename RealType1,
         size_t DerivativeOrder,
         typename ARG,
         typename = typename std::enable_if<!detail::is_expression<RealType2>::value>::type>
expr_pow_float_expr<RealType1, DerivativeOrder, ARG> pow(
    const expression<RealType1, DerivativeOrder, ARG> &arg, const RealType2 &v)
{
    return expr_pow_float_expr<RealType1, DerivativeOrder, ARG>(arg, static_cast<RealType1>(v));
};

template<typename RealType, size_t DerivativeOrder, typename ARG>
float_pow_expr_expr<RealType, DerivativeOrder, ARG> pow(const RealType &v, const expression<RealType,DerivativeOrder, ARG> &arg)
{
    return float_pow_expr_expr<RealType, DerivativeOrder, ARG>(arg, v);
};

template<typename RealType, size_t DerivativeOrder, typename ARG>
log_expr<RealType, DerivativeOrder, ARG> log(const expression<RealType,DerivativeOrder, ARG> &arg)
{
    return log_expr<RealType, DerivativeOrder, ARG>(arg, static_cast<RealType>(0.0));
};

template<typename RealType, size_t DerivativeOrder, typename ARG>
sqrt_expr<RealType, DerivativeOrder, ARG> sqrt(const expression<RealType,DerivativeOrder, ARG> &arg)
{
    return sqrt_expr<RealType, DerivativeOrder, ARG>(arg, static_cast<RealType>(0.0));
};

template<typename RealType, size_t DerivativeOrder, typename ARG>
auto frexp(const expression<RealType, DerivativeOrder, ARG> &arg, int *i)
{
    BOOST_MATH_STD_USING
    frexp(arg.evaluate(), i);
    return arg / pow(static_cast<RealType>(2.0), *i);
}

template<typename RealType, size_t DerivativeOrder, typename ARG>
auto ldexp(const expression<RealType, DerivativeOrder, ARG> &arg, const int &i)
{
    BOOST_MATH_STD_USING
    return arg * pow(static_cast<RealType>(2.0), i);
}

template<typename RealType, size_t DerivativeOrder, typename ARG>
cos_expr<RealType, DerivativeOrder, ARG> cos(const expression<RealType,DerivativeOrder, ARG> &arg)
{
    return cos_expr<RealType, DerivativeOrder, ARG>(arg, static_cast<RealType>(0.0));
};

template<typename RealType, size_t DerivativeOrder, typename ARG>
sin_expr<RealType, DerivativeOrder, ARG> sin(const expression<RealType,DerivativeOrder, ARG> &arg)
{
    return sin_expr<RealType, DerivativeOrder, ARG>(arg, static_cast<RealType>(0.0));
};

template<typename RealType, size_t DerivativeOrder, typename ARG>
tan_expr<RealType, DerivativeOrder, ARG> tan(const expression<RealType,DerivativeOrder, ARG> &arg)
{
    return tan_expr<RealType, DerivativeOrder, ARG>(arg, static_cast<RealType>(0.0));
};

template<typename RealType, size_t DerivativeOrder, typename ARG>
acos_expr<RealType, DerivativeOrder, ARG> acos(const expression<RealType,DerivativeOrder, ARG> &arg)
{
    return acos_expr<RealType, DerivativeOrder, ARG>(arg, static_cast<RealType>(0.0));
};

template<typename RealType, size_t DerivativeOrder, typename ARG>
asin_expr<RealType, DerivativeOrder, ARG> asin(const expression<RealType,DerivativeOrder, ARG> &arg)
{
    return asin_expr<RealType, DerivativeOrder, ARG>(arg, static_cast<RealType>(0.0));
};

template<typename RealType, size_t DerivativeOrder, typename ARG>
atan_expr<RealType, DerivativeOrder, ARG> atan(const expression<RealType,DerivativeOrder, ARG> &arg)
{
    return atan_expr<RealType, DerivativeOrder, ARG>(arg, static_cast<RealType>(0.0));
};

template<typename RealType, size_t DerivativeOrder, typename LHS, typename RHS>
atan2_expr<RealType, DerivativeOrder, LHS, RHS> atan2(const expression<RealType, DerivativeOrder, LHS> &lhs,
                                     const expression<RealType, DerivativeOrder, RHS> &rhs)
{
    return atan2_expr<RealType, DerivativeOrder, LHS, RHS>(lhs, rhs);
}

template<typename RealType, size_t DerivativeOrder, typename ARG>
atan2_right_float_expr<RealType, DerivativeOrder, ARG> atan2(const expression<RealType,DerivativeOrder, ARG> &arg, const RealType &v)
{
    return atan2_right_float_expr<RealType, DerivativeOrder, ARG>(arg, v);
};

template<typename RealType, size_t DerivativeOrder, typename ARG>
atan2_left_float_expr<RealType, DerivativeOrder, ARG> atan2(const RealType &v, const expression<RealType,DerivativeOrder, ARG> &arg)
{
    return atan2_left_float_expr<RealType, DerivativeOrder, ARG>(arg, v);
};

template<typename RealType, size_t DerivativeOrder, typename ARG>
trunc_expr<RealType, DerivativeOrder, ARG> trunc(const expression<RealType,DerivativeOrder, ARG> &arg)
{
    return trunc_expr<RealType, DerivativeOrder, ARG>(arg, static_cast<RealType>(0.0));
}

template<typename RealType, size_t DerivativeOrder, typename LHS, typename RHS>
auto fmod(const expression<RealType, DerivativeOrder, LHS> &lhs, const expression<RealType, DerivativeOrder, RHS> &rhs)
{
    return fmod_expr<RealType, DerivativeOrder, LHS, RHS>(lhs, rhs);
}

template<typename RealType, size_t DerivativeOrder, typename ARG>
auto fmod(const expression<RealType, DerivativeOrder, ARG> &lhs, const RealType rhs)
{
    return fmod_right_float_expr<RealType, DerivativeOrder, ARG>(lhs, rhs);
}

template<typename RealType, size_t DerivativeOrder, typename ARG>
auto fmod(const RealType lhs, const expression<RealType, DerivativeOrder, ARG> &rhs)
{
    return fmod_left_float_expr<RealType, DerivativeOrder, ARG>(rhs, lhs);
}

template<typename RealType, size_t DerivativeOrder, typename ARG>
round_expr<RealType, DerivativeOrder, ARG> round(const expression<RealType,DerivativeOrder, ARG> &arg)
{
    return round_expr<RealType, DerivativeOrder, ARG>(arg, static_cast<RealType>(0.0));
}

template<typename RealType, size_t DerivativeOrder, typename ARG>
int iround(const expression<RealType, DerivativeOrder, ARG> &arg)
{
    rvar<RealType, DerivativeOrder> tmp = arg.evaluate();
    return iround(tmp.item());
}
template<typename RealType, size_t DerivativeOrder, typename ARG>
long lround(const expression<RealType, DerivativeOrder, ARG> &arg)
{
    BOOST_MATH_STD_USING
    rvar<RealType, DerivativeOrder> tmp = arg.evaluate();
    return lround(tmp.item());
}

template<typename RealType, size_t DerivativeOrder, typename ARG>
long long llround(const expression<RealType, DerivativeOrder, ARG> &arg)
{
    rvar<RealType, DerivativeOrder> tmp = arg.evaluate();
    return llround(tmp.item());
}

template<typename RealType, size_t DerivativeOrder, typename ARG>
int itrunc(const expression<RealType, DerivativeOrder, ARG> &arg)
{
    rvar<RealType, DerivativeOrder> tmp = arg.evaluate();
    return itrunc(tmp.item());
}

template<typename RealType, size_t DerivativeOrder, typename ARG>
long ltrunc(const expression<RealType, DerivativeOrder, ARG> &arg)
{
    rvar<RealType, DerivativeOrder> tmp = arg.evaluate();
    return ltrunc(tmp.item());
}

template<typename RealType, size_t DerivativeOrder, typename ARG>
long long lltrunc(const expression<RealType, DerivativeOrder, ARG> &arg)
{
    rvar<RealType, DerivativeOrder> tmp = arg.evaluate();
    return lltrunc(tmp.item());
}

template<typename RealType, size_t DerivativeOrder, typename ARG>
sinh_expr<RealType, DerivativeOrder, ARG> sinh(const expression<RealType,DerivativeOrder, ARG> &arg)
{
    return sinh_expr<RealType, DerivativeOrder, ARG>(arg, static_cast<RealType>(0.0));
}

template<typename RealType, size_t DerivativeOrder, typename ARG>
cosh_expr<RealType, DerivativeOrder, ARG> cosh(const expression<RealType,DerivativeOrder, ARG> &arg)
{
    return cosh_expr<RealType, DerivativeOrder, ARG>(arg, static_cast<RealType>(0.0));
}

template<typename RealType, size_t DerivativeOrder, typename ARG>
tanh_expr<RealType, DerivativeOrder, ARG> tanh(const expression<RealType,DerivativeOrder, ARG> &arg)
{
    return tanh_expr<RealType, DerivativeOrder, ARG>(arg, static_cast<RealType>(0.0));
}

template<typename RealType, size_t DerivativeOrder, typename ARG>
log10_expr<RealType, DerivativeOrder, ARG> log10(const expression<RealType,DerivativeOrder, ARG> &arg)
{
    return log10_expr<RealType, DerivativeOrder, ARG>(arg, static_cast<RealType>(0.0));
}
template<typename RealType, size_t DerivativeOrder, typename ARG>
asinh_expr<RealType, DerivativeOrder, ARG> asinh(const expression<RealType,DerivativeOrder, ARG> &arg)
{
    return asinh_expr<RealType, DerivativeOrder, ARG>(arg, static_cast<RealType>(0.0));
}
template<typename RealType, size_t DerivativeOrder, typename ARG>
acosh_expr<RealType, DerivativeOrder, ARG> acosh(const expression<RealType,DerivativeOrder, ARG> &arg)
{
    return acosh_expr<RealType, DerivativeOrder, ARG>(arg, static_cast<RealType>(0.0));
}
template<typename RealType, size_t DerivativeOrder, typename ARG>
atanh_expr<RealType, DerivativeOrder, ARG> atanh(const expression<RealType,DerivativeOrder, ARG> &arg)
{
    return atanh_expr<RealType, DerivativeOrder, ARG>(arg, static_cast<RealType>(0.0));
}
} // namespace reverse_mode
} // namespace differentiation
} // namespace math
} // namespace boost
#endif
