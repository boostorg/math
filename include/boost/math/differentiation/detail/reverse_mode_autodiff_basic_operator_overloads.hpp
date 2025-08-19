//           Copyright Maksym Zhelyenzyakov 2025-2026.
// Distributed under the Boost Software License, Version 1.0.
//      (See accompanying file LICENSE_1_0.txt or copy at
//           https://www.boost.org/LICENSE_1_0.txt)
#ifndef REVERSE_MODE_AUTODOFF_BASIC_OPERATOR_OVERLOADS_HPP
#define REVERSE_MODE_AUTODOFF_BASIC_OPERATOR_OVERLOADS_HPP

#include <boost/math/differentiation/detail/reverse_mode_autodiff_expression_template_base.hpp>

namespace boost {
namespace math {
namespace differentiation {
namespace reverse_mode {
/****************************************************************************************************************/
template<typename T, size_t order, typename LHS, typename RHS>
struct add_expr
    : public abstract_binary_expression<T, order, LHS, RHS, add_expr<T, order, LHS, RHS>>
{
    /* @brief addition
   * rvar+rvar
   * */
    using inner_t    = rvar_t<T, order - 1>;
    // Explicitly define constructor to forward to base class
    explicit add_expr(const expression<T, order, LHS> &left_hand_expr,
                      const expression<T, order, RHS> &right_hand_expr)
        : abstract_binary_expression<T, order, LHS, RHS, add_expr<T, order, LHS, RHS>>(
              left_hand_expr, right_hand_expr)
    {}

    inner_t              evaluate() const { return this->lhs.evaluate() + this->rhs.evaluate(); }
    static const inner_t left_derivative(const inner_t & /*l*/,
                                         const inner_t & /*r*/,
                                         const inner_t & /*v*/)
    {
        return inner_t(1.0);
    }
    static const inner_t right_derivative(const inner_t & /*l*/,
                                          const inner_t & /*r*/,
                                          const inner_t & /*v*/)
    {
        return inner_t(1.0);
    }
};
template<typename T, size_t order, typename ARG>
struct add_const_expr
    : public abstract_unary_expression<T, order, ARG, add_const_expr<T, order, ARG>>
{
    /* @brief
   * rvar+float or float+rvar
   * */
    using inner_t    = rvar_t<T, order - 1>;
    explicit add_const_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, add_const_expr<T, order, ARG>>(arg_expr, v){};
    inner_t              evaluate() const { return this->arg.evaluate() + inner_t(this->constant); }
    static const inner_t derivative(const inner_t & /*argv*/,
                                    const inner_t & /*v*/,
                                    const T & /*constant*/)
    {
        return inner_t(1.0);
    }
};
/****************************************************************************************************************/
template<typename T, size_t order, typename LHS, typename RHS>
struct mult_expr
    : public abstract_binary_expression<T, order, LHS, RHS, mult_expr<T, order, LHS, RHS>>
{
    /* @brief multiplication
   * rvar * rvar
   * */
    using inner_t    = rvar_t<T, order - 1>;
    explicit mult_expr(const expression<T, order, LHS> &left_hand_expr,
                       const expression<T, order, RHS> &right_hand_expr)
        : abstract_binary_expression<T, order, LHS, RHS, mult_expr<T, order, LHS, RHS>>(
              left_hand_expr, right_hand_expr)
    {}

    inner_t              evaluate() const { return this->lhs.evaluate() * this->rhs.evaluate(); };
    static const inner_t left_derivative(const inner_t & /*l*/,
                                         const inner_t &r,
                                         const inner_t & /*v*/)
    {
        return r;
    };
    static const inner_t right_derivative(const inner_t &l,
                                          const inner_t & /*r*/,
                                          const inner_t & /*v*/)
    {
        return l;
    };
};
template<typename T, size_t order, typename ARG>
struct mult_const_expr
    : public abstract_unary_expression<T, order, ARG, mult_const_expr<T, order, ARG>>
{
    /* @brief
   * rvar+float or float+rvar
   * */
    using inner_t    = rvar_t<T, order - 1>;

    explicit mult_const_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, mult_const_expr<T, order, ARG>>(arg_expr, v){};

    inner_t              evaluate() const { return this->arg.evaluate() * inner_t(this->constant); }
    static const inner_t derivative(const inner_t & /*argv*/,
                                    const inner_t & /*v*/,
                                    const T &constant)
    {
        return inner_t(constant);
    }
};
/****************************************************************************************************************/
template<typename T, size_t order, typename LHS, typename RHS>
struct sub_expr
    : public abstract_binary_expression<T, order, LHS, RHS, sub_expr<T, order, LHS, RHS>>
{
    /* @brief addition
   * rvar-rvar
   * */
    using inner_t    = rvar_t<T, order - 1>;
    // Explicitly define constructor to forward to base class
    explicit sub_expr(const expression<T, order, LHS> &left_hand_expr,
                      const expression<T, order, RHS> &right_hand_expr)
        : abstract_binary_expression<T, order, LHS, RHS, sub_expr<T, order, LHS, RHS>>(
              left_hand_expr, right_hand_expr)
    {}

    inner_t              evaluate() const { return this->lhs.evaluate() - this->rhs.evaluate(); }
    static const inner_t left_derivative(const inner_t & /*l*/,
                                         const inner_t & /*r*/,
                                         const inner_t & /*v*/)
    {
        return inner_t(1.0);
    }
    static const inner_t right_derivative(const inner_t & /*l*/,
                                          const inner_t & /*r*/,
                                          const inner_t & /*v*/)
    {
        return inner_t(-1.0);
    }
};

/****************************************************************************************************************/
template<typename T, size_t order, typename LHS, typename RHS>
struct div_expr
    : public abstract_binary_expression<T, order, LHS, RHS, div_expr<T, order, LHS, RHS>>
{
    /* @brief multiplication
   * rvar / rvar
   * */
    using inner_t    = rvar_t<T, order - 1>;
    // Explicitly define constructor to forward to base class
    explicit div_expr(const expression<T, order, LHS> &left_hand_expr,
                      const expression<T, order, RHS> &right_hand_expr)
        : abstract_binary_expression<T, order, LHS, RHS, div_expr<T, order, LHS, RHS>>(
              left_hand_expr, right_hand_expr)
    {}

    inner_t              evaluate() const { return this->lhs.evaluate() / this->rhs.evaluate(); };
    static const inner_t left_derivative(const inner_t & /*l*/,
                                         const inner_t &r,
                                         const inner_t & /*v*/)
    {
        return static_cast<T>(1.0) / r;
    };
    static const inner_t right_derivative(const inner_t &l, const inner_t &r, const inner_t & /*v*/)
    {
        return -l / (r * r);
    };
};
template<typename T, size_t order, typename ARG>
struct div_by_const_expr
    : public abstract_unary_expression<T, order, ARG, div_by_const_expr<T, order, ARG>>
{
    /* @brief
   * rvar/float
   * */
    using inner_t    = rvar_t<T, order - 1>;

    explicit div_by_const_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, div_by_const_expr<T, order, ARG>>(arg_expr, v){};

    inner_t              evaluate() const { return this->arg.evaluate() / inner_t(this->constant); }
    static const inner_t derivative(const inner_t & /*argv*/,
                                    const inner_t & /*v*/,
                                    const T &constant)
    {
        return inner_t(1.0 / constant);
    }
};

template<typename T, size_t order, typename ARG>
struct const_div_by_expr
    : public abstract_unary_expression<T, order, ARG, const_div_by_expr<T, order, ARG>>
{
    /** @brief
    * float/rvar
    * */
    using inner_t    = rvar_t<T, order - 1>;

    explicit const_div_by_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, const_div_by_expr<T, order, ARG>>(arg_expr, v){};

    inner_t              evaluate() const { return inner_t(this->constant) / this->arg.evaluate(); }
    static const inner_t derivative(const inner_t &argv, const inner_t & /*v*/, const T &constant)
    {
        return -inner_t{constant} / (argv * argv);
    }
};
/****************************************************************************************************************/

template<typename T, size_t order, typename LHS, typename RHS>
mult_expr<T, order, LHS, RHS> operator*(const expression<T, order, LHS> &lhs,
                                        const expression<T, order, RHS> &rhs)
{
    return mult_expr<T, order, LHS, RHS>(lhs, rhs);
}

/** @brief type promotion is handled by casting the numeric type to
 *  the type inside expression. This is to avoid converting the
 *  entire tape in case you have something like double * rvar<float>
 *  */
template<typename U,
         typename T,
         size_t order,
         typename ARG,
         typename = typename std::enable_if<!detail::is_expression<U>::value>::type>
mult_const_expr<T, order, ARG> operator*(const expression<T, order, ARG> &arg, const U &v)
{
    return mult_const_expr<T, order, ARG>(arg, static_cast<T>(v));
}
template<typename U,
         typename T,
         size_t order,
         typename ARG,
         typename = typename std::enable_if<!detail::is_expression<U>::value>::type>
mult_const_expr<T, order, ARG> operator*(const U &v, const expression<T, order, ARG> &arg)
{
    return mult_const_expr<T, order, ARG>(arg, static_cast<T>(v));
}
/****************************************************************************************************************/
/* + */
template<typename T, size_t order, typename LHS, typename RHS>
add_expr<T, order, LHS, RHS> operator+(const expression<T, order, LHS> &lhs,
                                       const expression<T, order, RHS> &rhs)
{
    return add_expr<T, order, LHS, RHS>(lhs, rhs);
}
template<typename U,
         typename T,
         size_t order,
         typename ARG,
         typename = typename std::enable_if<!detail::is_expression<U>::value>::type>
add_const_expr<T, order, ARG> operator+(const expression<T, order, ARG> &arg, const U &v)
{
    return add_const_expr<T, order, ARG>(arg, static_cast<T>(v));
}
template<typename U,
         typename T,
         size_t order,
         typename ARG,
         typename = typename std::enable_if<!detail::is_expression<U>::value>::type>
add_const_expr<T, order, ARG> operator+(const U &v, const expression<T, order, ARG> &arg)
{
    return add_const_expr<T, order, ARG>(arg, static_cast<T>(v));
}
/****************************************************************************************************************/
/* - overload */
/** @brief
 *  negation (-1.0*rvar) */
template<typename T, size_t order, typename ARG>
mult_const_expr<T, order, ARG> operator-(const expression<T, order, ARG> &arg)
{
    return mult_const_expr<T, order, ARG>(arg, static_cast<T>(-1.0));
}

/** @brief
 *  subtraction rvar-rvar */
template<typename T, size_t order, typename LHS, typename RHS>
sub_expr<T, order, LHS, RHS> operator-(const expression<T, order, LHS> &lhs,
                                       const expression<T, order, RHS> &rhs)
{
    return sub_expr<T, order, LHS, RHS>(lhs, rhs);
}

/** @brief
 *  subtraction float - rvar */
template<typename U,
         typename T,
         size_t order,
         typename ARG,
         typename = typename std::enable_if<!detail::is_expression<U>::value>::type>
add_const_expr<T, order, ARG> operator-(const expression<T, order, ARG> &arg, const U &v)
{
    /* rvar - float = rvar + (-float) */
    return add_const_expr<T, order, ARG>(arg, static_cast<T>(-v));
}

/** @brief
 *   subtraction float - rvar
 *  @return add_expr<neg_expr<ARG>>
 */
template<typename U,
         typename T,
         size_t order,
         typename ARG,
         typename = typename std::enable_if<!detail::is_expression<U>::value>::type>
auto operator-(const U &v, const expression<T, order, ARG> &arg)
{
    auto neg = -arg;
    return neg + static_cast<T>(v);
}
/****************************************************************************************************************/
/* / */
template<typename T, size_t order, typename LHS, typename RHS>
div_expr<T, order, LHS, RHS> operator/(const expression<T, order, LHS> &lhs,
                                       const expression<T, order, RHS> &rhs)
{
    return div_expr<T, order, LHS, RHS>(lhs, rhs);
}

template<typename U,
         typename T,
         size_t order,
         typename ARG,
         typename = typename std::enable_if<!detail::is_expression<U>::value>::type>
const_div_by_expr<T, order, ARG> operator/(const U &v, const expression<T, order, ARG> &arg)
{
    return const_div_by_expr<T, order, ARG>(arg, static_cast<T>(v));
}

template<typename U,
         typename T,
         size_t order,
         typename ARG,
         typename = typename std::enable_if<!detail::is_expression<U>::value>::type>
div_by_const_expr<T, order, ARG> operator/(const expression<T, order, ARG> &arg, const U &v)
{
    return div_by_const_expr<T, order, ARG>(arg, static_cast<T>(v));
}

} // namespace reverse_mode
} // namespace differentiation
} // namespace math
} // namespace boost

#endif
