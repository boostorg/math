//           Copyright Maksym Zhelyenzyakov 2025-2026.
// Distributed under the Boost Software License, Version 1.0.
//      (See accompanying file LICENSE_1_0.txt or copy at
//           https://www.boost.org/LICENSE_1_0.txt)

#ifndef REVERSE_MODE_AUTODIFF_EXPRESSION_TEMPLATE_BASE_HPP
#define REVERSE_MODE_AUTODIFF_EXPRESSION_TEMPLATE_BASE_HPP
#include <stddef.h>
#include <type_traits>
namespace boost {
namespace math {
namespace differentiation {
namespace reverse_mode {

/* forward declarations for utitlity functions */
struct expression_base
{};

template<typename T, size_t order, class derived_expression>
struct expression;

template<typename T, size_t order>
struct rvar;

template<typename T, size_t order, typename LHS, typename RHS, typename concrete_binary_operation>
struct abstract_binary_expression;

template<typename T, size_t order, typename ARG, typename concrete_unary_operation>
struct abstract_unary_expression;

template<typename T, size_t order>
class gradient_node; // forward declaration for tape

namespace detail {
template<typename...>
using void_t = void;
// Check if T has a 'value_type' alias
template<typename T, typename Enable = void>
struct has_value_type : std::false_type
{};
template<typename T>
struct has_value_type<T, void_t<typename T::value_type>> : std::true_type
{};
template<typename T, typename Enable = void>
struct has_binary_sub_types : std::false_type
{};
template<typename T>
struct has_binary_sub_types<T, void_t<typename T::lhs_type, typename T::rhs_type>> : std::true_type
{};
template<typename T, typename Enable = void>
struct has_unary_sub_type : std::false_type
{};
template<typename T>
struct has_unary_sub_type<T, void_t<typename T::arg_type>> : std::true_type
{};

template<typename T, size_t order, typename Enable = void>
struct count_rvar_impl
{
    static constexpr std::size_t value = 0;
};
template<typename U, size_t order>
struct count_rvar_impl<rvar<U, order>, order>
{
    static constexpr std::size_t value = 1;
};

template<typename T, std::size_t order>
struct count_rvar_impl<T,
                       order,
                       std::enable_if_t<has_binary_sub_types<T>::value
                                        && !std::is_same<T, rvar<typename T::value_type, order>>::value
                                        && !has_unary_sub_type<T>::value>>
{
    static constexpr std::size_t value = count_rvar_impl<typename T::lhs_type, order>::value
                                         + count_rvar_impl<typename T::rhs_type, order>::value;
};

template<typename T, size_t order>
struct count_rvar_impl<
    T,
    order,
    typename std::enable_if_t<has_unary_sub_type<T>::value
                              && !std::is_same<T, rvar<typename T::value_type, order>>::value
                              && !has_binary_sub_types<T>::value>>
{
    static constexpr std::size_t value = count_rvar_impl<typename T::arg_type, order>::value;
};
template<typename T, size_t order>
constexpr std::size_t count_rvars = detail::count_rvar_impl<T, order>::value;

template<typename T>
struct is_expression : std::is_base_of<expression_base, typename std::decay<T>::type>
{};

template<typename T, size_t N>
struct rvar_type_impl
{
    using type = rvar<T, N>;
};

template<typename T>
struct rvar_type_impl<T, 0>
{
    using type = T;
};

} // namespace detail

template<typename T, size_t N>
using rvar_t = typename detail::rvar_type_impl<T, N>::type;

template<typename T, size_t order, class derived_expression>
struct expression : expression_base
{
    /* @brief
   * base expression class
   * */

    using value_type                = T;
    static constexpr size_t order_v = order;
    using derived_type              = derived_expression;

    static constexpr size_t num_literals = 0;
    using inner_t                        = rvar_t<T, order - 1>;
    inner_t evaluate() const { return static_cast<const derived_expression *>(this)->evaluate(); }

    template<size_t arg_index>
    void propagatex(gradient_node<T, order> *node, inner_t adj) const
    {
        return static_cast<const derived_expression *>(this)->template propagatex<arg_index>(node,
                                                                                             adj);
    };
};

template<typename T, size_t order, typename LHS, typename RHS, typename concrete_binary_operation>
struct abstract_binary_expression
    : public expression<T,
                        order,
                        abstract_binary_expression<T, order, LHS, RHS, concrete_binary_operation>>
{
    using lhs_type   = LHS;
    using rhs_type   = RHS;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;
    const lhs_type lhs;
    const rhs_type rhs;

    explicit abstract_binary_expression(const expression<T, order, LHS> &left_hand_expr,
                                        const expression<T, order, RHS> &right_hand_expr)
        : lhs(static_cast<const LHS &>(left_hand_expr))
        , rhs(static_cast<const RHS &>(right_hand_expr)){};

    inner_t evaluate() const
    {
        return static_cast<const concrete_binary_operation *>(this)->evaluate();
    };

    template<size_t arg_index>
    void propagatex(gradient_node<T, order> *node, inner_t adj) const
    {
        inner_t lv        = lhs.evaluate();
        inner_t rv        = rhs.evaluate();
        inner_t v         = evaluate();
        inner_t partial_l = concrete_binary_operation::left_derivative(lv, rv, v);
        inner_t partial_r = concrete_binary_operation::right_derivative(lv, rv, v);

        constexpr size_t num_lhs_args = detail::count_rvars<LHS, order>;
        constexpr size_t num_rhs_args = detail::count_rvars<RHS, order>;

        propagate_lhs<num_lhs_args, arg_index>(node, adj * partial_l);
        propagate_rhs<num_rhs_args, arg_index + num_lhs_args>(node, adj * partial_r);
    }

private:
    /* everything here just emulates c++17 if constexpr */

    template<std::size_t num_args,
             std::size_t arg_index_,
             typename std::enable_if<(num_args > 0), int>::type = 0>
    void propagate_lhs(gradient_node<T, order> *node, inner_t adj) const
    {
        lhs.template propagatex<arg_index_>(node, adj);
    }

    template<std::size_t num_args,
             std::size_t arg_index_,
             typename std::enable_if<(num_args == 0), int>::type = 0>
    void propagate_lhs(gradient_node<T, order> *, inner_t) const
    {}

    template<std::size_t num_args,
             std::size_t arg_index_,
             typename std::enable_if<(num_args > 0), int>::type = 0>
    void propagate_rhs(gradient_node<T, order> *node, inner_t adj) const
    {
        rhs.template propagatex<arg_index_>(node, adj);
    }

    template<std::size_t num_args,
             std::size_t arg_index_,
             typename std::enable_if<(num_args == 0), int>::type = 0>
    void propagate_rhs(gradient_node<T, order> *, inner_t) const
    {}
};
template<typename T, size_t order, typename ARG, typename concrete_unary_operation>
struct abstract_unary_expression
    : public expression<T, order, abstract_unary_expression<T, order, ARG, concrete_unary_operation>>
{
    using arg_type   = ARG;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;
    const arg_type arg;
    const T        constant;
    explicit abstract_unary_expression(const expression<T, order, ARG> &arg_expr, const T &constant)
        : arg(static_cast<const ARG &>(arg_expr))
        , constant(constant){};
    // explicit abstract_unary_expression(const ARG& arg_expr, const T& constant)
    //    : arg(arg_expr)
    //    , constant(constant)
    //{}
    inner_t evaluate() const
    {
        return static_cast<const concrete_unary_operation *>(this)->evaluate();
    };

    template<size_t arg_index>
    void propagatex(gradient_node<T, order> *node, inner_t adj) const
    {
        inner_t argv        = arg.evaluate();
        inner_t v           = evaluate();
        inner_t partial_arg = concrete_unary_operation::derivative(argv, v, constant);

        arg.template propagatex<arg_index>(node, adj * partial_arg);
    }
};
} // namespace reverse_mode
} // namespace differentiation
} // namespace math
} // namespace boost

#endif // REVERSE_MODE_AUTODIFF_EXPRESSION_TEMPLATE_BASE_HPP
