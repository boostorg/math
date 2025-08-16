//           Copyright Maksym Zhelyenzyakov 2025-2026.
// Distributed under the Boost Software License, Version 1.0.
//      (See accompanying file LICENSE_1_0.txt or copy at
//           https://www.boost.org/LICENSE_1_0.txt)

#ifndef REVERSE_MODE_AUTODIFF_COMPARISON_OPERATOR_OVERLOADS_HPP
#define REVERSE_MODE_AUTODIFF_COMPARISON_OPERATOR_OVERLOADS_HPP
#include <boost/math/differentiation/detail/reverse_mode_autodiff_expression_template_base.hpp>
namespace boost {
namespace math {
namespace differentiation {
namespace reverse_mode {
template<typename T, size_t order_1, size_t order_2, class E, class F>
bool operator==(const expression<T, order_1, E> &lhs, const expression<T, order_2, F> &rhs)
{
    return lhs.evaluate() == rhs.evaluate();
}

template<typename U,
         typename T,
         size_t order,
         class E,
         typename = typename std::enable_if<!detail::is_expression<U>::value>::type>
bool operator==(const expression<T, order, E> &lhs, const U &rhs)
{
    return lhs.evaluate() == static_cast<T>(rhs);
}
template<typename U,
         typename T,
         size_t order,
         class E,
         typename = typename std::enable_if<!detail::is_expression<U>::value>::type>
bool operator==(const U &lhs, const expression<T, order, E> &rhs)
{
    return lhs == rhs.evaluate();
}

template<typename T, size_t order_1, size_t order_2, class E, class F>
bool operator!=(const expression<T, order_1, E> &lhs, const expression<T, order_2, F> &rhs)
{
    return lhs.evaluate() != rhs.evaluate();
}
template<typename U,
         typename T,
         size_t order,
         class E,
         typename = typename std::enable_if<!detail::is_expression<U>::value>::type>
bool operator!=(const expression<T, order, E> &lhs, const U &rhs)
{
    return lhs.evaluate() != rhs;
}
template<typename U,
         typename T,
         size_t order,
         class E,
         typename = typename std::enable_if<!detail::is_expression<U>::value>::type>
bool operator!=(const U &lhs, const expression<T, order, E> &rhs)
{
    return lhs != rhs.evaluate();
}

template<typename T, size_t order_1, size_t order_2, class E, class F>
bool operator<(const expression<T, order_1, E> &lhs, const expression<T, order_2, F> &rhs)
{
    return lhs.evaluate() < rhs.evaluate();
}
template<typename U,
         typename T,
         size_t order,
         class E,
         typename = typename std::enable_if<!detail::is_expression<U>::value>::type>
bool operator<(const expression<T, order, E> &lhs, const U &rhs)
{
    return lhs.evaluate() < rhs;
}
template<typename U,
         typename T,
         size_t order,
         class E,
         typename = typename std::enable_if<!detail::is_expression<U>::value>::type>
bool operator<(const U &lhs, const expression<T, order, E> &rhs)
{
    return lhs < rhs.evaluate();
}

template<typename T, size_t order_1, size_t order_2, class E, class F>
bool operator>(const expression<T, order_1, E> &lhs, const expression<T, order_2, F> &rhs)
{
    return lhs.evaluate() > rhs.evaluate();
}
template<typename U,
         typename T,
         size_t order,
         class E,
         typename = typename std::enable_if<!detail::is_expression<U>::value>::type>
bool operator>(const expression<T, order, E> &lhs, const U &rhs)
{
    return lhs.evaluate() > rhs;
}
template<typename U,
         typename T,
         size_t order,
         class E,
         typename = typename std::enable_if<!detail::is_expression<U>::value>::type>
bool operator>(const U &lhs, const expression<T, order, E> &rhs)
{
    return lhs > rhs.evaluate();
}

template<typename T, size_t order_1, size_t order_2, class E, class F>
bool operator<=(const expression<T, order_1, E> &lhs, const expression<T, order_2, F> &rhs)
{
    return lhs.evaluate() <= rhs.evaluate();
}
template<typename U,
         typename T,
         size_t order,
         class E,
         typename = typename std::enable_if<!detail::is_expression<U>::value>::type>
bool operator<=(const expression<T, order, E> &lhs, const U &rhs)
{
    return lhs.evaluate() <= rhs;
}
template<typename U,
         typename T,
         size_t order,
         class E,
         typename = typename std::enable_if<!detail::is_expression<U>::value>::type>
bool operator<=(const U &lhs, const expression<T, order, E> &rhs)
{
    return lhs <= rhs.evaluate();
}

template<typename T, size_t order_1, size_t order_2, class E, class F>
bool operator>=(const expression<T, order_1, E> &lhs, const expression<T, order_2, F> &rhs)
{
    return lhs.evaluate() >= rhs.evaluate();
}
template<typename U,
         typename T,
         size_t order,
         class E,
         typename = typename std::enable_if<!detail::is_expression<U>::value>::type>
bool operator>=(const expression<T, order, E> &lhs, const U &rhs)
{
    return lhs.evaluate() >= rhs;
}
template<typename U,
         typename T,
         size_t order,
         class E,
         typename = typename std::enable_if<!detail::is_expression<U>::value>::type>
bool operator>=(const U &lhs, const expression<T, order, E> &rhs)
{
    return lhs >= rhs.evaluate();
}
} // namespace reverse_mode
} // namespace differentiation
} // namespace math
} // namespace boost
#endif // REVERSE_MODE_AUTODIFF_COMPARISON_OPERATOR_OVERLOADS_HPP
