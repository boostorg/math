//           Copyright Maksym Zhelyenzyakov 2025-2026.
// Distributed under the Boost Software License, Version 1.0.
//      (See accompanying file LICENSE_1_0.txt or copy at
//           https://www.boost.org/LICENSE_1_0.txt)
#ifndef BOOST_MATH_DIFFERENTIATION_AUTODIFF_HPP
#define BOOST_MATH_DIFFERENTIATION_AUTODIFF_HPP

#include <boost/cstdfloat.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/differentiation/detail/reverse_mode_autodiff_basic_operator_overloads.hpp>
#include <boost/math/differentiation/detail/reverse_mode_autodiff_comparison_operator_overloads.hpp>
#include <boost/math/differentiation/detail/reverse_mode_autodiff_erf_overloads.hpp>
#include <boost/math/differentiation/detail/reverse_mode_autodiff_expression_template_base.hpp>
#include <boost/math/differentiation/detail/reverse_mode_autodiff_memory_management.hpp>
#include <boost/math/differentiation/detail/reverse_mode_autodiff_stl_overloads.hpp>
#include <boost/math/special_functions/acosh.hpp>
#include <boost/math/special_functions/asinh.hpp>
#include <boost/math/special_functions/atanh.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/special_functions/lambert_w.hpp>
#include <boost/math/special_functions/polygamma.hpp>
#include <boost/math/special_functions/round.hpp>
#include <boost/math/special_functions/trunc.hpp>
#include <boost/math/tools/config.hpp>
#include <boost/math/tools/promotion.hpp>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <type_traits>
#include <vector>
#define BOOST_MATH_BUFFER_SIZE 65536

namespace boost {
namespace math {
namespace differentiation {
namespace reverse_mode {

/* forward declarations for utitlity functions */
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
// manages nodes in computational graph
template<typename T, size_t order, size_t buffer_size = BOOST_MATH_BUFFER_SIZE>
class gradient_tape
{
    /** @brief tape (graph) management class for autodiff
   *  holds all the data structures for autodiff */
private:
    /* type decays to order - 1 to support higher order derivatives */
    using inner_t = rvar_t<T, order - 1>;
    /* adjoints are the overall derivative, and derivatives are the "local"
   * derivative */
    detail::flat_linear_allocator<inner_t, buffer_size>                   adjoints_;
    detail::flat_linear_allocator<inner_t, buffer_size>                   derivatives_;
    detail::flat_linear_allocator<gradient_node<T, order>, buffer_size>   gradient_nodes_;
    detail::flat_linear_allocator<gradient_node<T, order> *, buffer_size> argument_nodes_;

    // compile time check if emplace_back calls on zero
    template<size_t n>
    gradient_node<T, order> *fill_node_at_compile_time(std::true_type,
                                                       gradient_node<T, order> *node_ptr)
    {
        node_ptr->derivatives_    = derivatives_.template emplace_back_n<n>();
        node_ptr->argument_nodes_ = argument_nodes_.template emplace_back_n<n>();
        return node_ptr;
    }

    template<size_t n>
    gradient_node<T, order> *fill_node_at_compile_time(std::false_type,
                                                       gradient_node<T, order> *node_ptr)
    {
        node_ptr->derivatives_       = nullptr;
        node_ptr->argument_adjoints_ = nullptr;
        node_ptr->argument_nodes_    = nullptr;
        return node_ptr;
    }

public:
    /* gradient node stores iterators to its data memebers
   * (adjoint/derivative/arguments) so that in case flat linear allocator
   * reaches its block boundary and needs more memory for that node, the
   * iterator can be invoked to access it */
    using adjoint_iterator = typename detail::flat_linear_allocator<inner_t, buffer_size>::iterator;
    using derivatives_iterator =
        typename detail::flat_linear_allocator<inner_t, buffer_size>::iterator;
    using gradient_nodes_iterator =
        typename detail::flat_linear_allocator<gradient_node<T, order>, buffer_size>::iterator;
    using argument_nodes_iterator =
        typename detail::flat_linear_allocator<gradient_node<T, order> *, buffer_size>::iterator;

    gradient_tape() { clear(); };

    gradient_tape(const gradient_tape &)            = delete;
    gradient_tape &operator=(const gradient_tape &) = delete;
    gradient_tape(gradient_tape &&other)            = delete;
    gradient_tape operator=(gradient_tape &&other)  = delete;
    ~gradient_tape() noexcept { clear(); }
    void clear() noexcept
    {
        adjoints_.clear();
        derivatives_.clear();
        gradient_nodes_.clear();
        argument_nodes_.clear();
    }

    // no derivatives or arguments
    gradient_node<T, order> *emplace_leaf_node()
    {
        gradient_node<T, order> *node = &*gradient_nodes_.emplace_back();
        node->adjoint_                = adjoints_.emplace_back();
        node->derivatives_            = derivatives_iterator();    // nullptr;
        node->argument_nodes_         = argument_nodes_iterator(); // nullptr;

        return node;
    };

    // single argument, single derivative
    gradient_node<T, order> *emplace_active_unary_node()
    {
        gradient_node<T, order> *node = &*gradient_nodes_.emplace_back();
        node->n_                      = 1;
        node->adjoint_                = adjoints_.emplace_back();
        node->derivatives_            = derivatives_.emplace_back();

        return node;
    };

    // arbitrary number of arguments/derivatives (compile time)
    template<size_t n>
    gradient_node<T, order> *emplace_active_multi_node()
    {
        gradient_node<T, order> *node = &*gradient_nodes_.emplace_back();
        node->n_                      = n;
        node->adjoint_                = adjoints_.emplace_back();
        // emulate if constexpr
        return fill_node_at_compile_time<n>(std::integral_constant<bool, (n > 0)>{}, node);
    };

    // same as above at runtime
    gradient_node<T, order> *emplace_active_multi_node(size_t n)
    {
        gradient_node<T, order> *node = &*gradient_nodes_.emplace_back();
        node->n_                      = n;
        node->adjoint_                = adjoints_.emplace_back();
        if (n > 0) {
            node->derivatives_    = derivatives_.emplace_back_n(n);
            node->argument_nodes_ = argument_nodes_.emplace_back_n(n);
        }
        return node;
    };
    /* manual reset button for all adjoints */
    void zero_grad()
    {
        const T zero = T(0.0);
        adjoints_.fill(zero);
    }

    // return type is an iterator
    auto begin() { return gradient_nodes_.begin(); }
    auto end() { return gradient_nodes_.end(); }
    auto find(gradient_node<T, order> *node) { return gradient_nodes_.find(node); };
    void add_checkpoint()
    {
        gradient_nodes_.add_checkpoint();
        adjoints_.add_checkpoint();
        derivatives_.add_checkpoint();
        argument_nodes_.add_checkpoint();
    };

    auto last_checkpoint() { return gradient_nodes_.last_checkpoint(); };
    auto first_checkpoint() { return gradient_nodes_.last_checkpoint(); };
    auto checkpoint_at(size_t index) { return gradient_nodes_.get_checkpoint_at(index); };
    void rewind_to_last_checkpoint()
    {
        gradient_nodes_.rewind_to_last_checkpoint();
        adjoints_.rewind_to_last_checkpoint();
        derivatives_.rewind_to_last_checkpoint();
        argument_nodes_.rewind_to_last_checkpoint();
    };
    void rewind_to_checkpoint_at(size_t index) // index is "checkpoint" index. so
                                               // order which checkpoint was set
    {
        gradient_nodes_.rewind_to_checkpoint_at(index);
        adjoints_.rewind_to_checkpoint_at(index);
        derivatives_.rewind_to_checkpoint_at(index);
        argument_nodes_.rewind_to_checkpoint_at(index);
    }

    // rewind to beginning of computational graph
    void rewind()
    {
        gradient_nodes_.rewind();
        adjoints_.rewind();
        derivatives_.rewind();
        argument_nodes_.rewind();
    }

    // random acces
    gradient_node<T, order>       &operator[](size_t i) { return gradient_nodes_[i]; }
    const gradient_node<T, order> &operator[](size_t i) const { return gradient_nodes_[i]; }
};
// class rvar;
template<typename T, size_t order> // no CRTP, just storage
class gradient_node
{
    /*
   * @brief manages adjoints, derivatives, and stores points to argument
   * adjoints pointers to arguments aren't needed here
   * */
public:
    using adjoint_iterator        = typename gradient_tape<T, order>::adjoint_iterator;
    using derivatives_iterator    = typename gradient_tape<T, order>::derivatives_iterator;
    using argument_nodes_iterator = typename gradient_tape<T, order>::argument_nodes_iterator;

private:
    size_t n_;
    using inner_t = rvar_t<T, order - 1>;
    /* these are iterators in case
   * flat linear allocator is at capacity, and needs to allocate a new block of
   * memory. */
    adjoint_iterator        adjoint_;
    derivatives_iterator    derivatives_;
    argument_nodes_iterator argument_nodes_;

public:
    friend class gradient_tape<T, order>;
    friend class rvar<T, order>;

    gradient_node() = default;
    explicit gradient_node(const size_t n)
        : n_(n)
        , adjoint_(nullptr)
        , derivatives_(nullptr)
    {}
    explicit gradient_node(const size_t n, T *adjoint, T *derivatives, rvar<T, order> **arguments)
        : n_(n)
        , adjoint_(adjoint)
        , derivatives_(derivatives)
    {}

    inner_t get_adjoint_v() const { return *adjoint_; }
    inner_t get_derivative_v(size_t arg_id) const { return derivatives_[arg_id]; };
    inner_t get_argument_adjoint_v(size_t arg_id) const
    {
        return *argument_nodes_[arg_id]->adjoint_;
    }

    adjoint_iterator get_adjoint_ptr() { return adjoint_; }
    adjoint_iterator get_adjoint_ptr() const { return adjoint_; };
    void             update_adjoint_v(inner_t value) { *adjoint_ = value; };
    void update_derivative_v(size_t arg_id, inner_t value) { derivatives_[arg_id] = value; };
    void update_argument_adj_v(size_t arg_id, inner_t value)
    {
        argument_nodes_[arg_id]->update_adjoint_v(value);
    };
    void update_argument_ptr_at(size_t arg_id, gradient_node<T, order> *node_ptr)
    {
        argument_nodes_[arg_id] = node_ptr;
    }

    void backward()
    {
        if (!n_) // leaf node
            return;

        using boost::math::differentiation::reverse_mode::fabs;
        using std::fabs;
        if (!adjoint_ || fabs(*adjoint_) < 2 * std::numeric_limits<T>::epsilon())
            return;

        if (!argument_nodes_) // no arguments
            return;

        if (!derivatives_) // no derivatives
            return;

        for (size_t i = 0; i < n_; ++i) {
            auto adjoint          = get_adjoint_v();
            auto derivative       = get_derivative_v(i);
            auto argument_adjoint = get_argument_adjoint_v(i);
            update_argument_adj_v(i, argument_adjoint + derivative * adjoint);
        }
    }
};

/****************************************************************************************************************/
template<typename T, size_t order>
inline gradient_tape<T, order, BOOST_MATH_BUFFER_SIZE> &get_active_tape()
{
    static BOOST_MATH_THREAD_LOCAL gradient_tape<T, order, BOOST_MATH_BUFFER_SIZE> tape;
    return tape;
}

template<typename T, size_t order = 1>
class rvar : public expression<T, order, rvar<T, order>>
{
private:
    using inner_t = rvar_t<T, order - 1>;
    friend class gradient_node<T, order>;
    inner_t                  value_;
    gradient_node<T, order> *node_ = nullptr;
    template<typename, size_t>
    friend class rvar;
    /*****************************************************************************************/
    /**
     * @brief implementation helpers for get_value_at
     */
    template<size_t target_order, size_t current_order>
    struct get_value_at_impl
    {
        static_assert(target_order <= current_order, "Requested depth exceeds variable order.");

        /** @return value_ at rvar_t<T,current_order - 1>
         */
        static auto &get(rvar<T, current_order> &v)
        {
            return get_value_at_impl<target_order, current_order - 1>::get(v.get_value());
        }
        /** @return const value_ at rvar_t<T,current_order - 1>
         */
        static const auto &get(const rvar<T, current_order> &v)
        {
            return get_value_at_impl<target_order, current_order - 1>::get(v.get_value());
        }
    };

    /** @brief base case specialization for target_order == current order
     */
    template<size_t target_order>
    struct get_value_at_impl<target_order, target_order>
    {
        /** @return value_ at rvar_t<T,target_order>
         */
        static auto       &get(rvar<T, target_order> &v) { return v; }
        /** @return const value_ at rvar_t<T,target_order>
         */
        static const auto &get(const rvar<T, target_order> &v) { return v; }
    };
    /*****************************************************************************************/
    void make_leaf_node()
    {
        gradient_tape<T, order, BOOST_MATH_BUFFER_SIZE> &tape = get_active_tape<T, order>();
        node_                                      = tape.emplace_leaf_node();
    }

    void make_unary_node()
    {
        gradient_tape<T, order, BOOST_MATH_BUFFER_SIZE> &tape = get_active_tape<T, order>();
        node_                                      = tape.emplace_active_unary_node();
    }

    void make_multi_node(size_t n)
    {
        gradient_tape<T, order, BOOST_MATH_BUFFER_SIZE> &tape = get_active_tape<T, order>();
        node_                                      = tape.emplace_active_multi_node(n);
    }

    template<size_t n>
    void make_multi_node()
    {
        gradient_tape<T, order, BOOST_MATH_BUFFER_SIZE> &tape = get_active_tape<T, order>();
        node_                                      = tape.template emplace_active_multi_node<n>();
    }

    template<typename E>
    void make_rvar_from_expr(const expression<T, order, E> &expr)
    {
        make_multi_node<detail::count_rvars<E, order>>();
        expr.template propagatex<0>(node_, inner_t(1.0));
    }
    T get_item_impl(std::true_type) const
    {
        return value_.get_item_impl(std::integral_constant<bool, (order - 1 > 1)>{});
    }

    T get_item_impl(std::false_type) const { return value_; }

public:
    using value_type                = T;
    static constexpr size_t order_v = order;
    rvar()
        : value_()
    {
        make_leaf_node();
    }
    rvar(const T value)
        : value_(inner_t{value})
    {
        make_leaf_node();
    }

    rvar &operator=(T v)
    {
        value_ = inner_t(v);
        if (node_ == nullptr) {
            make_leaf_node();
        }
        return *this;
    }
    rvar(const rvar<T, order> &other)            = default;
    rvar &operator=(const rvar<T, order> &other) = default;

    template<size_t arg_index>
    void propagatex(gradient_node<T, order> *node, inner_t adj) const
    {
        node->update_derivative_v(arg_index, adj);
        node->update_argument_ptr_at(arg_index, node_);
    }

    template<class E>
    rvar(const expression<T, order, E> &expr)
    {
        value_ = expr.evaluate();
        make_rvar_from_expr(expr);
    }
    template<class E>
    rvar &operator=(const expression<T, order, E> &expr)
    {
        value_ = expr.evaluate();
        make_rvar_from_expr(expr);
        return *this;
    }
    /***************************************************************************************************/
    template<class E>
    rvar<T, order> &operator+=(const expression<T, order, E> &expr)
    {
        *this = *this + expr;
        return *this;
    }

    template<class E>
    rvar<T, order> &operator*=(const expression<T, order, E> &expr)
    {
        *this = *this * expr;
        return *this;
    }

    template<class E>
    rvar<T, order> &operator-=(const expression<T, order, E> &expr)
    {
        *this = *this - expr;
        return *this;
    }

    template<class E>
    rvar<T, order> &operator/=(const expression<T, order, E> &expr)
    {
        *this = *this / expr;
        return *this;
    }
    /***************************************************************************************************/
    rvar<T, order> &operator+=(const T &v)
    {
        *this = *this + v;
        return *this;
    }

    rvar<T, order> &operator*=(const T &v)
    {
        *this = *this * v;
        return *this;
    }

    rvar<T, order> &operator-=(const T &v)
    {
        *this = *this - v;
        return *this;
    }

    rvar<T, order> &operator/=(const T &v)
    {
        *this = *this / v;
        return *this;
    }

    /***************************************************************************************************/
    const inner_t &adjoint() const { return *node_->get_adjoint_ptr(); }
    inner_t       &adjoint() { return *node_->get_adjoint_ptr(); }

    const inner_t &evaluate() const { return value_; };
    inner_t       &get_value() { return value_; };

    explicit operator T() const { return item(); }

    explicit       operator int() const { return static_cast<int>(item()); }
    explicit       operator long() const { return static_cast<long>(item()); }
    explicit       operator long long() const { return static_cast<long long>(item()); }

    /**
     *  @brief same as evaluate but returns proper depth for higher order derivatives
     *  @return value_ at depth N
     */
    template<size_t N>
    auto &get_value_at()
    {
        static_assert(N <= order, "Requested depth exceeds variable order.");
        return get_value_at_impl<N, order>::get(*this);
    }
    /** @brief same as above but const
     */
    template<size_t N>
    const auto &get_value_at() const
    {
        static_assert(N <= order, "Requested depth exceeds variable order.");
        return get_value_at_impl<N, order>::get(*this);
    }

    T    item() const { return get_item_impl(std::integral_constant<bool, (order > 1)>{}); }

    void backward()
    {
        gradient_tape<T, order, BOOST_MATH_BUFFER_SIZE> &tape = get_active_tape<T, order>();
        auto                                  it   = tape.find(node_);
        it->update_adjoint_v(inner_t(1.0));
        while (it != tape.begin()) {
            it->backward();
            --it;
        }
        it->backward();
    }
};

template<typename T, size_t order>
std::ostream &operator<<(std::ostream &os, const rvar<T, order> var)
{
    os << "rvar<" << order << ">(" << var.item() << "," << var.adjoint() << ")";
    return os;
}

template<typename T, size_t order, typename E>
std::ostream &operator<<(std::ostream &os, const expression<T, order, E> &expr)
{
    rvar<T, order> tmp = expr;
    os << "rvar<" << order << ">(" << tmp.item() << "," << tmp.adjoint() << ")";
    return os;
}

template<typename T, size_t order>
rvar<T, order> make_rvar(const T v)
{
    static_assert(order > 0, "rvar order must be >= 1");
    return rvar<T, order>(v);
}
template<typename T, size_t order, typename E>
rvar<T, order> make_rvar(const expression<T, order, E> &expr)
{
    static_assert(order > 0, "rvar order must be >= 1");
    return rvar<T, order>(expr);
}

namespace detail {

/** @brief helper overload for grad implementation.
 *  @return vector<rvar<T,order-1> of gradients of the autodiff graph.
 *  specialization for autodiffing through autodiff. i.e. being able to
 *  compute higher order grads
*/
template<typename T, size_t order>
struct grad_op_impl
{
    std::vector<rvar<T, order - 1>> operator()(rvar<T, order> &f, std::vector<rvar<T, order> *> &x)
    {
        auto &tape = get_active_tape<T, order>();
        tape.zero_grad();
        f.backward();

        std::vector<rvar<T, order - 1>> gradient_vector;
        gradient_vector.reserve(x.size());

        for (auto xi : x) {
            // make a new rvar<T,order-1> holding the adjoint value
            gradient_vector.emplace_back(xi->adjoint());
        }
        return gradient_vector;
    }
    /*
    std::vector<rvar_t<T, order - 1> *> operator()(rvar<T, order>                &f,
                                                   std::vector<rvar<T, order> *> &x)
    {
        gradient_tape<T, order, BOOST_MATH_BUFFER_SIZE> &tape = get_active_tape<T, order>();
        tape.zero_grad();
        f.backward();
        std::vector<rvar_t<T, order - 1> *> gradient_vector;
        for (auto xi : x) {
            gradient_vector.push_back(&(xi->adjoint()));
        }
        return gradient_vector;
    }*/
};
/** @brief helper overload for grad implementation.
 *  @return vector<T> of gradients of the autodiff graph.
 *          base specialization for order 1 autodiff
*/
template<typename T>
struct grad_op_impl<T, 1>
{
    std::vector<T> operator()(rvar<T, 1> &f, std::vector<rvar<T, 1> *> &x)
    {
        gradient_tape<T, 1, BOOST_MATH_BUFFER_SIZE> &tape = get_active_tape<T, 1>();
        tape.zero_grad();
        f.backward();
        std::vector<T> gradient_vector;
        for (auto xi : x) {
            gradient_vector.push_back(xi->adjoint());
        }
        return gradient_vector;
    }
};

/** @brief helper overload for higher order autodiff
 *  @return nested vector representing N-d tensor of
 *      higher order derivatives
 */
template<size_t N, typename T, size_t order_1, size_t order_2, typename Enable = void>
struct grad_nd_impl
{
    auto operator()(rvar<T, order_1> &f, std::vector<rvar<T, order_2> *> &x)
    {
        static_assert(N > 1, "N must be greater than 1 for this template");

        auto current_grad = grad(f, x); // vector<rvar<T,order_1-1>> or vector<T>

        std::vector<decltype(grad_nd_impl<N - 1, T, order_1 - 1, order_2>()(current_grad[0], x))>
            result;
        result.reserve(current_grad.size());

        for (auto &g_i : current_grad) {
            result.push_back(grad_nd_impl<N - 1, T, order_1 - 1, order_2>()(g_i, x));
        }
        return result;
    }
    /*
    auto operator()(rvar<T, order_1> &f, std::vector<rvar<T, order_2> *> &x)
    {
        static_assert(N > 1, "N must be greater than 1 for this template");
        auto current_grad = grad(f, x);
        std::vector<decltype(grad_nd_impl<N - 1, T, order_1 - 1, order_2>()(*current_grad[0], x))>
            result;
        for (auto &g_i : current_grad) {
            result.push_back(grad_nd_impl<N - 1, T, order_1 - 1, order_2>()(*g_i, x));
        }
        return result;
    }*/
};
/** @brief spcialization for order = 1,
 *  @return vector<rvar<T,order_1-1>> gradients */
template<typename T, size_t order_1, size_t order_2>
struct grad_nd_impl<1, T, order_1, order_2>
{
    auto operator()(rvar<T, order_1> &f, std::vector<rvar<T, order_2> *> &x) { return grad(f, x); }
};

template<typename ptr>
struct rvar_order;

template<typename T, size_t order>
struct rvar_order<rvar<T, order> *>
{
    static constexpr size_t value = order;
};

} // namespace detail

/**
 * @brief grad computes gradient with respect to vector of pointers x
 * @param f -> computational graph
 * @param x -> variables gradients to record. Note ALL gradients of the graph
 *             are computed simultaneously, only the ones w.r.t. x are returned
 * @return vector<rvar<T,order_1 - 1> of gradients. in the case of order_1 = 1
 *            rvar<T,order_1-1> decays to T
 *
 * safe to call recursively with grad(grad(grad...
 */
template<typename T, size_t order_1, size_t order_2>
auto grad(rvar<T, order_1> &f, std::vector<rvar<T, order_2> *> &x)
{
    static_assert(order_1 <= order_2,
                  "variable differentiating w.r.t. must have order >= function order");
    std::vector<rvar<T, order_1> *> xx;
    for (auto xi : x)
        xx.push_back(&(xi->template get_value_at<order_1>()));
    return detail::grad_op_impl<T, order_1>{}(f, xx);
}
/** @brief variadic overload of above
 */
template<typename T, size_t order_1, typename First, typename... Other>
auto grad(rvar<T, order_1> &f, First first, Other... other)
{
    constexpr size_t order_2 = detail::rvar_order<First>::value;
    static_assert(order_1 <= order_2,
                  "variable differentiating w.r.t. must have order >= function order");
    std::vector<rvar<T, order_2> *> x_vec = {first, other...};
    return grad(f, x_vec);
}

/** @brief computes hessian matrix of computational graph w.r.t.
 *         vector of variables x.
 *  @return std::vector<std::vector<rvar<T,order_1-2>> hessian matrix
 *          rvar<T,2> decays to T
 *
 *  NOT recursion safe, cannot do hess(hess(
 */
template<typename T, size_t order_1, size_t order_2>
auto hess(rvar<T, order_1> &f, std::vector<rvar<T, order_2> *> &x)
{
    return detail::grad_nd_impl<2, T, order_1, order_2>{}(f, x);
}
/** @brief variadic overload of above
 */
template<typename T, size_t order_1, typename First, typename... Other>
auto hess(rvar<T, order_1> &f, First first, Other... other)
{
    constexpr size_t                order_2 = detail::rvar_order<First>::value;
    std::vector<rvar<T, order_2> *> x_vec   = {first, other...};
    return hess(f, x_vec);
}

/** @brief comput N'th gradient of computational graph w.r.t. x
 *  @return vector<vector<.... up N nestings representing tensor
 *          of gradients of order N
 *
 *  NOT recursively safe, cannot do grad_nd(grad_nd(... etc...
 */
template<size_t N, typename T, size_t order_1, size_t order_2>
auto grad_nd(rvar<T, order_1> &f, std::vector<rvar<T, order_2> *> &x)
{
    static_assert(order_1 >= N, "Function order must be at least N");
    static_assert(order_2 >= order_1, "Variable order must be at least function order");

    return detail::grad_nd_impl<N, T, order_1, order_2>()(f, x);
}

/** @brief variadic overload of above
 */
template<size_t N, typename ftype, typename First, typename... Other>
auto grad_nd(ftype &f, First first, Other... other)
{
    using T                                 = typename ftype::value_type;
    constexpr size_t                order_1 = detail::rvar_order<ftype *>::value;
    constexpr size_t                order_2 = detail::rvar_order<First>::value;
    std::vector<rvar<T, order_2> *> x_vec   = {first, other...};
    return detail::grad_nd_impl<N, T, order_1, order_1>{}(f, x_vec);
}
} // namespace reverse_mode
} // namespace differentiation
} // namespace math
} // namespace boost
namespace std {

// copied from forward mode
template<typename T, size_t order>
class numeric_limits<boost::math::differentiation::reverse_mode::rvar<T, order>>
    : public numeric_limits<
          typename boost::math::differentiation::reverse_mode::rvar<T, order>::value_type>
{};
} // namespace std
#endif
