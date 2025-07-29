#ifndef BOOST_MATH_DIFFERENTIATION_AUTODIFF_HPP
#define BOOST_MATH_DIFFERENTIATION_AUTODIFF_HPP

#include <array>
#include <boost/cstdfloat.hpp>
#include <boost/math/constants/constants.hpp>
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
#include <iterator>
#include <memory>
#include <type_traits>
#include <vector>
#define BUFFER_SIZE 16
//constexpr size_t MAX_DEPTH = 21;

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
struct count_rvar_impl<boost::math::differentiation::reverse_mode::rvar<U, order>, order>
{
    static constexpr std::size_t value = 1;
};

template<typename T, std::size_t order>
struct count_rvar_impl<T,
                       order,
                       std::enable_if_t<has_binary_sub_types<T>::value
                                        && !std::is_same<T,
                                                         boost::math::differentiation::reverse_mode::
                                                             rvar<typename T::value_type, order>>::value
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
                              && !std::is_same<T,
                                               boost::math::differentiation::reverse_mode::
                                                   rvar<typename T::value_type, order>>::value
                              && !has_binary_sub_types<T>::value>>
{
    static constexpr std::size_t value = count_rvar_impl<typename T::arg_type, order>::value;
};
template<typename T, size_t order>
constexpr std::size_t count_rvars = detail::count_rvar_impl<T, order>::value;

/* memory management helps for tape */
template<typename allocator_type, size_t buffer_size>
class flat_linear_allocator_iterator
{
    /**
   * @brief enables iterating over linear allocator with
   * c++ iterators
   */
public:
    using raw_allocator_type   = std::remove_const_t<allocator_type>;
    using value_type           = typename allocator_type::value_type;
    using pointer              = typename allocator_type::value_type*;
    using const_ptr_type       = const value_type*;
    using reference            = typename allocator_type::value_type&;
    using const_reference_type = const value_type&;
    using iterator_category    = std::random_access_iterator_tag;
    using difference_type      = ptrdiff_t;

private:
    size_t                begin_   = 0;
    size_t                index_   = 0;
    size_t                end_     = 0;
    const allocator_type* storage_ = nullptr;

public:
    explicit flat_linear_allocator_iterator(allocator_type* storage, size_t index)
        : storage_(storage)
        , index_(index)
        , begin_(0)
        , end_(storage->size())
    {}

    explicit flat_linear_allocator_iterator(allocator_type* storage,
                                            size_t          index,
                                            size_t          begin,
                                            size_t          end)
        : storage_(storage)
        , index_(index)
        , begin_(begin)
        , end_(end)
    {}

    explicit flat_linear_allocator_iterator(const allocator_type* storage, size_t index)
        : storage_(storage)
        , index_(index)
        , begin_(0)
        , end_(storage->size())
    {}

    explicit flat_linear_allocator_iterator(const allocator_type* storage,
                                            size_t                index,
                                            size_t                begin,
                                            size_t                end)
        : storage_(storage)
        , index_(index)
        , begin_(begin)
        , end_(end)
    {}
    reference operator*()
    {
        assert(index_ >= begin_ && index_ < end_);
        return (*storage_->data_[index_ / buffer_size])[index_ % buffer_size];
    }

    const_reference_type operator*() const
    {
        assert(index_ >= begin_ && index_ < end_);
        return (*storage_->data_[index_ / buffer_size])[index_ % buffer_size];
    }

    pointer operator->()
    {
        assert(index_ >= begin_ && index_ < end_);
        return &operator*();
    }

    const_ptr_type operator->() const
    {
        assert(index_ >= begin_ && index_ < end_);
        return &operator*();
    }
    flat_linear_allocator_iterator& operator++()
    {
        ++index_;
        return *this;
    }

    flat_linear_allocator_iterator operator++(int)
    {
        auto tmp = *this;
        ++(*this);
        return tmp;
    }

    flat_linear_allocator_iterator& operator--()
    {
        --index_;
        return *this;
    }

    flat_linear_allocator_iterator operator--(int)
    {
        auto tmp = *this;
        --(*this);
        return tmp;
    }

    bool operator==(const flat_linear_allocator_iterator& other) const
    {
        return index_ == other.index_ && storage_ == other.storage_;
    }

    bool operator!=(const flat_linear_allocator_iterator& other) const { return !(*this == other); }

    flat_linear_allocator_iterator operator+(difference_type n) const
    {
        return flat_linear_allocator_iterator(storage_, index_ + n, begin_, end_);
    }

    flat_linear_allocator_iterator& operator+=(difference_type n)
    {
        index_ += n;
        return *this;
    }

    flat_linear_allocator_iterator operator-(difference_type n) const
    {
        return flat_linear_allocator_iterator(storage_, index_ - n, begin_, end_);
    }
    flat_linear_allocator_iterator& operator-=(difference_type n)
    {
        index_ -= n;
        return *this;
    }

    difference_type operator-(const flat_linear_allocator_iterator& other) const
    {
        return static_cast<difference_type>(index_) - static_cast<difference_type>(other.index_);
    }

    reference operator[](difference_type n) { return *(*this + n); }

    const_reference_type operator[](difference_type n) const { return *(*this + n); }

    bool operator<(const flat_linear_allocator_iterator& other) const
    {
        return index_ < other.index_;
    }

    bool operator>(const flat_linear_allocator_iterator& other) const
    {
        return index_ > other.index_;
    }

    bool operator<=(const flat_linear_allocator_iterator& other) const
    {
        return index_ <= other.index_;
    }

    bool operator>=(const flat_linear_allocator_iterator& other) const
    {
        return index_ >= other.index_;
    }
};
/* memory management helps for tape */
template<typename T, size_t buffer_size>
class flat_linear_allocator
{
    /*
     * @brief basically a vector<array<T*, size>> 
     * intended to work like a vector that allocates memory in chunks
     * and doesn't invalidate references
     * */
public:
    // store vector of unique pointers to arrays
    // to avoid vector reference invalidation
    using buffer_type = std::array<T, buffer_size>;
    using buffer_ptr  = std::unique_ptr<buffer_type>;

private:
    std::vector<buffer_ptr> data_;
    size_t                  total_size_ = 0;
    std::vector<size_t>     checkpoints_; //{0};

public:
    friend class flat_linear_allocator_iterator<flat_linear_allocator<T, buffer_size>, buffer_size>;
    friend class flat_linear_allocator_iterator<const flat_linear_allocator<T, buffer_size>,
                                                buffer_size>;
    using value_type = T;
    using iterator
        = flat_linear_allocator_iterator<flat_linear_allocator<T, buffer_size>, buffer_size>;
    using const_iterator
        = flat_linear_allocator_iterator<const flat_linear_allocator<T, buffer_size>, buffer_size>;

    size_t buffer_id() const { return total_size_ / buffer_size; }
    size_t item_id() const { return total_size_ % buffer_size; }

private:
    void allocate_buffer() { data_.emplace_back(std::make_unique<buffer_type>()); }

public:
    flat_linear_allocator() { allocate_buffer(); }
    flat_linear_allocator(const flat_linear_allocator&)            = delete;
    flat_linear_allocator& operator=(const flat_linear_allocator&) = delete;
    flat_linear_allocator(flat_linear_allocator&&)                 = default;
    flat_linear_allocator& operator=(flat_linear_allocator&&)      = default;
    ~flat_linear_allocator()                                       = default;

    void clear()
    {
        data_.clear();
        total_size_ = 0;
        checkpoints_.clear();
        allocate_buffer();
    }

    // doesn't delete anything, only sets the current index to zero
    void reset() { total_size_ = 0; }
    void rewind() { total_size_ = 0; };
    // adds current index as a checkpoint to be able to walk back to
    void add_checkpoint()
    {
        if (total_size_ > 0) {
            checkpoints_.push_back(total_size_ - 1);
        }
    };

    /* @brief clears all checkpoints
     * */
    void reset_checkpoints() { checkpoints_.clear(); }

    void rewind_to_last_checkpoint() { total_size_ = checkpoints_.back(); }
    void rewind_to_checkpoint_at(size_t index) { total_size_ = checkpoints_[index]; }

    void fill(const T& val)
    {
        for (size_t i = 0; i < total_size_; ++i) {
            size_t bid         = i / buffer_size;
            size_t iid         = i % buffer_size;
            (*data_[bid])[iid] = val;
        }
    }

    /* @brief emplaces back object at the end of the
     * data structure, calls default constructor */
    iterator emplace_back()
    {
        if (item_id() == 0 && total_size_ != 0) {
            allocate_buffer();
        }
        size_t bid = buffer_id();
        size_t iid = item_id();

        T* ptr = &(*data_[bid])[iid];
        new (ptr) T();
        ++total_size_;
        return iterator(this, total_size_ - 1);
    };

    /* @brief, emplaces back object at end of data structure,
     * passes arguments to constructor */
    template<typename... Args>
    iterator emplace_back(Args&&... args)
    {
        if (item_id() == 0 && total_size_ != 0) {
            allocate_buffer();
        }
        assert(buffer_id() < data_.size());
        assert(item_id() < buffer_size);
        T* ptr = &(*data_[buffer_id()])[item_id()];
        new (ptr) T(std::forward<Args>(args)...);
        ++total_size_;
        return iterator(this, total_size_ - 1);
    };
    /* @brief default constructs n objects at end of
     * data structure, n known at compile time */
    template<size_t n>
    iterator emplace_back_n()
    {
        size_t bid = buffer_id();
        size_t iid = item_id();
        if (iid + n < buffer_size) {
            T* ptr = &(*data_[bid])[iid];
            for (size_t i = 0; i < n; ++i) {
                new (ptr + i) T();
            }
            total_size_ += n;
            return iterator(this, total_size_ - n, total_size_ - n, total_size_);
        } else {
            size_t allocs_in_curr_buffer = buffer_size - iid;
            size_t allocs_in_next_buffer = n - (buffer_size - iid);
            T*     ptr                   = &(*data_[bid])[iid];
            for (size_t i = 0; i < allocs_in_curr_buffer; ++i) {
                new (ptr + i) T();
            }
            allocate_buffer();
            bid = data_.size() - 1;
            iid = 0;
            total_size_ += n;

            T* ptr2 = &(*data_[bid])[iid];
            for (size_t i = 0; i < allocs_in_next_buffer; i++) {
                new (ptr2 + i) T();
            }
            return iterator(this, total_size_ - n, total_size_ - n, total_size_);
        }
    }
    /* @brief default constructs n objects at end of
     * data structure, n known at run time
     */
    iterator emplace_back_n(size_t n)
    {
        size_t bid = buffer_id();
        size_t iid = item_id();
        if (iid + n < buffer_size) {
            T* ptr = &(*data_[bid])[iid];
            for (size_t i = 0; i < n; ++i) {
                new (ptr + i) T();
            }
            total_size_ += n;
            return iterator(this, total_size_ - n, total_size_ - n, total_size_);
        } else {
            size_t allocs_in_curr_buffer = buffer_size - iid;
            size_t allocs_in_next_buffer = n - (buffer_size - iid);
            T*     ptr                   = &(*data_[bid])[iid];
            for (size_t i = 0; i < allocs_in_curr_buffer; ++i) {
                new (ptr + i) T();
            }
            allocate_buffer();
            bid = data_.size() - 1;
            iid = 0;
            total_size_ += n;

            T* ptr2 = &(*data_[bid])[iid];
            for (size_t i = 0; i < allocs_in_next_buffer; i++) {
                new (ptr2 + i) T();
            }
            return iterator(this, total_size_ - n, total_size_ - n, total_size_);
        }
    }

    size_t size() const { return total_size_; }
    size_t capacity() const { return data_.size() * buffer_size; }

    iterator       begin() { return iterator(this, 0); }
    iterator       end() { return iterator(this, total_size_); }
    const_iterator begin() const { return const_iterator(this, 0); }
    const_iterator end() const { return const_iterator(this, total_size_); }

    iterator last_checkpoint() { return iterator(this, checkpoints_.back(), 0, total_size_); }
    iterator first_checkpoint() { return iterator(this, checkpoints_[0], 0, total_size_); };
    iterator checkpoint_at(size_t index)
    {
        return iterator(this, checkpoints_[index], 0, total_size_);
    };

    iterator find(const T* const item)
    {
        return std::find_if(begin(), end(), [&](const T& val) { return &val == item; });
    }
    T& operator[](std::size_t i)
    {
        assert(i < total_size_);
        return (*data_[i / buffer_size])[i % buffer_size];
    }
    const T& operator[](std::size_t i) const
    {
        assert(i < total_size_);
        return (*data_[i / buffer_size])[i % buffer_size];
    }
};

}; // namespace detail
namespace detail {
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

/* rvar_t<T,0> decays to T
 * otherwise its just identity rvar_t<T,N> = rvar<T,N>
 * */
template<typename T, size_t N>
using rvar_t = typename detail::rvar_type_impl<T, N>::type;

// manages nodes in computational graph
template<typename T, size_t order, size_t buffer_size = BUFFER_SIZE>
class gradient_tape
{
private:
    using inner_t = rvar_t<T, order - 1>;
    detail::flat_linear_allocator<inner_t, buffer_size>                  adjoints_;
    detail::flat_linear_allocator<inner_t, buffer_size>                  derivatives_;
    detail::flat_linear_allocator<gradient_node<T, order>, buffer_size>  gradient_nodes_;
    detail::flat_linear_allocator<gradient_node<T, order>*, buffer_size> argument_nodes_;

    // compile time check if emplace_back calls on zero
    template<size_t n>
    gradient_node<T, order>* fill_node_at_compile_time(std::true_type,
                                                       gradient_node<T, order>* node_ptr)
    {
        node_ptr->derivatives_       = &*derivatives_.template emplace_back_n<n>();
        node_ptr->argument_nodes_    = &*argument_nodes_.template emplace_back_n<n>();
        return node_ptr;
    }

    template<size_t n>
    gradient_node<T, order>* fill_node_at_compile_time(std::false_type,
                                                       gradient_node<T, order>* node_ptr)
    {
        node_ptr->derivatives_       = nullptr;
        node_ptr->argument_adjoints_ = nullptr;
        node_ptr->argument_nodes_    = nullptr;
        return node_ptr;
    }

public:
    gradient_tape() { clear(); };
    //size_t depth() const override { return depth_; };

    gradient_tape(const gradient_tape&)            = delete;
    gradient_tape& operator=(const gradient_tape&) = delete;
    gradient_tape(gradient_tape&& other)
        : adjoints_(std::move(other.adjoints_))
        , derivatives_(std::move(other.derivatives_))
        , gradient_nodes_(std::move(other.gradient_nodes_))
        , argument_nodes_(std::move(other.argument_nodes_))

    {
        other.clear();
    }
    gradient_tape operator=(gradient_tape&& other)
    {
        if (this != &other) {
            adjoints_          = std::move(other.adjoints_);
            derivatives_       = std::move(other.derivatives_);
            gradient_nodes_    = std::move(other.gradient_nodes_);
            argument_nodes_    = std::move(other.argument_nodes_);
        }
        return *this;
    }
    ~gradient_tape() noexcept {
        clear();
    }
    void clear() noexcept
    {
        adjoints_.clear();
        derivatives_.clear();
        gradient_nodes_.clear();
        argument_nodes_.clear();
    }

    // no derivatives or arguments
    gradient_node<T, order>* emplace_leaf_node()
    {
        gradient_node<T, order>* node = &*gradient_nodes_.emplace_back();
        node->adjoint_           = &*adjoints_.emplace_back();
        node->derivatives_       = nullptr;
        node->argument_adjoints_ = nullptr;
        node->argument_nodes_    = nullptr;

        return node;
    };

    // single argument, single derivative
    gradient_node<T, order>* emplace_active_unary_node()
    {
        gradient_node<T, order>* node = &*gradient_nodes_.emplace_back();
        node->n_                 = 1;
        node->adjoint_           = &*adjoints_.emplace_back();
        node->derivatives_       = &*derivatives_.emplace_back();

        return node;
    };

    // arbitrary number of arguments/derivatives (compile time)
    template<size_t n>
    gradient_node<T, order>* emplace_active_multi_node()
    {
        gradient_node<T, order>* node = &*gradient_nodes_.emplace_back();
        node->n_               = n;
        node->adjoint_         = &*adjoints_.emplace_back();
        // emulate if constexpr
        return fill_node_at_compile_time<n>(std::integral_constant<bool, (n > 0)>{}, node);
    };

    // same as above at runtime
    gradient_node<T, order>* emplace_active_multi_node(size_t n)
    {
        gradient_node<T, order>* node = &*gradient_nodes_.emplace_back();
        node->n_               = n;
        node->adjoint_         = &*adjoints_.emplace_back();
        if (n > 0) {
            node->derivatives_       = &*derivatives_.emplace_back_n(n);
            node->argument_nodes_    = &*argument_nodes_.emplace_back_n(n);
        }
        return node;
    };

    void zero_grad()
    {
        const T zero = T(0.0);
        adjoints_.fill(zero);
    }

    // return type is an iterator
    auto begin() { return gradient_nodes_.begin(); }
    auto end() { return gradient_nodes_.end(); }
    auto find(gradient_node<T, order>* node) { return gradient_nodes_.find(node); };
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
    void rewind_to_checkpoint_at(
        size_t index) // index is "checkpoint" index. so order which checkpoint was set
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
    gradient_node<T, order>&       operator[](size_t i) { return gradient_nodes_[i]; }
    const gradient_node<T, order>& operator[](size_t i) const { return gradient_nodes_[i]; }
};
//class rvar;
template<typename T, size_t order> // no CRTP, just storage
class gradient_node
{
    /*
     * @brief manages adjoints, derivatives, and stores points to argument adjoints
     * pointers to arguments aren't needed here
     * */
private:
    size_t n_;
    using inner_t = rvar_t<T, order - 1>;
    inner_t*                  adjoint_;
    inner_t*                  derivatives_;
    inner_t**                 argument_adjoints_;
    gradient_node<T, order>** argument_nodes_;

public:
    friend class gradient_tape<T, order>;
    friend class rvar<T, order>;

    gradient_node() = default;
    explicit gradient_node(const size_t n)
        : n_(n)
        , adjoint_(nullptr)
        , derivatives_(nullptr)
        , argument_adjoints_(nullptr)
    {}
    explicit gradient_node(const size_t     n,
                           T*               adjoint,
                           T*               derivatives,
                           T**              argument_adjoints,
                           rvar<T, order>** arguments)
        : n_(n)
        , adjoint_(adjoint)
        , derivatives_(derivatives)
        , argument_adjoints_(argument_adjoints)
    {}

    inner_t get_adjoint_v() const { return *adjoint_; }
    inner_t get_derivative_v(size_t arg_id) const { return derivatives_[arg_id]; };
    inner_t get_argument_adjoint_v(size_t arg_id) const
    {
        return *argument_nodes_[arg_id]->adjoint_;
    }

    inner_t* get_adjoint_ptr() const { return adjoint_; };
    inner_t* get_derivative_ptr() const { return derivatives_; };
    inner_t* get_argument_adjoint_ptr(size_t arg_id) const
    {
        return argument_nodes_[arg_id]->adjoint_;
    };

    void update_adjoint_v(inner_t value) { *adjoint_ = value; };
    void update_derivative_v(size_t arg_id, inner_t value) { derivatives_[arg_id] = value; };
    void update_argument_adj_v(size_t arg_id, inner_t value)
    {
        argument_nodes_[arg_id]->update_adjoint_v(value);
    };
    void update_argument_ptr_at(size_t arg_id, gradient_node<T, order>* node_ptr)
    {
        argument_nodes_[arg_id] = node_ptr;
    }

    void backward()
    {
        if (!n_)
            return;

        if (!adjoint_ || fabs(*adjoint_) < std::numeric_limits<T>::epsilon())
            return;

        if (!argument_nodes_)
            return;

        if (!derivatives_)
            return;

        for (size_t i = 0; i < n_; ++i) {
            auto adjoint          = get_adjoint_v();
            auto derivative       = get_derivative_v(i);
            auto argument_adjoint = get_argument_adjoint_v(i);
            update_argument_adj_v(i, argument_adjoint + derivative * adjoint);
        }
    }
};
template<typename T, size_t order, class derived_expression>
struct expression
{
    /* @brief
     * base expression class
     * */
    static constexpr size_t num_literals = 0;
    using inner_t                        = rvar_t<T, order - 1>;
    inner_t evaluate() const { return static_cast<const derived_expression*>(this)->evaluate(); }

    template<size_t arg_index>
    void propagatex(gradient_node<T, order>* node, inner_t adj) const
    {
        return static_cast<const derived_expression*>(this)->template propagatex<arg_index>(node,
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

    explicit abstract_binary_expression(const expression<T, order, LHS>& left_hand_expr,
                                        const expression<T, order, RHS>& right_hand_expr)
        : lhs(static_cast<const LHS&>(left_hand_expr))
        , rhs(static_cast<const RHS&>(right_hand_expr)){};

    inner_t evaluate() const
    {
        return static_cast<const concrete_binary_operation*>(this)->evaluate();
    };

    template<size_t arg_index>
    void propagatex(gradient_node<T, order>* node, inner_t adj) const
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
    void propagate_lhs(gradient_node<T, order>* node, inner_t adj) const
    {
        lhs.template propagatex<arg_index_>(node, adj);
    }

    template<std::size_t num_args,
             std::size_t arg_index_,
             typename std::enable_if<(num_args == 0), int>::type = 0>
    void propagate_lhs(gradient_node<T, order>*, inner_t) const
    {}

    template<std::size_t num_args,
             std::size_t arg_index_,
             typename std::enable_if<(num_args > 0), int>::type = 0>
    void propagate_rhs(gradient_node<T, order>* node, inner_t adj) const
    {
        rhs.template propagatex<arg_index_>(node, adj);
    }

    template<std::size_t num_args,
             std::size_t arg_index_,
             typename std::enable_if<(num_args == 0), int>::type = 0>
    void propagate_rhs(gradient_node<T, order>*, inner_t) const
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
    const T constant;
    explicit abstract_unary_expression(const expression<T, order, ARG>& arg_expr, const T& constant)

        : arg(static_cast<const ARG&>(arg_expr)), constant(constant) {};

    inner_t evaluate() const
    {
        return static_cast<const concrete_unary_operation*>(this)->evaluate();
    };

    template<size_t arg_index>
    void propagatex(gradient_node<T, order>* node, inner_t adj) const
    {
        inner_t argv        = arg.evaluate();
        inner_t v         = evaluate();
        inner_t partial_arg = concrete_unary_operation::derivative(argv, v, constant);

        constexpr size_t num_arg_args = detail::count_rvars<ARG, order>;

        arg.template propagatex<arg_index>(node, adj);
    }

};
/****************************************************************************************************************/
template<typename T, size_t order, typename LHS, typename RHS>
struct add_expr
    : public abstract_binary_expression<T, order, LHS, RHS, add_expr<T, order, LHS, RHS>>
{
    /* @brief addition
     * rvar+rvar
     * */
    using lhs_type   = LHS;
    using rhs_type   = RHS;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;
    // Explicitly define constructor to forward to base class
    explicit add_expr(const expression<T, order, LHS>& left_hand_expr,
                      const expression<T, order, RHS>& right_hand_expr)
        : abstract_binary_expression<T, order, LHS, RHS, add_expr<T, order, LHS, RHS>>(
              left_hand_expr, right_hand_expr)
    {}

    inner_t              evaluate() const { return this->lhs.evaluate() + this->rhs.evaluate(); }
    static const inner_t left_derivative(const inner_t& l, const inner_t& r, const inner_t& v)
    {
        return inner_t(1.0);
    }
    static const inner_t right_derivative(const inner_t& l, const inner_t& r, const inner_t& v)
    {
        return inner_t(1.0);
    }
};
template<typename T, size_t order, typename ARG>
struct add_const_expr : public abstract_unary_expression<T, order, ARG, add_const_expr<T, order, ARG>>
{
    /* @brief
     * rvar+float or float+rvar
     * */
    using arg_type = ARG;
    using value_type = T;
    using inner_t = rvar_t<T,order-1>;
    explicit add_const_expr(const expression<T, order, ARG>& arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, add_const_expr<T, order, ARG>>(arg_expr, v)
    {};

    inner_t evaluate() const {return this->arg.evaluate() + inner_t(this->constant);}
    static const inner_t derivative(const inner_t& argv, const inner_t& v, const T& constant)
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
    using lhs_type   = LHS;
    using rhs_type   = RHS;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;
    // Explicitly define constructor to forward to base class
    explicit mult_expr(const expression<T, order, LHS>& left_hand_expr,
                       const expression<T, order, RHS>& right_hand_expr)
        : abstract_binary_expression<T, order, LHS, RHS, mult_expr<T, order, LHS, RHS>>(
              left_hand_expr, right_hand_expr)
    {}

    inner_t              evaluate() const { return this->lhs.evaluate() * this->rhs.evaluate(); };
    static const inner_t left_derivative(const inner_t& l, const inner_t& r, const inner_t& v)
    {
        return r;
    };
    static const inner_t right_derivative(const inner_t& l, const inner_t& r, const inner_t& v)
    {
        return l;
    };
};
template<typename T, size_t order, typename ARG>
struct mult_const_expr : public abstract_unary_expression<T, order, ARG, mult_const_expr<T, order, ARG>>
{
    /* @brief
     * rvar+float or float+rvar
     * */
    using arg_type = ARG;
    using value_type = T;
    using inner_t = rvar_t<T,order-1>;

    explicit mult_const_expr(const expression<T, order, ARG>& arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, mult_const_expr<T, order, ARG>>(arg_expr, v)
    {};

    inner_t evaluate() const {return this->arg.evaluate() * inner_t(this->constant);}
    static const inner_t derivative(const inner_t& argv, const inner_t& v, const T& constant)
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
    using lhs_type   = LHS;
    using rhs_type   = RHS;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;
    // Explicitly define constructor to forward to base class
    explicit sub_expr(const expression<T, order, LHS>& left_hand_expr,
                      const expression<T, order, RHS>& right_hand_expr)
        : abstract_binary_expression<T, order, LHS, RHS, sub_expr<T, order, LHS, RHS>>(
              left_hand_expr, right_hand_expr)
    {}

    inner_t              evaluate() const { return this->lhs.evaluate() - this->rhs.evaluate(); }
    static const inner_t left_derivative(const inner_t& l, const inner_t& r, const inner_t& v)
    {
        return inner_t(1.0);
    }
    static const inner_t right_derivative(const inner_t& l, const inner_t& r, const inner_t& v)
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
    using lhs_type   = LHS;
    using rhs_type   = RHS;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;
    // Explicitly define constructor to forward to base class
    explicit div_expr(const expression<T, order, LHS>& left_hand_expr,
                       const expression<T, order, RHS>& right_hand_expr)
        : abstract_binary_expression<T, order, LHS, RHS, div_expr<T, order, LHS, RHS>>(
              left_hand_expr, right_hand_expr)
    {}

    inner_t              evaluate() const { return this->lhs.evaluate() / this->rhs.evaluate(); };
    static const inner_t left_derivative(const inner_t& l, const inner_t& r, const inner_t& v)
    {
        return 1.0 / r;
    };
    static const inner_t right_derivative(const inner_t& l, const inner_t& r, const inner_t& v)
    {
        return - l / (r*r);
    };
};
template<typename T, size_t order, typename ARG>
struct inv_expr : public abstract_unary_expression<T, order, ARG, inv_expr<T, order, ARG>>
{
    /* @brief
     * 1.0 / rvar
     * */
    using arg_type = ARG;
    using value_type = T;
    using inner_t = rvar_t<T,order-1>;

    explicit inv_expr(const expression<T, order, ARG>& arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, inv_expr<T, order, ARG>>(arg_expr, v)
    {};

    inner_t evaluate() const {return 1.0/this->arg.evaluate();}
    static const inner_t derivative(const inner_t& argv, const inner_t& v, const T& constant)
    {
        return inner_t(constant);
    }

};
/****************************************************************************************************************/
template<typename T, size_t order, typename LHS, typename RHS>
mult_expr<T, order, LHS, RHS> operator*(const expression<T, order, LHS>& lhs,
                                        const expression<T, order, RHS>& rhs)
{
    return mult_expr<T, order, LHS, RHS>(lhs, rhs);
}
template<typename T, size_t order, typename ARG>
mult_const_expr<T,order,ARG> operator*(const expression<T,order,ARG>& arg, const T& v)
{
    return mult_const_expr<T,order,ARG>(arg, v);
}
template<typename T, size_t order, typename ARG>
mult_const_expr<T,order,ARG> operator*(const T& v, const expression<T,order,ARG>& arg)
{
    return mult_const_expr<T,order,ARG>(arg, v);
}
/****************************************************************************************************************/
/* + */
template<typename T, size_t order, typename LHS, typename RHS>
add_expr<T, order, LHS, RHS> operator+(const expression<T, order, LHS>& lhs,
                                       const expression<T, order, RHS>& rhs)
{
    return add_expr<T, order, LHS, RHS>(lhs, rhs);
}
template<typename T, size_t order, typename ARG>
add_const_expr<T,order,ARG> operator+(const expression<T,order,ARG>& arg, const T& v)
{
    return add_const_expr<T,order,ARG>(arg, v);
}
template<typename T, size_t order, typename ARG>
add_const_expr<T,order,ARG> operator+(const T& v, const expression<T,order,ARG>& arg)
{
    return add_const_expr<T,order,ARG>(arg, v);
}
/****************************************************************************************************************/
/* - */

/* @brief negation (-1.0*rvar) */
template<typename T, size_t order, typename ARG>
mult_const_expr<T,order,ARG> operator-(const expression<T,order,ARG>& arg)
{
    return mult_const_expr<T,order,ARG>(arg, -1.0);
}

/* @brief subtraction rvar-rvar */
template<typename T, size_t order, typename LHS, typename RHS>
sub_expr<T,order,LHS,RHS> operator-(const expression<T, order, LHS>& lhs,
                                       const expression<T, order, RHS>& rhs)
{
    return sub_expr<T,order,LHS,RHS>(lhs, rhs);
}

/* @brief subtraction float - rvar */
template<typename T, size_t order, typename ARG>
add_const_expr<T,order,ARG> operator-(const expression<T,order,ARG>& arg, const T& v)
{
    /* rvar - float = rvar + (-float) */
    return add_const_expr<T,order,ARG>(arg, -1.0*v);
}

/* @brief subtraction rvar - float */
template<typename T, size_t order, typename ARG>
add_const_expr<T,order,ARG> operator-(const T& v, const expression<T,order,ARG>& arg)
{
    /* float - rvar = -rvar + float */
    return add_const_expr<T,order,ARG>(-arg, v);
}
/****************************************************************************************************************/
/* / */

/* @brief negation (-1.0*rvar) */
template<typename T, size_t order, typename ARG>
mult_const_expr<T,order,ARG> operator-(const expression<T,order,ARG>& arg)
{
    return mult_const_expr<T,order,ARG>(arg, -1.0);
}

/* @brief subtraction rvar-rvar */
template<typename T, size_t order, typename LHS, typename RHS>
sub_expr<T,order,LHS,RHS> operator-(const expression<T, order, LHS>& lhs,
                                       const expression<T, order, RHS>& rhs)
{
    return sub_expr<T,order,LHS,RHS>(lhs, rhs);
}

/* @brief subtraction float - rvar */
template<typename T, size_t order, typename ARG>
add_const_expr<T,order,ARG> operator-(const expression<T,order,ARG>& arg, const T& v)
{
    /* rvar - float = rvar + (-float) */
    return add_const_expr<T,order,ARG>(arg, -1.0*v);
}

/* @brief subtraction rvar - float */
template<typename T, size_t order, typename ARG>
add_const_expr<T,order,ARG> operator-(const T& v, const expression<T,order,ARG>& arg)
{
    /* float - rvar = -rvar + float */
    return add_const_expr<T,order,ARG>(-arg, v);
}
/****************************************************************************************************************/
template<typename T, size_t order>
inline gradient_tape<T, order, BUFFER_SIZE>& get_active_tape()
{
    //static thread_local gradient_tape<T, order, BUFFER_SIZE> tape;
    /* something strange happens with thread local storage destruction. this fix seems to work but i don't
     * quite understand it yet */
    static thread_local gradient_tape<T,order, BUFFER_SIZE>* tape = new gradient_tape<T,order,BUFFER_SIZE>();
    return *tape;
}

template<typename T, size_t order = 1>
class rvar : public expression<T, order, rvar<T, order>>
{
private:
    using inner_t = rvar_t<T, order - 1>;
    friend class gradient_node<T, order>;
    inner_t                  value_;
    gradient_node<T, order>* node_;
    template<typename, size_t>
    friend class rvar;

    void make_leaf_node()
    {
        gradient_tape<T, order, BUFFER_SIZE>& tape = get_active_tape<T, order>();
        node_                               = tape.emplace_leaf_node();
    }

    void make_unary_node()
    {
        gradient_tape<T, order, BUFFER_SIZE>& tape = get_active_tape<T, order>();
        node_                                      = tape.emplace_active_unary_node();
    }

    void make_multi_node(size_t n)
    {
        gradient_tape<T, order, BUFFER_SIZE>& tape = get_active_tape<T, order>();
        node_                                      = tape.emplace_active_multi_node(n);
    }

    template<size_t n>
    void make_multi_node()
    {
        gradient_tape<T, order, BUFFER_SIZE>& tape = get_active_tape<T, order>();
        node_                                      = tape.template emplace_active_multi_node<n>();
    }

    template<typename E>
    void make_rvar_from_expr(const expression<T, order, E>& expr)
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
    using value_type = T;

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
    rvar& operator=(T v)
    {
        value_ = inner_t(v);
        make_leaf_node();
    }

    template<size_t arg_index>
    void propagatex(gradient_node<T, order>* node, inner_t adj) const
    {
        node->update_derivative_v(arg_index, adj);
        node->update_argument_ptr_at(arg_index, node_);
    }

    template<class E>
    rvar(const expression<T, order, E>& expr)
    {
        value_ = expr.evaluate();
        make_rvar_from_expr(expr);
    }
    template<class E>
    rvar& operator=(const expression<T, order, E>& expr)
    {
        value_ = expr.evaluate();
        make_rvar_from_expr(expr);
    }

    const inner_t& adjoint() const { return *node_->get_adjoint_ptr(); }
    const inner_t& evaluate() const { return value_; };
    T              item() { return get_item_impl(std::integral_constant<bool, (order > 1)>{}); }

    void print_der_info()
    {
        for (size_t i = 0; i < node_->n_; i++) {
            std::cout << node_->derivatives_[i] << std::endl;
        }
    }
    void backward()
    {
        gradient_tape<T, order, BUFFER_SIZE>& tape = get_active_tape<T, order>();
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
rvar<T, order> make_rvar(const T v)
{
    static_assert(order > 0, "rvar order must be >= 1");
    return rvar<T, order>(v);
}
template<typename T, size_t order, typename E>
rvar<T, order> make_rvar(const expression<T, order, E>& expr)
{
    static_assert(order > 0, "rvar order must be >= 1");
    return rvar<T, order>(expr);
}
} // namespace reverse_mode
} // namespace differentiation
} // namespace math
} // namespace boost

#endif
