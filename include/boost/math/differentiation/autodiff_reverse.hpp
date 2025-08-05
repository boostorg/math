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
#define BUFFER_SIZE 60000
// constexpr size_t MAX_DEPTH = 21;

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
struct is_expression : std::false_type
{};

template<typename T, size_t order, typename Derived>
struct is_expression<expression<T, order, Derived>> : std::true_type
{};
template<typename T, size_t order>
struct is_expression<rvar<T, order>> : std::true_type
{};

template<typename T, typename... Rest>
struct all_same;

template<typename T>
struct all_same<T> : std::true_type
{};

template<typename T, typename U, typename... Rest>
struct all_same<T, U, Rest...>
    : std::integral_constant<bool, std::is_same<T, U>::value && all_same<T, Rest...>::value>
{};

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
    using pointer              = typename allocator_type::value_type *;
    using const_ptr_type       = const value_type *;
    using reference            = typename allocator_type::value_type &;
    using const_reference_type = const value_type &;
    using iterator_category    = std::random_access_iterator_tag;
    using difference_type      = ptrdiff_t;

private:
    size_t                begin_   = 0;
    size_t                index_   = 0;
    size_t                end_     = 0;
    const allocator_type *storage_ = nullptr;

public:
    flat_linear_allocator_iterator() = default;

    explicit flat_linear_allocator_iterator(allocator_type *storage, size_t index)
        : storage_(storage)
        , index_(index)
        , begin_(0)
        , end_(storage->size())
    {}

    explicit flat_linear_allocator_iterator(allocator_type *storage,
                                            size_t          index,
                                            size_t          begin,
                                            size_t          end)
        : storage_(storage)
        , index_(index)
        , begin_(begin)
        , end_(end)
    {}

    explicit flat_linear_allocator_iterator(const allocator_type *storage, size_t index)
        : storage_(storage)
        , index_(index)
        , begin_(0)
        , end_(storage->size())
    {}

    explicit flat_linear_allocator_iterator(const allocator_type *storage,
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
    flat_linear_allocator_iterator &operator++()
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

    flat_linear_allocator_iterator &operator--()
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

    bool operator==(const flat_linear_allocator_iterator &other) const
    {
        return index_ == other.index_ && storage_ == other.storage_;
    }

    bool operator!=(const flat_linear_allocator_iterator &other) const { return !(*this == other); }

    flat_linear_allocator_iterator operator+(difference_type n) const
    {
        return flat_linear_allocator_iterator(storage_, index_ + n, begin_, end_);
    }

    flat_linear_allocator_iterator &operator+=(difference_type n)
    {
        index_ += n;
        return *this;
    }

    flat_linear_allocator_iterator operator-(difference_type n) const
    {
        return flat_linear_allocator_iterator(storage_, index_ - n, begin_, end_);
    }
    flat_linear_allocator_iterator &operator-=(difference_type n)
    {
        index_ -= n;
        return *this;
    }

    difference_type operator-(const flat_linear_allocator_iterator &other) const
    {
        return static_cast<difference_type>(index_) - static_cast<difference_type>(other.index_);
    }

    reference            operator[](difference_type n) { return *(*this + n); }

    const_reference_type operator[](difference_type n) const { return *(*this + n); }

    bool                 operator<(const flat_linear_allocator_iterator &other) const
    {
        return index_ < other.index_;
    }

    bool operator>(const flat_linear_allocator_iterator &other) const
    {
        return index_ > other.index_;
    }

    bool operator<=(const flat_linear_allocator_iterator &other) const
    {
        return index_ <= other.index_;
    }

    bool operator>=(const flat_linear_allocator_iterator &other) const
    {
        return index_ >= other.index_;
    }

    bool operator!() const { return storage_ == nullptr; }
};
/* memory management helps for tape */
template<typename T, size_t buffer_size>
class flat_linear_allocator
{
    /** @brief basically a vector<array<T*, size>>
   * intended to work like a vector that allocates memory in chunks
   * and doesn't invalidate references
   * */
public:
    // store vector of unique pointers to arrays
    // to avoid vector reference invalidation
    using buffer_type = std::array<T, buffer_size>;
    using buffer_ptr  = std::unique_ptr<std::array<T, buffer_size>>;

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
    void allocate_buffer()
    {
        data_.emplace_back(std::make_unique<buffer_type>());
        // buffer_ptr new_buffer = new buffer_type();
        // data_.push_back(new_buffer);
    }

public:
    flat_linear_allocator() { allocate_buffer(); }
    flat_linear_allocator(const flat_linear_allocator &)            = delete;
    flat_linear_allocator &operator=(const flat_linear_allocator &) = delete;
    flat_linear_allocator(flat_linear_allocator &&)                 = delete;
    flat_linear_allocator &operator=(flat_linear_allocator &&)      = delete;
    ~flat_linear_allocator() { data_.clear(); }

    /** @brief
   * helper functions to clear tape and create block in tape
   */
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

    /** @brief clears all checkpoints
   * */
    void reset_checkpoints() { checkpoints_.clear(); }

    void rewind_to_last_checkpoint() { total_size_ = checkpoints_.back(); }
    void rewind_to_checkpoint_at(size_t index) { total_size_ = checkpoints_[index]; }

    void fill(const T &val)
    {
        for (size_t i = 0; i < total_size_; ++i) {
            size_t bid         = i / buffer_size;
            size_t iid         = i % buffer_size;
            (*data_[bid])[iid] = val;
        }
    }

    /** @brief emplaces back object at the end of the
   * data structure, calls default constructor */
    iterator emplace_back()
    {
        if (item_id() == 0 && total_size_ != 0) {
            allocate_buffer();
        }
        size_t bid = buffer_id();
        size_t iid = item_id();

        T     *ptr = &(*data_[bid])[iid];
        new (ptr) T();
        ++total_size_;
        return iterator(this, total_size_ - 1);
    };

    /** @brief, emplaces back object at end of data structure,
   * passes arguments to constructor */
    template<typename... Args>
    iterator emplace_back(Args &&...args)
    {
        if (item_id() == 0 && total_size_ != 0) {
            allocate_buffer();
        }
        assert(buffer_id() < data_.size());
        assert(item_id() < buffer_size);
        T *ptr = &(*data_[buffer_id()])[item_id()];
        new (ptr) T(std::forward<Args>(args)...);
        ++total_size_;
        return iterator(this, total_size_ - 1);
    };
    /** @brief default constructs n objects at end of
   * data structure, n known at compile time */
    template<size_t n>
    iterator emplace_back_n()
    {
        size_t bid = buffer_id();
        size_t iid = item_id();
        if (iid + n < buffer_size) {
            T *ptr = &(*data_[bid])[iid];
            for (size_t i = 0; i < n; ++i) {
                new (ptr + i) T();
            }
            total_size_ += n;
            return iterator(this, total_size_ - n, total_size_ - n, total_size_);
        } else {
            size_t allocs_in_curr_buffer = buffer_size - iid;
            size_t allocs_in_next_buffer = n - (buffer_size - iid);
            T     *ptr                   = &(*data_[bid])[iid];
            for (size_t i = 0; i < allocs_in_curr_buffer; ++i) {
                new (ptr + i) T();
            }
            allocate_buffer();
            bid          = data_.size() - 1;
            iid          = 0;
            total_size_ += n;

            T *ptr2      = &(*data_[bid])[iid];
            for (size_t i = 0; i < allocs_in_next_buffer; i++) {
                new (ptr2 + i) T();
            }
            return iterator(this, total_size_ - n, total_size_ - n, total_size_);
        }
    }
    /** @brief default constructs n objects at end of
   * data structure, n known at run time
   */
    iterator emplace_back_n(size_t n)
    {
        size_t bid = buffer_id();
        size_t iid = item_id();
        if (iid + n < buffer_size) {
            T *ptr = &(*data_[bid])[iid];
            for (size_t i = 0; i < n; ++i) {
                new (ptr + i) T();
            }
            total_size_ += n;
            return iterator(this, total_size_ - n, total_size_ - n, total_size_);
        } else {
            size_t allocs_in_curr_buffer = buffer_size - iid;
            size_t allocs_in_next_buffer = n - (buffer_size - iid);
            T     *ptr                   = &(*data_[bid])[iid];
            for (size_t i = 0; i < allocs_in_curr_buffer; ++i) {
                new (ptr + i) T();
            }
            allocate_buffer();
            bid          = data_.size() - 1;
            iid          = 0;
            total_size_ += n;

            T *ptr2      = &(*data_[bid])[iid];
            for (size_t i = 0; i < allocs_in_next_buffer; i++) {
                new (ptr2 + i) T();
            }
            return iterator(this, total_size_ - n, total_size_ - n, total_size_);
        }
    }

    /** @brief number of elements */
    size_t         size() const { return total_size_; }

    /** @brief total capacity */
    size_t         capacity() const { return data_.size() * buffer_size; }

    /** @brief iterator helpers */
    iterator       begin() { return iterator(this, 0); }
    iterator       end() { return iterator(this, total_size_); }
    const_iterator begin() const { return const_iterator(this, 0); }
    const_iterator end() const { return const_iterator(this, total_size_); }

    iterator       last_checkpoint() { return iterator(this, checkpoints_.back(), 0, total_size_); }
    iterator       first_checkpoint() { return iterator(this, checkpoints_[0], 0, total_size_); };
    iterator       checkpoint_at(size_t index)
    {
        return iterator(this, checkpoints_[index], 0, total_size_);
    };

    /** @brief searches for item in allocator
   *  only used to find gradient nodes for propagation */
    iterator find(const T *const item)
    {
        return std::find_if(begin(), end(), [&](const T &val) { return &val == item; });
    }
    /** @brief vector like access,
   *  currently unused anywhere but very useful for debugging
   */
    T &operator[](std::size_t i)
    {
        assert(i < total_size_);
        return (*data_[i / buffer_size])[i % buffer_size];
    }
    const T &operator[](std::size_t i) const
    {
        assert(i < total_size_);
        return (*data_[i / buffer_size])[i % buffer_size];
    }
};

/** @brief template metafunctions that
 *  construct rvar<T, N> if n > 0 and
 *  T if n == 0 */
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

/** @brief rvar_t<T,0> decays to T
 * otherwise its just identity rvar_t<T,N> = rvar<T,N>
 * */
template<typename T, size_t N>
using rvar_t = typename detail::rvar_type_impl<T, N>::type;
// manages nodes in computational graph
template<typename T, size_t order, size_t buffer_size = BUFFER_SIZE>
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
    // size_t depth() const override { return depth_; };

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

        if (!adjoint_ || fabs(*adjoint_) < std::numeric_limits<T>::epsilon()) //
                                                                              // zero
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
template<typename T, size_t order, class derived_expression>
struct expression
{
    /* @brief
   * base expression class
   * */
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
        , rhs(static_cast<const RHS &>(right_hand_expr)) {};

    inner_t evaluate() const
    {
        return static_cast<const concrete_binary_operation *>(this)->evaluate();
    };

    template<size_t arg_index>
    void propagatex(gradient_node<T, order> *node, inner_t adj) const
    {
        inner_t          lv           = lhs.evaluate();
        inner_t          rv           = rhs.evaluate();
        inner_t          v            = evaluate();
        inner_t          partial_l    = concrete_binary_operation::left_derivative(lv, rv, v);
        inner_t          partial_r    = concrete_binary_operation::right_derivative(lv, rv, v);

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
        , constant(constant) {};
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
    explicit add_expr(const expression<T, order, LHS> &left_hand_expr,
                      const expression<T, order, RHS> &right_hand_expr)
        : abstract_binary_expression<T, order, LHS, RHS, add_expr<T, order, LHS, RHS>>(
              left_hand_expr, right_hand_expr)
    {}

    inner_t              evaluate() const { return this->lhs.evaluate() + this->rhs.evaluate(); }
    static const inner_t left_derivative(const inner_t &l, const inner_t &r, const inner_t &v)
    {
        return inner_t(1.0);
    }
    static const inner_t right_derivative(const inner_t &l, const inner_t &r, const inner_t &v)
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
    using arg_type   = ARG;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;
    explicit add_const_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, add_const_expr<T, order, ARG>>(arg_expr, v) {};
    // explicit add_const_expr(const ARG& arg_expr, const T v)
    //     : abstract_unary_expression<T, order, ARG, add_const_expr<T, order,
    //     ARG>>(arg_expr, v)
    //{}
    inner_t              evaluate() const { return this->arg.evaluate() + inner_t(this->constant); }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
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
    explicit mult_expr(const expression<T, order, LHS> &left_hand_expr,
                       const expression<T, order, RHS> &right_hand_expr)
        : abstract_binary_expression<T, order, LHS, RHS, mult_expr<T, order, LHS, RHS>>(
              left_hand_expr, right_hand_expr)
    {}

    inner_t              evaluate() const { return this->lhs.evaluate() * this->rhs.evaluate(); };
    static const inner_t left_derivative(const inner_t &l, const inner_t &r, const inner_t &v)
    {
        return r;
    };
    static const inner_t right_derivative(const inner_t &l, const inner_t &r, const inner_t &v)
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
    using arg_type   = ARG;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;

    explicit mult_const_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, mult_const_expr<T, order, ARG>>(arg_expr, v) {};

    inner_t              evaluate() const { return this->arg.evaluate() * inner_t(this->constant); }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
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
    explicit sub_expr(const expression<T, order, LHS> &left_hand_expr,
                      const expression<T, order, RHS> &right_hand_expr)
        : abstract_binary_expression<T, order, LHS, RHS, sub_expr<T, order, LHS, RHS>>(
              left_hand_expr, right_hand_expr)
    {}

    inner_t              evaluate() const { return this->lhs.evaluate() - this->rhs.evaluate(); }
    static const inner_t left_derivative(const inner_t &l, const inner_t &r, const inner_t &v)
    {
        return inner_t(1.0);
    }
    static const inner_t right_derivative(const inner_t &l, const inner_t &r, const inner_t &v)
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
    explicit div_expr(const expression<T, order, LHS> &left_hand_expr,
                      const expression<T, order, RHS> &right_hand_expr)
        : abstract_binary_expression<T, order, LHS, RHS, div_expr<T, order, LHS, RHS>>(
              left_hand_expr, right_hand_expr)
    {}

    inner_t              evaluate() const { return this->lhs.evaluate() / this->rhs.evaluate(); };
    static const inner_t left_derivative(const inner_t &l, const inner_t &r, const inner_t &v)
    {
        return 1.0 / r;
    };
    static const inner_t right_derivative(const inner_t &l, const inner_t &r, const inner_t &v)
    {
        return -l / (r * r);
    };
};
template<typename T, size_t order, typename ARG>
struct div_by_const_expr
    : public abstract_unary_expression<T, order, ARG, div_by_const_expr<T, order, ARG>>
{
    /* @brief
   * rvar+float or float+rvar
   * */
    using arg_type   = ARG;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;

    explicit div_by_const_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, div_by_const_expr<T, order, ARG>>(arg_expr, v) {};

    inner_t              evaluate() const { return this->arg.evaluate() / inner_t(this->constant); }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
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
    using arg_type   = ARG;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;

    explicit const_div_by_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, const_div_by_expr<T, order, ARG>>(arg_expr, v) {};

    inner_t              evaluate() const { return inner_t(this->constant) / this->arg.evaluate(); }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        return -inner_t{constant} / (argv * argv);
    }
};
/****************************************************************************************************************/
/* stl support : expressions */

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
    using arg_type   = ARG;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;

    explicit fabs_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, fabs_expr<T, order, ARG>>(arg_expr, v) {};

    inner_t evaluate() const
    {
        using std::fabs;
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
    using arg_type   = ARG;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;

    explicit ceil_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, ceil_expr<T, order, ARG>>(arg_expr, v) {};

    inner_t evaluate() const
    {
        using std::ceil;
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
    using arg_type   = ARG;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;

    explicit floor_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, floor_expr<T, order, ARG>>(arg_expr, v) {};

    inner_t evaluate() const
    {
        using std::floor;
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
    using arg_type   = ARG;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;

    explicit trunc_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, trunc_expr<T, order, ARG>>(arg_expr, v) {};

    inner_t evaluate() const
    {
        using boost::math::trunc;
        using std::trunc;
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
    using arg_type   = ARG;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;

    explicit exp_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, exp_expr<T, order, ARG>>(arg_expr, v) {};

    inner_t evaluate() const
    {
        using std::exp;
        return exp(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
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
    using lhs_type   = LHS;
    using rhs_type   = RHS;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;
    // Explicitly define constructor to forward to base class
    explicit pow_expr(const expression<T, order, LHS> &left_hand_expr,
                      const expression<T, order, RHS> &right_hand_expr)
        : abstract_binary_expression<T, order, LHS, RHS, pow_expr<T, order, LHS, RHS>>(
              left_hand_expr, right_hand_expr)
    {}

    inner_t evaluate() const
    {
        using std::pow;
        return pow(this->lhs.evaluate(), this->rhs.evaluate());
    };
    static const inner_t left_derivative(const inner_t &l, const inner_t &r, const inner_t &v)
    {
        using std::pow;
        return r * pow(l, r - 1);
    };
    static const inner_t right_derivative(const inner_t &l, const inner_t &r, const inner_t &v)
    {
        using std::log;
        using std::pow;
        return pow(l, r) * log(l);
    };
};

template<typename T, size_t order, typename ARG>
struct expr_pow_float_expr
    : public abstract_unary_expression<T, order, ARG, expr_pow_float_expr<T, order, ARG>>
{
    /** @brief pow(rvar,float)
     *
      * */
    using arg_type   = ARG;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;

    explicit expr_pow_float_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, expr_pow_float_expr<T, order, ARG>>(arg_expr,
                                                                                       v) {};

    inner_t evaluate() const
    {
        using std::pow;
        return pow(this->arg.evaluate(), this->constant);
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        using std::pow;
        return inner_t{constant} * pow(argv, inner_t{constant - 1});
    }
};

template<typename T, size_t order, typename ARG>
struct float_pow_expr_expr
    : public abstract_unary_expression<T, order, ARG, float_pow_expr_expr<T, order, ARG>>
{
    /** @brief pow(float, rvar)
      * */
    using arg_type   = ARG;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;

    explicit float_pow_expr_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, float_pow_expr_expr<T, order, ARG>>(arg_expr,
                                                                                       v) {};

    inner_t evaluate() const
    {
        using std::pow;
        return pow(this->constant, this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        using std::pow;
        return pow(constant, argv) * log(constant);
    }
};

template<typename T, size_t order, typename ARG>
struct sqrt_expr : public abstract_unary_expression<T, order, ARG, sqrt_expr<T, order, ARG>>
{
    /** @brief  sqrt(x)
     *  d/dx sqrt(x) = 1/(2 sqrt(x))
    * */
    using arg_type   = ARG;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;

    explicit sqrt_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, sqrt_expr<T, order, ARG>>(arg_expr, v) {};

    inner_t evaluate() const
    {
        using std::sqrt;
        return sqrt(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        using std::sqrt;
        return 1.0 / (2 * sqrt(argv));
    }
};

template<typename T, size_t order, typename ARG>
struct log_expr : public abstract_unary_expression<T, order, ARG, log_expr<T, order, ARG>>
{
    /** @brief log(x)
     *  d/dx log(x) = 1/x
      * */
    using arg_type   = ARG;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;

    explicit log_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, log_expr<T, order, ARG>>(arg_expr, v) {};

    inner_t evaluate() const
    {
        using std::log;
        return log(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        return 1.0 / argv;
    }
};

template<typename T, size_t order, typename ARG>
struct cos_expr : public abstract_unary_expression<T, order, ARG, cos_expr<T, order, ARG>>
{
    /** @brief log(x)
     *  d/dx log(x) = 1/x
      * */
    using arg_type   = ARG;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;

    explicit cos_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, cos_expr<T, order, ARG>>(arg_expr, v) {};

    inner_t evaluate() const
    {
        using std::cos;
        return cos(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        using std::sin;
        return -sin(argv);
    }
};

template<typename T, size_t order, typename ARG>
struct sin_expr : public abstract_unary_expression<T, order, ARG, sin_expr<T, order, ARG>>
{
    /** @brief log(x)
     *  d/dx log(x) = 1/x
      * */
    using arg_type   = ARG;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;

    explicit sin_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, sin_expr<T, order, ARG>>(arg_expr, v) {};

    inner_t evaluate() const
    {
        using std::sin;
        return sin(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        using std::cos;
        return cos(argv);
    }
};

template<typename T, size_t order, typename ARG>
struct tan_expr : public abstract_unary_expression<T, order, ARG, tan_expr<T, order, ARG>>
{
    /** @brief log(x)
     *  d/dx log(x) = 1/x
      * */
    using arg_type   = ARG;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;

    explicit tan_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, tan_expr<T, order, ARG>>(arg_expr, v) {};

    inner_t evaluate() const
    {
        using std::tan;
        return tan(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        using std::cos;
        return 1.0 / (cos(argv) * cos(argv));
    }
};

template<typename T, size_t order, typename ARG>
struct acos_expr : public abstract_unary_expression<T, order, ARG, acos_expr<T, order, ARG>>
{
    /** @brief log(x)
     *  d/dx log(x) = 1/x
      * */
    using arg_type   = ARG;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;

    explicit acos_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, acos_expr<T, order, ARG>>(arg_expr, v) {};

    inner_t evaluate() const
    {
        using std::acos;
        return acos(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        using std::sqrt;
        return -1.0 / sqrt(1 - argv * argv);
    }
};

template<typename T, size_t order, typename ARG>
struct asin_expr : public abstract_unary_expression<T, order, ARG, asin_expr<T, order, ARG>>
{
    /** @brief log(x)
     *  d/dx log(x) = 1/x
      * */
    using arg_type   = ARG;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;

    explicit asin_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, asin_expr<T, order, ARG>>(arg_expr, v) {};

    inner_t evaluate() const
    {
        using std::asin;
        return asin(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        using std::sqrt;
        return 1.0 / sqrt(1 - argv * argv);
    }
};

template<typename T, size_t order, typename ARG>
struct atan_expr : public abstract_unary_expression<T, order, ARG, atan_expr<T, order, ARG>>
{
    /** @brief log(x)
     *  d/dx log(x) = 1/x
      * */
    using arg_type   = ARG;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;

    explicit atan_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, atan_expr<T, order, ARG>>(arg_expr, v) {};

    inner_t evaluate() const
    {
        using std::atan;
        return atan(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        return 1.0 / (1 + argv * argv);
    }
};
template<typename T, size_t order, typename LHS, typename RHS>
struct atan2_expr
    : public abstract_binary_expression<T, order, LHS, RHS, atan2_expr<T, order, LHS, RHS>>
{
    /** @brief pow(x,y)
     *  d/dx pow(x,y) = y pow (x, y-1)
     *  d/dy pow(x,y) = pow(x,y) log(x)
    * */
    using lhs_type   = LHS;
    using rhs_type   = RHS;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;
    // Explicitly define constructor to forward to base class
    explicit atan2_expr(const expression<T, order, LHS> &left_hand_expr,
                        const expression<T, order, RHS> &right_hand_expr)
        : abstract_binary_expression<T, order, LHS, RHS, atan2_expr<T, order, LHS, RHS>>(
              left_hand_expr, right_hand_expr)
    {}

    inner_t evaluate() const
    {
        using std::atan2;
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
    /** @brief pow(rvar,float)
     *
      * */
    using arg_type   = ARG;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;

    explicit atan2_left_float_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, atan2_left_float_expr<T, order, ARG>>(arg_expr,
                                                                                         v) {};

    inner_t evaluate() const
    {
        using std::atan2;
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
    /** @brief pow(float, rvar)
      * */
    using arg_type   = ARG;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;

    explicit atan2_right_float_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, atan2_right_float_expr<T, order, ARG>>(arg_expr,
                                                                                          v) {};

    inner_t evaluate() const
    {
        using std::atan2;
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
    /** @brief
     *  d/dx
      * */
    using arg_type   = ARG;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;

    explicit round_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, round_expr<T, order, ARG>>(arg_expr, v) {};

    inner_t evaluate() const
    {
        using boost::math::round;
        using std::round;
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
    /** @brief log(x)
     *  d/dx log(x) = 1/x
      * */
    using arg_type   = ARG;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;

    explicit sinh_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, sinh_expr<T, order, ARG>>(arg_expr, v) {};

    inner_t evaluate() const
    {
        using std::sinh;
        return sinh(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        using std::cosh;
        return cosh(argv);
    }
};

template<typename T, size_t order, typename ARG>
struct cosh_expr : public abstract_unary_expression<T, order, ARG, cosh_expr<T, order, ARG>>
{
    /** @brief log(x)
     *  d/dx log(x) = 1/x
      * */
    using arg_type   = ARG;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;

    explicit cosh_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, cosh_expr<T, order, ARG>>(arg_expr, v) {};

    inner_t evaluate() const
    {
        using std::cosh;
        return cosh(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        using std::sinh;
        return sinh(argv);
    }
};
template<typename T, size_t order, typename ARG>
struct tanh_expr : public abstract_unary_expression<T, order, ARG, tanh_expr<T, order, ARG>>
{
    /** @brief log(x)
     *  d/dx log(x) = 1/x
      * */
    using arg_type   = ARG;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;

    explicit tanh_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, tanh_expr<T, order, ARG>>(arg_expr, v) {};

    inner_t evaluate() const
    {
        using std::tanh;
        return tanh(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        using std::cosh;
        return 1.0 / (cosh(argv) * cosh(argv));
    }
};

template<typename T, size_t order, typename ARG>
struct log10_expr : public abstract_unary_expression<T, order, ARG, log10_expr<T, order, ARG>>
{
    /** @brief log(x)
     *  d/dx log(x) = 1/x
      * */
    using arg_type   = ARG;
    using value_type = T;
    using inner_t    = rvar_t<T, order - 1>;

    explicit log10_expr(const expression<T, order, ARG> &arg_expr, const T v)
        : abstract_unary_expression<T, order, ARG, log10_expr<T, order, ARG>>(arg_expr, v) {};

    inner_t evaluate() const
    {
        using std::log10;
        return log10(this->arg.evaluate());
    }
    static const inner_t derivative(const inner_t &argv, const inner_t &v, const T &constant)
    {
        using std::log;
        return 1.0 / (argv * log(10.0));
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
         typename = typename std::enable_if<std::is_arithmetic<U>::value>::type>
mult_const_expr<T, order, ARG> operator*(const expression<T, order, ARG> &arg, const U &v)
{
    return mult_const_expr<T, order, ARG>(arg, static_cast<T>(v));
}
template<typename U,
         typename T,
         size_t order,
         typename ARG,
         typename = typename std::enable_if<std::is_arithmetic<U>::value>::type>
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
         typename = typename std::enable_if<std::is_arithmetic<U>::value>::type>
add_const_expr<T, order, ARG> operator+(const expression<T, order, ARG> &arg, const U &v)
{
    return add_const_expr<T, order, ARG>(arg, static_cast<T>(v));
}
template<typename U,
         typename T,
         size_t order,
         typename ARG,
         typename = typename std::enable_if<std::is_arithmetic<U>::value>::type>
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
    return mult_const_expr<T, order, ARG>(arg, -1.0);
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
         typename = typename std::enable_if<std::is_arithmetic<U>::value>::type>
add_const_expr<T, order, ARG> operator-(const expression<T, order, ARG> &arg, const U &v)
{
    /* rvar - float = rvar + (-float) */
    return add_const_expr<T, order, ARG>(arg, static_cast<T>(-1.0 * v));
}

/** @brief
 *   subtraction float - rvar
 *  @return add_expr<neg_expr<ARG>>
 */
template<typename U,
         typename T,
         size_t order,
         typename ARG,
         typename = typename std::enable_if<std::is_arithmetic<U>::value>::type>
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
         typename = typename std::enable_if<std::is_arithmetic<U>::value>::type>
const_div_by_expr<T, order, ARG> operator/(const U &v, const expression<T, order, ARG> &arg)
{
    return const_div_by_expr<T, order, ARG>(arg, static_cast<T>(v));
}

template<typename U,
         typename T,
         size_t order,
         typename ARG,
         typename = typename std::enable_if<std::is_arithmetic<U>::value>::type>
div_by_const_expr<T, order, ARG> operator/(const expression<T, order, ARG> &arg, const U &v)
{
    return div_by_const_expr<T, order, ARG>(arg, static_cast<T>(v));
}
/****************************************************************************************************************/
template<typename T, size_t order_1, size_t order_2, class E, class F>
bool operator==(const expression<T, order_1, E> &lhs, const expression<T, order_2, F> &rhs)
{
    return lhs.evaluate() == rhs.evaluate();
}

template<typename U,
         typename T,
         size_t order,
         class E,
         typename = typename std::enable_if<std::is_arithmetic<U>::value>::type>
bool operator==(const expression<T, order, E> &lhs, const U &rhs)
{
    return lhs.evaluate() == static_cast<T>(rhs);
}
template<typename U,
         typename T,
         size_t order,
         class E,
         typename = typename std::enable_if<std::is_arithmetic<U>::value>::type>
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
         typename = typename std::enable_if<std::is_arithmetic<U>::value>::type>
bool operator!=(const expression<T, order, E> &lhs, const U &rhs)
{
    return lhs.evaluate() != rhs;
}
template<typename U,
         typename T,
         size_t order,
         class E,
         typename = typename std::enable_if<std::is_arithmetic<U>::value>::type>
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
         typename = typename std::enable_if<std::is_arithmetic<U>::value>::type>
bool operator<(const expression<T, order, E> &lhs, const U &rhs)
{
    return lhs.evaluate() < rhs;
}
template<typename U,
         typename T,
         size_t order,
         class E,
         typename = typename std::enable_if<std::is_arithmetic<U>::value>::type>
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
         typename = typename std::enable_if<std::is_arithmetic<U>::value>::type>
bool operator>(const expression<T, order, E> &lhs, const U &rhs)
{
    return lhs.evaluate() > rhs;
}
template<typename U,
         typename T,
         size_t order,
         class E,
         typename = typename std::enable_if<std::is_arithmetic<U>::value>::type>
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
         typename = typename std::enable_if<std::is_arithmetic<U>::value>::type>
bool operator<=(const expression<T, order, E> &lhs, const U &rhs)
{
    return lhs.evaluate() <= rhs;
}
template<typename U,
         typename T,
         size_t order,
         class E,
         typename = typename std::enable_if<std::is_arithmetic<U>::value>::type>
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
         typename = typename std::enable_if<std::is_arithmetic<U>::value>::type>
bool operator>=(const expression<T, order, E> &lhs, const U &rhs)
{
    return lhs.evaluate() >= rhs;
}
template<typename U,
         typename T,
         size_t order,
         class E,
         typename = typename std::enable_if<std::is_arithmetic<U>::value>::type>
bool operator>=(const U &lhs, const expression<T, order, E> &rhs)
{
    return lhs >= rhs.evaluate();
}
/****************************************************************************************************************/
/* stl ops */
template<typename T, size_t order, typename ARG>
fabs_expr<T, order, ARG> fabs(const expression<T, order, ARG> &arg)
{
    return fabs_expr<T, order, ARG>(arg, 0.0);
}
template<typename T, size_t order, typename ARG>
auto abs(const expression<T, order, ARG> &arg)
{
    return fabs(arg);
}
template<typename T, size_t order, typename ARG>
ceil_expr<T, order, ARG> ceil(const expression<T, order, ARG> &arg)
{
    return ceil_expr<T, order, ARG>(arg, 0.0);
}
template<typename T, size_t order, typename ARG>
floor_expr<T, order, ARG> floor(const expression<T, order, ARG> &arg)
{
    return floor_expr<T, order, ARG>(arg, 0.0);
}

template<typename T, size_t order, typename ARG>
exp_expr<T, order, ARG> exp(const expression<T, order, ARG> &arg)
{
    return exp_expr<T, order, ARG>(arg, 0.0);
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
         typename = typename std::enable_if<std::is_arithmetic<U>::value>::type>
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
    return log_expr<T, order, ARG>(arg, 0.0);
};

template<typename T, size_t order, typename ARG>
sqrt_expr<T, order, ARG> sqrt(const expression<T, order, ARG> &arg)
{
    return sqrt_expr<T, order, ARG>(arg, 0.0);
};

template<typename T, size_t order, typename ARG>
auto frexp(const expression<T, order, ARG> &arg, int *i)
{
    using std::frexp;
    using std::pow;
    T tmp = frexp(arg.evaluate(), i);
    return arg / pow(2.0, *i);
}

template<typename T, size_t order, typename ARG>
auto ldexp(const expression<T, order, ARG> &arg, const int &i)
{
    using std::ldexp;
    using std::pow;
    return arg * pow(2.0, i);
}

template<typename T, size_t order, typename ARG>
cos_expr<T, order, ARG> cos(const expression<T, order, ARG> &arg)
{
    return cos_expr<T, order, ARG>(arg, 0.0);
};

template<typename T, size_t order, typename ARG>
sin_expr<T, order, ARG> sin(const expression<T, order, ARG> &arg)
{
    return sin_expr<T, order, ARG>(arg, 0.0);
};

template<typename T, size_t order, typename ARG>
tan_expr<T, order, ARG> tan(const expression<T, order, ARG> &arg)
{
    return tan_expr<T, order, ARG>(arg, 0.0);
};

template<typename T, size_t order, typename ARG>
acos_expr<T, order, ARG> acos(const expression<T, order, ARG> &arg)
{
    return acos_expr<T, order, ARG>(arg, 0.0);
};

template<typename T, size_t order, typename ARG>
asin_expr<T, order, ARG> asin(const expression<T, order, ARG> &arg)
{
    return asin_expr<T, order, ARG>(arg, 0.0);
};

template<typename T, size_t order, typename ARG>
atan_expr<T, order, ARG> atan(const expression<T, order, ARG> &arg)
{
    return atan_expr<T, order, ARG>(arg, 0.0);
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
    return trunc_expr<T, order, ARG>(arg, 0.0);
}

template<typename T, size_t order, typename LHS, typename RHS>
auto fmod(const expression<T, order, LHS> &lhs, const expression<T, order, RHS> &rhs)
{
    return lhs - trunc(lhs / rhs) * rhs;
}

template<typename T, size_t order, typename ARG>
auto fmod(const expression<T, order, ARG> &lhs, const T rhs)
{
    return lhs - trunc(lhs / rhs) * rhs;
}

template<typename T, size_t order, typename ARG>
auto fmod(const T lhs, const expression<T, order, ARG> &rhs)
{
    return lhs - trunc(lhs / rhs) * rhs;
}

template<typename T, size_t order, typename ARG>
round_expr<T, order, ARG> round(const expression<T, order, ARG> &arg)
{
    return round_expr<T, order, ARG>(arg, 0.0);
}

template<typename T, size_t order, typename ARG>
int iround(const expression<T, order, ARG> &arg)
{
    using boost::math::iround;
    rvar<T, order> tmp = arg.evaluate();
    return iround(tmp.item());
}
template<typename T, size_t order, typename ARG>
long lround(const expression<T, order, ARG> &arg)
{
    using boost::math::lround;
    rvar<T, order> tmp = arg.evaluate();
    return iround(tmp.item());
}

template<typename T, size_t order, typename ARG>
long long llround(const expression<T, order, ARG> &arg)
{
    using boost::math::llround;
    rvar<T, order> tmp = arg.evaluate();
    return iround(tmp.item());
}

template<typename T, size_t order, typename ARG>
int itrunc(const expression<T, order, ARG> &arg)
{
    using boost::math::itrunc;
    rvar<T, order> tmp = arg.evaluate();
    return itrunc(tmp.item());
}

template<typename T, size_t order, typename ARG>
long ltrunc(const expression<T, order, ARG> &arg)
{
    using boost::math::ltrunc;
    rvar<T, order> tmp = arg.evaluate();
    return ltrunc(tmp.item());
}

template<typename T, size_t order, typename ARG>
long long lltrunc(const expression<T, order, ARG> &arg)
{
    using boost::math::lltrunc;
    rvar<T, order> tmp = arg.evaluate();
    return lltrunc(tmp.item());
}

template<typename T, size_t order, typename ARG>
sinh_expr<T, order, ARG> sinh(const expression<T, order, ARG> &arg)
{
    return sinh_expr<T, order, ARG>(arg, 0.0);
}

template<typename T, size_t order, typename ARG>
cosh_expr<T, order, ARG> cosh(const expression<T, order, ARG> &arg)
{
    return cosh_expr<T, order, ARG>(arg, 0.0);
}

template<typename T, size_t order, typename ARG>
tanh_expr<T, order, ARG> tanh(const expression<T, order, ARG> &arg)
{
    return tanh_expr<T, order, ARG>(arg, 0.0);
}

template<typename T, size_t order, typename ARG>
log10_expr<T, order, ARG> log10(const expression<T, order, ARG> &arg)
{
    return log10_expr<T, order, ARG>(arg, 0.0);
}
/****************************************************************************************************************/
template<typename T, size_t order>
inline gradient_tape<T, order, BUFFER_SIZE> &get_active_tape()
{
    static thread_local gradient_tape<T, order, BUFFER_SIZE> tape;
    return tape;
}
/****************************************************************************************************************/
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
        gradient_tape<T, order, BUFFER_SIZE> &tape = get_active_tape<T, order>();
        node_                                      = tape.emplace_leaf_node();
    }

    void make_unary_node()
    {
        gradient_tape<T, order, BUFFER_SIZE> &tape = get_active_tape<T, order>();
        node_                                      = tape.emplace_active_unary_node();
    }

    void make_multi_node(size_t n)
    {
        gradient_tape<T, order, BUFFER_SIZE> &tape = get_active_tape<T, order>();
        node_                                      = tape.emplace_active_multi_node(n);
    }

    template<size_t n>
    void make_multi_node()
    {
        gradient_tape<T, order, BUFFER_SIZE> &tape = get_active_tape<T, order>();
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

    explicit       operator T() const { return item(); }

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
        gradient_tape<T, order, BUFFER_SIZE> &tape = get_active_tape<T, order>();
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
 *  		specialization for autodiffing through autodiff. i.e. being able to
 *  		compute higher order grads
*/
template<typename T, size_t order>
struct grad_op_impl
{
    std::vector<rvar_t<T, order - 1> *> operator()(rvar<T, order>                &f,
                                                   std::vector<rvar<T, order> *> &x)
    {
        gradient_tape<T, order, BUFFER_SIZE> &tape = get_active_tape<T, order>();
        tape.zero_grad();
        f.backward();
        std::vector<rvar_t<T, order - 1> *> gradient_vector;
        for (auto xi : x) {
            gradient_vector.push_back(&(xi->adjoint()));
        }
        return gradient_vector;
    }
};
/** @brief helper overload for grad implementation.
 *  @return vector<T> of gradients of the autodiff graph.
 *  		base specialization for order 1 autodiff
*/
template<typename T>
struct grad_op_impl<T, 1>
{
    std::vector<T> operator()(rvar<T, 1> &f, std::vector<rvar<T, 1> *> &x)
    {
        gradient_tape<T, 1, BUFFER_SIZE> &tape = get_active_tape<T, 1>();
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
 *  	higher order derivatives
 */
template<size_t N, typename T, size_t order_1, size_t order_2, typename Enable = void>
struct grad_nd_impl
{
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
    }
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
 * 			are computed simultaneously, only the ones w.r.t. x are returned
 * @return vector<rvar<T,order_1 - 1> of gradients. in the case of order_1 = 1
 * 		   rvar<T,order_1-1> decays to T
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
 *  	   vector of variables x.
 *  @return std::vector<std::vector<rvar<T,order_1-2>> hessian matrix
 *  		rvar<T,2> decays to T
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
 *  		of gradients of order N
 *
 *  NOT recursively safe, cannot do grad_nd(grad_nd(... etc..
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
