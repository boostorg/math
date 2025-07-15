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
#include <map>
#include <memory>
#include <stack>
#include <type_traits>
#include <vector>
#define BUFFER_SIZE 65536

namespace boost {
namespace math {
namespace differentiation {
namespace reverse_mode {
namespace detail {
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
    using ptr_type             = typename allocator_type::value_type*;
    using const_ptr_type       = const value_type*;
    using reference_type       = typename allocator_type::value_type&;
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
    reference_type operator*()
    {
        assert(index_ >= begin_ && index_ < end_);
        return (*storage_->data_[index_ / buffer_size])[index_ % buffer_size];
    }

    const_reference_type operator*() const
    {
        assert(index_ >= begin_ && index_ < end_);
        return (*storage_->data_[index_ / buffer_size])[index_ % buffer_size];
    }

    ptr_type operator->()
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

    reference_type operator[](difference_type n) { return *(*this + n); }

    const_reference_type operator[](difference_type n) const { return *(*this + n); }

    bool                 operator<(const flat_linear_allocator_iterator& other) const
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

    // resets checkpoint
    void reset_checkpoints()
    {
        checkpoints_.clear();
    }

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

        T*     ptr = &(*data_[bid])[iid];
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

    std::size_t    size() const { return total_size_; }
    std::size_t    capacity() const { return data_.size() * buffer_size; }

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

    /*
    iterator       find(const T* const item)
    {
        iterator start = begin();
        iterator it    = end();
        while (it != start) {
            --it;
            if (it.operator->() == item) {
                return it;
            }
        }
        if (&*it == item) {
            return it;
        };
        return end();
    }
	*/
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
template<typename T>
class gradient_node; // forward declaration for tape

// base class to be able to store tapes for different T's
// in case of higher order autodiff
class gradient_tape_base
{
public:
    virtual ~gradient_tape_base() = default;
    virtual size_t depth() const  = 0;
};

// manages nodes in computational graph
template<typename T, size_t buffer_size = BUFFER_SIZE>
class gradient_tape : public gradient_tape_base
{
private:
    size_t                                               depth_ = 0;
    detail::flat_linear_allocator<T, buffer_size>                adjoints_;
    detail::flat_linear_allocator<T, buffer_size>                derivatives_;
    detail::flat_linear_allocator<T*, buffer_size>               argument_adjoints_;
    detail::flat_linear_allocator<gradient_node<T>, buffer_size> gradient_nodes_;

    // compile time check if emplace_back_multi calls on zero
    template<size_t n>
    gradient_node<T>* fill_node_at_compile_time(std::true_type, gradient_node<T>* node_ptr)
    {
        node_ptr->derivatives_       = &*derivatives_.template emplace_back_n<n>();
        node_ptr->argument_adjoints_ = &*argument_adjoints_.template emplace_back_n<n>();
        return node_ptr;
    }

    template<size_t n>
    gradient_node<T>* fill_node_at_compile_time(std::false_type, gradient_node<T>* node_ptr)
    {
        node_ptr->derivatives_       = nullptr;
        node_ptr->argument_adjoints_ = nullptr;
        return node_ptr;
    }

public:
    gradient_tape(size_t depth)
        : depth_(depth)
    {
        clear();
    };
    size_t depth() const override { return depth_; };

    gradient_tape(const gradient_tape&)            = delete;
    gradient_tape& operator=(const gradient_tape&) = delete;
    gradient_tape(gradient_tape&& other)
        : adjoints_(std::move(other.adjoints_))
        , derivatives_(std::move(other.derivatives_))
        , argument_adjoints_(std::move(other.argument_adjoints_))
        , gradient_nodes_(std::move(other.gradient_nodes_))
    {
        other.clear();
    }
    gradient_tape operator=(gradient_tape&& other)
    {
        if (this != &other) {
            adjoints_          = std::move(other.adjoints_);
            derivatives_       = std::move(other.derivatives_);
            argument_adjoints_ = std::move(other.argument_adjoints_);
            gradient_nodes_    = std::move(other.gradient_nodes_);
        }
        return *this;
    }
    void clear()
    {
        adjoints_.clear();
        derivatives_.clear();
        argument_adjoints_.clear();
        gradient_nodes_.clear();
    }

    // no derivatives or arguments
    gradient_node<T>* emplace_leaf_node()
    {
        gradient_node<T>* node   = &*gradient_nodes_.emplace_back();
        node->adjoint_           = &*adjoints_.emplace_back();
        node->derivatives_       = nullptr;
        node->argument_adjoints_ = nullptr;
        return node;
    };

    // single argument, single derivative
    gradient_node<T>* emplace_active_unary_node()
    {
        gradient_node<T>* node   = &*gradient_nodes_.emplace_back();
        node->n_                 = 1;
        node->adjoint_           = &*adjoints_.emplace_back();
        node->derivatives_       = &*derivatives_.emplace_back();
        node->argument_adjoints_ = &*argument_adjoints_.emplace_back();
        return node;
    };

    // arbitrary number of arguments/derivatives (compile time)
    template<size_t n>
    gradient_node<T>* emplace_active_multi_node()
    {
        gradient_node<T>* node = &*gradient_nodes_.emplace_back();
        node->n_               = n;
        node->adjoint_         = &*adjoints_.emplace_back();
        // emulate if constexpr
        return fill_node_at_compile_time<n>(std::integral_constant<bool, (n > 0)>{}, node);
    };

    // same as above at runtime
    gradient_node<T>* emplace_active_multi_node(size_t n)
    {
        gradient_node<T>* node = &*gradient_nodes_.emplace_back();
        node->n_               = n;
        node->adjoint_         = &*adjoints_.emplace_back();
        if (n > 0) {
            node->derivatives_       = &*derivatives_.emplace_back_multi(n);
            node->argument_adjoints_ = &*argument_adjoints_.emplace_back_multi(n);
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
    auto find(gradient_node<T>* node) { return gradient_nodes_->find(node); };
    void add_checkpoint()
    {
        gradient_nodes_.add_checkpoint();
        adjoints_.add_checkpoint();
        derivatives_.add_checkpoint();
        argument_adjoints_.add_checkpoint();
    };

    auto last_checkpoint() { return gradient_nodes_.last_checkpoint(); };
    auto first_checkpoint() { return gradient_nodes_.last_checkpoint(); };
    auto checkpoint_at(size_t index) { return gradient_nodes_.get_checkpoint_at(index); };
    void rewind_to_last_checkpoint()
    {
        gradient_nodes_.rewind_to_last_checkpoint();
        adjoints_.rewind_to_last_checkpoint();
        derivatives_.rewind_to_last_checkpoint();
        argument_adjoints_.rewind_to_last_checkpoint();
    };
    void rewind_to_checkpoint_at(
        size_t index) // index is "checkpoint" index. so order which checkpoint was set
    {
        gradient_nodes_.rewind_to_checkpoint_at(index);
        adjoints_.rewind_to_checkpoint_at(index);
        derivatives_.rewind_to_checkpoint_at(index);
        argument_adjoints_.rewind_to_checkpoint_at(index);
    }

    // rewind to beginning of computational graph
    void rewind()
    {
        gradient_nodes_.rewind();
        adjoints_.rewind();
        derivatives_.rewind();
        argument_adjoints_.rewind();
    }

    // randoma acces
    T&       operator[](std::size_t i) { return gradient_nodes_[i]; }
    const T& operator[](std::size_t i) const { return gradient_nodes_[i]; }
};

enum arg_index { lhs_id = 0, rhs_id = 1, arg_id = 0 };

template<typename T> // no CRTP, just storage, maybe needs size_t template in the future
class gradient_node
{
    /*
     * @brief manages adjoints, derivatives, and stores points to argument adjoints
     * pointers to arguments aren't needed here
     * */
public:
    size_t n_;
    T*     adjoint_;
    T*     derivatives_;
    T**    argument_adjoints_;

public:
    friend class gradient_tape<T>;
    gradient_node() = default;
    explicit gradient_node(const size_t n)
        : n_(n)
        , adjoint_(nullptr)
        , derivatives_(nullptr)
        , argument_adjoints_(nullptr){};
    explicit gradient_node(const size_t n, T* adjoint, T* derivatives, T** argument_adjoints)
        : n_(n)
        , adjoint_(adjoint)
        , derivatives_(derivatives)
        , argument_adjoints_(argument_adjoints){};

    void update_adjoint_value(T value) { *adjoint_ = value; };
    void update_derivative_value_at(size_t id, T value) { derivatives_[id] = value; };

    void set_adjoint_pointer(T* adj) { adjoint_ = adj; };
    void set_argument_adjoint_pointer_at(size_t id, T* adj_ptr)
    {
        argument_adjoints_[id] = adj_ptr;
    };
    T adjoint() const { return *adjoint_; }
    T derivative() const { return derivatives_[arg_index::arg_id]; };
    T lhs_derivative() const { return derivatives_[arg_index::lhs_id]; };
    T rhs_derivative() const { return derivatives_[arg_index::rhs_id]; };
    T derivative_at(size_t id) const { return derivatives_[id]; };

    void backward()
    {
        if (!n_)
            return;

        if (!adjoint_)
            return;

        if (!argument_adjoints_)
            return;

        for (size_t i = 0; i < n_; ++i) {
            *(argument_adjoints_[i]) += derivatives_[i] * (*adjoint_);
        }
    }
};
template<typename T, class derived_expression>
struct expression
{
    /* @brief
     * base expression class
     * */
    static constexpr size_t num_literals = 0;
    T evaluate() const { return static_cast<const derived_expression*>(this)->evaluate(); }
    template<size_t literal_index>
    void propagatex(gradient_node<T>* node, T adj) const
    {
        return static_cast<const derived_expression*>(this)->template propagatex<literal_index>(node,
                                                                                                adj);
    }
};

template<typename T, typename LHS, typename RHS, typename concrete_binary_operation>
struct abstract_binary_expression
    : public expression<T, abstract_binary_expression<T, LHS, RHS, concrete_binary_operation>>
{
    const LHS               lhs;
    const RHS               rhs;
    static constexpr size_t num_literals = LHS::num_literals + RHS::num_literals;

    explicit abstract_binary_expression(const expression<T, LHS>& left_hand_expr,
                                        const expression<T, RHS>& right_hand_expr)
        : lhs(static_cast<const LHS&>(left_hand_expr))
        , rhs(static_cast<const RHS&>(right_hand_expr)){};

    T evaluate() const { return static_cast<const concrete_binary_operation*>(this)->evaluate(); };

    template<size_t literal_index>
    void propagatex(gradient_node<T>* node, T adj) const
    {
        T lv        = lhs.evaluate();
        T rv        = rhs.evaluate();
        T v         = evaluate();
        T partial_l = concrete_binary_operation::left_derivative(lv, rv, v);
        T partial_r = concrete_binary_operation::right_derivative(lv, rv, v);
        if (LHS::num_literals) {
            lhs.template propagatex<literal_index>(node, adj * partial_l);
        }
        if (RHS::num_literals) {
            rhs.template propagatex<literal_index + LHS::num_literals>(node, adj * partial_r);
        }
    }
};
template<typename T, typename LHS, typename RHS>
struct add_expr : public abstract_binary_expression<T, LHS, RHS, add_expr<T, LHS, RHS>>
{
    // Explicitly define constructor to forward to base class
    explicit add_expr(const expression<T, LHS>& left_hand_expr,
                      const expression<T, RHS>& right_hand_expr)
        : abstract_binary_expression<T, LHS, RHS, add_expr<T, LHS, RHS>>(left_hand_expr,
                                                                         right_hand_expr)
    {}

    T              evaluate() const { return this->lhs.evaluate() + this->rhs.evaluate(); }
    static const T left_derivative(const T& l, const T& r, const T& v) { return T(1.0); }
    static const T right_derivative(const T& l, const T& r, const T& v) { return T(1.0); }
};

template<typename T, typename LHS, typename RHS>
struct mult_expr : public abstract_binary_expression<T, LHS, RHS, mult_expr<T, LHS, RHS>>
{
    // Explicitly define constructor to forward to base class
    explicit mult_expr(const expression<T, LHS>& left_hand_expr,
                       const expression<T, RHS>& right_hand_expr)
        : abstract_binary_expression<T, LHS, RHS, mult_expr<T, LHS, RHS>>(left_hand_expr,
                                                                          right_hand_expr)
    {}

    T              evaluate() const { return this->lhs.evaluate() * this->rhs.evaluate(); };
    static const T left_derivative(const T& l, const T& r, const T& v) { return r; };
    static const T right_derivative(const T& l, const T& r, const T& v) { return l; };
};

template<typename T, typename LHS, typename RHS>
mult_expr<T, LHS, RHS> operator*(const expression<T, LHS>& lhs, const expression<T, RHS>& rhs)
{
    return mult_expr<T, LHS, RHS>(lhs, rhs);
}
template<typename T, typename LHS, typename RHS>
add_expr<T, LHS, RHS> operator+(const expression<T, LHS>& lhs, const expression<T, RHS>& rhs)
{
    return add_expr<T, LHS, RHS>(lhs, rhs);
}

template<typename T>
class rvar;          // forward declaration
template<typename T> //counts variable depth, so rvar<double> = 0, rvar<rvar<double>> = 1 ...
struct type_depth
{
    static constexpr size_t value = 0;
};

template<typename T>
struct type_depth<rvar<T>>
{
    static constexpr size_t value = type_depth<T>::value + 1;
};

thread_local std::map<size_t, std::stack<gradient_tape_base*>> active_tapes; // tape manager

inline gradient_tape_base* get_active_tape(size_t depth) // returns currently active tape
{
    auto it = active_tapes.find(depth);
    if (it != active_tapes.end() && !it->second.empty()) {
        return it->second.top();
    }
    return nullptr;
}
class scoped_tape_context // tape context
{
private:
    size_t depth;

public:
    scoped_tape_context(gradient_tape_base& tape)
        : depth(tape.depth())
    {
        active_tapes[depth].push(&tape);
    }
    ~scoped_tape_context() { active_tapes[depth].pop(); }
};

template<typename T>
class rvar : public expression<T, rvar<T>>
{
private:
    friend class gradient_node<T>;
    T  value_;
    gradient_node<T>*       node_;
    static constexpr size_t depth_ = type_depth<T>::value;

    template<typename U>
    U* safe_cast(gradient_tape_base* base)
    {
        return dynamic_cast<U*>(base);
    }

    void make_leaf_node()
    {
        if (gradient_tape<T, BUFFER_SIZE>* tape = safe_cast<gradient_tape<T, BUFFER_SIZE>>(
                get_active_tape(depth_))) {
            node_ = tape->emplace_leaf_node();
        }
    }

    void make_unary_node()
    {
        if (gradient_tape<T, BUFFER_SIZE>* tape = safe_cast<gradient_tape<T, BUFFER_SIZE>>(
                get_active_tape(depth_))) {
            node_ = tape->emplace_active_unary_node();
        }
    }

    void make_multi_node(size_t n)
    {
        if (gradient_tape<T, BUFFER_SIZE>* tape = safe_cast<gradient_tape<T, BUFFER_SIZE>>(
                get_active_tape(depth_))) {
            node_ = tape->emplace_active_multi_node(n);
        }
    }

    template<size_t n>
    void make_multi_node()
    {
        if (gradient_tape<T, BUFFER_SIZE>* tape = safe_cast<gradient_tape<T, BUFFER_SIZE>>(
                get_active_tape(depth_))) {
            node_ = tape->template emplace_active_multi_node<n>();
        }
    }

    template<typename E>
    void make_rvar(expression<T, E>& expr)
    {
        make_multi_node<E::num_literals>();
        expr.template propagatex<0>(node_, T(1.0));
    }

public:
    static constexpr size_t num_literals = 1;
    rvar()
        : value_(0.0)
    {
        make_leaf_node();
    }
    explicit rvar(const T value)
        : value_(value)
    {
        make_leaf_node();
    }
    rvar operator=(const T value)
    {
        value_ = value;
        make_leaf_node();
    }

    template<size_t literal_index>
    void propagatex(gradient_node<T>* node, T adj) const
    {
        //node->update_derivative_value_at(literal_index, adj);
        //node->set_argument_adjoint_pointer_at(literal_index, node_->adjoint_);
        node->derivatives_[literal_index] = adj;
    }
    template<class E>
    rvar(expression<T, E>& expr)
    {
        value_ = expr.evaluate();
        make_rvar(expr);
    }

    T adjoint() const { return node_->adjoint(); };
    T evaluate() const { return value_; }; // TODO: these will be different
    T item() const { return value_; };

    void print_der_info()
    {
        for (size_t i = 0; i < node_->n_; i++) {
            std::cout << node_->derivatives_[i] << std::endl;
        }
    }
};
} // namespace reverse_mode
} // namespace differentiation
} // namespace math
} // namespace boost

#endif
