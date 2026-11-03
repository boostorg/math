//           Copyright Maksym Zhelyeznyakov 2025-2026.
// Distributed under the Boost Software License, Version 1.0.
//      (See accompanying file LICENSE_1_0.txt or copy at
//           https://www.boost.org/LICENSE_1_0.txt)

#ifndef REVERSE_MODE_AUTODIFF_MEMORY_MANAGEMENT_HPP
#define REVERSE_MODE_AUTODIFF_MEMORY_MANAGEMENT_HPP

#ifndef BOOST_MATH_BUILD_MODULE
#include <algorithm>
#include <array>
#endif
#include <boost/math/tools/assert.hpp>
#ifndef BOOST_MATH_BUILD_MODULE
#include <cassert>
#include <cstddef>
#include <iterator>
#include <memory>
#include <type_traits>
#include <vector>
#endif
BOOST_MATH_NAMESPACE_BEGIN
namespace differentiation {
namespace reverse_mode {
namespace detail {
template<typename allocator_type, size_t buffer_size>
class flat_linear_allocator_iterator
{
  /**
   * @brief enables iterating over linear allocator with
   * c++ iterators
   */
public:
  using raw_allocator_type = std::remove_const_t<allocator_type>;
  using value_type = typename allocator_type::value_type;
  using pointer = typename allocator_type::value_type*;
  using const_ptr_type = const value_type*;
  using reference = typename allocator_type::value_type&;
  using const_reference_type = const value_type&;
  using iterator_category = std::random_access_iterator_tag;
  using difference_type = ptrdiff_t;

private:
  const allocator_type* storage_ = nullptr;
  size_t index_ = 0;
  size_t begin_ = 0;
  size_t end_ = 0;

public:
  flat_linear_allocator_iterator() = default;

  explicit flat_linear_allocator_iterator(allocator_type* storage, size_t index)
    : storage_(storage)
    , index_(index)
    , begin_(0)
    , end_(storage->size())
  {
  }

  explicit flat_linear_allocator_iterator(allocator_type* storage,
                                          size_t index,
                                          size_t begin,
                                          size_t end)
    : storage_(storage)
    , index_(index)
    , begin_(begin)
    , end_(end)
  {
  }

  explicit flat_linear_allocator_iterator(const allocator_type* storage,
                                          size_t index)
    : storage_(storage)
    , index_(index)
    , begin_(0)
    , end_(storage->size())
  {
  }

  explicit flat_linear_allocator_iterator(const allocator_type* storage,
                                          size_t index,
                                          size_t begin,
                                          size_t end)
    : storage_(storage)
    , index_(index)
    , begin_(begin)
    , end_(end)
  {
  }
  reference operator*()
  {
    BOOST_MATH_ASSERT(index_ >= begin_ && index_ < end_);
    return (*storage_->data_[index_ / buffer_size])[index_ % buffer_size];
  }

  const_reference_type operator*() const
  {
    BOOST_MATH_ASSERT(index_ >= begin_ && index_ < end_);
    return (*storage_->data_[index_ / buffer_size])[index_ % buffer_size];
  }

  pointer operator->()
  {
    BOOST_MATH_ASSERT(index_ >= begin_ && index_ < end_);
    return &operator*();
  }

  const_ptr_type operator->() const
  {
    BOOST_MATH_ASSERT(index_ >= begin_ && index_ < end_);
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

  bool operator!=(const flat_linear_allocator_iterator& other) const
  {
    return !(*this == other);
  }

  flat_linear_allocator_iterator operator+(difference_type n) const
  {
    return flat_linear_allocator_iterator(
      storage_, index_ + static_cast<size_t>(n), begin_, end_);
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
    return static_cast<difference_type>(index_) -
           static_cast<difference_type>(other.index_);
  }

  reference operator[](difference_type n) { return *(*this + n); }

  const_reference_type operator[](difference_type n) const
  {
    return *(*this + n);
  }

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

  bool operator!() const noexcept { return storage_ == nullptr; }
};
/* memory management helps for tape */
template<typename RealType, size_t buffer_size>
class flat_linear_allocator
{
  /** @brief basically a vector<array<T*, size>>
   * intended to work like a vector that allocates memory in chunks
   * and doesn't invalidate references
   * */
public:
  // store vector of unique pointers to arrays
  // to avoid vector reference invalidation
  // Uninitialized storage: elements are constructed by emplace_back. Default-constructing a whole block
  // would, for the rvar<T, 1> elements of a second-order tape, put buffer_size leaves on the first-order tape.
  struct buffer_type
  {
    alignas(RealType) unsigned char bytes[sizeof(RealType) * buffer_size];
    RealType& operator[](size_t i) noexcept { return *reinterpret_cast<RealType*>(bytes + sizeof(RealType) * i); }
    const RealType& operator[](size_t i) const noexcept { return *reinterpret_cast<const RealType*>(bytes + sizeof(RealType) * i); }
  };
  using buffer_ptr = std::unique_ptr<buffer_type>;

private:
  std::vector<buffer_ptr> data_;
  size_t total_size_ = 0;
  std::vector<size_t> checkpoints_; //{0};

public:
  friend class flat_linear_allocator_iterator<
    flat_linear_allocator<RealType, buffer_size>,
    buffer_size>;
  friend class flat_linear_allocator_iterator<
    const flat_linear_allocator<RealType, buffer_size>,
    buffer_size>;
  using value_type = RealType;
  using iterator =
    flat_linear_allocator_iterator<flat_linear_allocator<RealType, buffer_size>,
                                   buffer_size>;
  using const_iterator = flat_linear_allocator_iterator<
    const flat_linear_allocator<RealType, buffer_size>,
    buffer_size>;

  size_t buffer_id() const noexcept { return total_size_ / buffer_size; }
  size_t item_id() const noexcept { return total_size_ % buffer_size; }

private:
  // The storage for element i; blocks are allocated on first use and kept, so a rewound tape reuses them.
  RealType* slot(size_t i)
  {
    while (i / buffer_size >= data_.size()) {
      data_.emplace_back(buffer_ptr(new buffer_type));
    }
    return &(*data_[i / buffer_size])[i % buffer_size];
  }
  // Destroys the elements from n on; their blocks are kept. Checkpoints past n are dropped.
  void truncate(size_t n)
  {
    BOOST_MATH_ASSERT(n <= total_size_);
    while (!checkpoints_.empty() && checkpoints_.back() > n) {
      checkpoints_.pop_back();
    }
    for (size_t i = n; i < total_size_; ++i) {
      (*this)[i].~RealType();
    }
    total_size_ = n;
  }

public:
  flat_linear_allocator() { slot(0); }
  flat_linear_allocator(const flat_linear_allocator&) = delete;
  flat_linear_allocator& operator=(const flat_linear_allocator&) = delete;
  flat_linear_allocator(flat_linear_allocator&&) = delete;
  flat_linear_allocator& operator=(flat_linear_allocator&&) = delete;
  ~flat_linear_allocator()
  {
    destroy_all();
    data_.clear();
  }

  void destroy_all()
  {
    for (size_t i = 0; i < total_size_; ++i) {
      size_t bid = i / buffer_size;
      size_t iid = i % buffer_size;
      (*data_[bid])[iid].~RealType();
    }
  }
  /** @brief
   * helper functions to clear tape and create block in tape
   */
  void clear()
  {
    destroy_all();
    data_.clear();
    total_size_ = 0;
    checkpoints_.clear();
    slot(0);
  }

  // doesn't delete anything, only sets the current index to zero
  void reset() { truncate(0); }
  void rewind() { truncate(0); };

  // adds current index as a checkpoint to be able to walk back to
  void add_checkpoint()
  {
    checkpoints_.push_back(total_size_);
    /*
     * if (total_size_ > 0) {
      checkpoints_.push_back(total_size_ - 1);
    } else {
      checkpoints_.push_back(0);
    }*/
  };

  /** @brief clears all checkpoints
   * */
  void reset_checkpoints() { checkpoints_.clear(); }

  void rewind_to_last_checkpoint() { truncate(checkpoints_.back()); }
  void rewind_to_checkpoint_at(size_t index)
  {
    truncate(checkpoints_[index]);
  }

  void fill(const RealType& val)
  {
    for (size_t i = 0; i < total_size_; ++i) {
      size_t bid = i / buffer_size;
      size_t iid = i % buffer_size;
      (*data_[bid])[iid] = val;
    }
  }

  /** @brief emplaces back object at end of data structure,
   * passes arguments to constructor */
  template<typename... Args>
  iterator emplace_back(Args&&... args)
  {
    new (slot(total_size_)) RealType(std::forward<Args>(args)...);
    ++total_size_;
    return iterator(this, total_size_ - 1);
  }
  /** @brief default constructs n objects at end of
   * data structure, n known at compile time */
  template<size_t n>
  iterator emplace_back_n()
  {
    return emplace_back_n(n);
  }
  /** @brief default constructs n objects at end of
   * data structure, n known at run time
   */
  iterator emplace_back_n(size_t n)
  {
    size_t first = total_size_;
    for (size_t i = 0; i < n; ++i) {
      emplace_back();
    }
    return iterator(this, first, first, total_size_);
  }

  /** @brief number of elements */
  size_t size() const { return total_size_; }

  /** @brief total capacity */
  size_t capacity() const { return data_.size() * buffer_size; }

  /** @brief iterator helpers */
  iterator begin() { return iterator(this, 0); }
  iterator end() { return iterator(this, total_size_); }
  const_iterator begin() const { return const_iterator(this, 0); }
  const_iterator end() const { return const_iterator(this, total_size_); }

  iterator last_checkpoint()
  {
    return iterator(this, checkpoints_.back(), 0, total_size_);
  }
  iterator first_checkpoint()
  {
    return iterator(this, checkpoints_[0], 0, total_size_);
  };
  iterator checkpoint_at(size_t index)
  {
    return iterator(this, checkpoints_[index], 0, total_size_);
  };

  /** @brief searches for item in allocator
   *  only used to find gradient nodes for propagation */
  iterator find(const RealType* const item)
  {
    return std::find_if(
      begin(), end(), [&](const RealType& val) { return &val == item; });
  }
  /** @brief vector like access,
   *  currently unused anywhere but very useful for debugging
   */
  RealType& operator[](std::size_t i)
  {
    BOOST_MATH_ASSERT(i < total_size_);
    return (*data_[i / buffer_size])[i % buffer_size];
  }
  const RealType& operator[](std::size_t i) const
  {
    BOOST_MATH_ASSERT(i < total_size_);
    return (*data_[i / buffer_size])[i % buffer_size];
  }
};
} // namespace detail
} // namespace reverse_mode
} // namespace differentiation
BOOST_MATH_NAMESPACE_END

#endif
