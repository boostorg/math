//           Copyright Maksym Zhelyeznyakov 2025-2026.
// Distributed under the Boost Software License, Version 1.0.
//      (See accompanying file LICENSE_1_0.txt or copy at
//           https://www.boost.org/LICENSE_1_0.txt)
#include "test_autodiff_reverse.hpp"
#include <vector>

template <class T>
void flat_linear_allocator_constructors()
{
  size_t buffer_size = 16;
  RandomSample<T> rng{ -1, 1 };
  rdiff_detail::flat_linear_allocator<T, 16> float_allocator{};
  for (size_t i = 0; i < 2 * buffer_size - buffer_size / 2; i++) {
    float_allocator.emplace_back(rng.next());
  }

  CHECK_EQUAL(float_allocator.size(),
                    static_cast<size_t>(2 * buffer_size - buffer_size / 2));
  CHECK_EQUAL(float_allocator.capacity(),
                    static_cast<size_t>(2 * buffer_size));

  float_allocator.clear();
  CHECK_EQUAL(float_allocator.size(), static_cast<size_t>(0));
  CHECK_EQUAL(float_allocator.capacity(), buffer_size);

  for (size_t i = 0; i < 2 * buffer_size - buffer_size / 2; i++) {
    float_allocator.emplace_back(rng.next());
  }
  float_allocator.reset();
  CHECK_EQUAL(float_allocator.size(), static_cast<size_t>(0));
  CHECK_EQUAL(float_allocator.capacity(), 2 * buffer_size);

  for (size_t i = 0; i < 2 * buffer_size - buffer_size / 2; i++) {
    float_allocator.emplace_back(rng.next());
  }
  T fill_value = T(0.25);
  float_allocator.fill(fill_value);
  for (size_t i = 0; i < float_allocator.size(); i++) {
    CHECK_EQUAL(float_allocator[i], fill_value);
  }
}

template <class T>
void flat_linear_allocator_test_emplace()
{
  size_t buffer_size = 16;
  RandomSample<T> rng{ -1, 1 };
  rdiff_detail::flat_linear_allocator<T, 16> float_allocator{};
  std::vector<T> test_vector;

  for (size_t i = 0; i < 2 * buffer_size - 1; i++) {
    test_vector.push_back(rng.next());
    float_allocator.emplace_back(test_vector[i]);
  }

  auto it1 = float_allocator.template emplace_back_n<4>();
  for (size_t i = 0; i < 4; i++) {
    T literal = rng.next();
    test_vector.push_back(literal);
    *(it1 + i) = literal;
  }

  auto it2 = float_allocator.emplace_back_n(buffer_size);
  for (size_t i = 0; i < buffer_size; i++) {
    T literal = rng.next();
    test_vector.push_back(literal);
    *(it2 + i) = literal;
  }

  auto vit = test_vector.begin();
  auto alloc_it = float_allocator.begin();
  for (; vit != test_vector.end(); vit++, alloc_it++) {
    CHECK_EQUAL(*vit,
                      *alloc_it); // should be ok to check floats like this
                                  // since they are expected to be the same.
  }

  for (size_t i = 0; i < test_vector.size(); i++) {
    CHECK_EQUAL(test_vector[i],
                      float_allocator[i]); // check random access aswell;
  }

  CHECK_EQUAL(test_vector.size(), float_allocator.size());
}

template <class T>
void flat_linear_allocator_test_checkpointing()
{
  constexpr size_t buffer_size = 16;
  RandomSample<T> rng{ -1, 1 };
  rdiff_detail::flat_linear_allocator<T, buffer_size> float_allocator{};
  std::vector<size_t> checkpoint_indices{ 2, 11, 15, 16, 17, 28 };
  std::vector<T> expected_value_at_checkpoint;

  size_t ckp_id = 0;
  for (size_t i = 0; i < 2 * buffer_size; ++i) {
    T literal = rng.next();
    float_allocator.emplace_back(literal);

    if (ckp_id < checkpoint_indices.size() &&
        (i + 1) == checkpoint_indices[ckp_id]) {
      float_allocator.add_checkpoint();
      expected_value_at_checkpoint.push_back(literal);
      ++ckp_id;
    }
  }

  for (size_t i = 0; i < checkpoint_indices.size(); ++i) {
    auto it = float_allocator.checkpoint_at(i);
    CHECK_TRUE(it != float_allocator.begin());
    --it;
    CHECK_EQUAL(*it, expected_value_at_checkpoint[i]);
  }

  auto first_ckp = float_allocator.first_checkpoint();
  CHECK_TRUE(first_ckp != float_allocator.begin());
  --first_ckp;
  CHECK_EQUAL(*first_ckp, expected_value_at_checkpoint.front());

  auto last_ckp = float_allocator.last_checkpoint();
  CHECK_TRUE(last_ckp != float_allocator.begin());
  --last_ckp;
  CHECK_EQUAL(*last_ckp, expected_value_at_checkpoint.back());

  float_allocator.rewind_to_last_checkpoint();
  CHECK_EQUAL(float_allocator.size(), checkpoint_indices.back());
  CHECK_EQUAL(float_allocator.capacity(), 2 * buffer_size);
}

int main()
{
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { flat_linear_allocator_constructors<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { flat_linear_allocator_test_emplace<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { flat_linear_allocator_test_checkpointing<typename decltype(t)::type>(); });
   return boost::math::test::report_errors();
}
