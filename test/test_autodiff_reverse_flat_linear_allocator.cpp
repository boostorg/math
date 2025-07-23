#include "test_autodiff_reverse.hpp"
#include <vector>
BOOST_AUTO_TEST_SUITE(test_flat_linear_allocator)

BOOST_AUTO_TEST_CASE_TEMPLATE(flat_linear_allocator_constructors, T, all_float_types)
{
    size_t                              buffer_size = 16;
    RandomSample<T>                     rng{-1, 1};
    rdiff::flat_linear_allocator<T, 16> float_allocator{};
    for (size_t i = 0; i < 2 * buffer_size - buffer_size / 2; i++) {
        float_allocator.emplace_back(rng.next());
    }

    BOOST_CHECK_EQUAL(float_allocator.size(), 2 * buffer_size - buffer_size / 2);
    BOOST_CHECK_EQUAL(float_allocator.capacity(), 2 * buffer_size);

    float_allocator.clear();
    BOOST_CHECK_EQUAL(float_allocator.size(), 0);
    BOOST_CHECK_EQUAL(float_allocator.capacity(), buffer_size);

    for (size_t i = 0; i < 2 * buffer_size - buffer_size / 2; i++) {
        float_allocator.emplace_back(rng.next());
    }
    float_allocator.reset();
    BOOST_CHECK_EQUAL(float_allocator.size(), 0);
    BOOST_CHECK_EQUAL(float_allocator.capacity(), 2 * buffer_size);

    for (size_t i = 0; i < 2 * buffer_size - buffer_size / 2; i++) {
        float_allocator.emplace_back(rng.next());
    }
    T fill_value = T(0.25);
    float_allocator.fill(fill_value);
    for (size_t i = 0; i < float_allocator.size(); i++) {
        BOOST_CHECK_EQUAL(float_allocator[i], fill_value);
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(flat_linear_allocator_test_emplace, T, all_float_types)
{
    size_t                              buffer_size = 16;
    RandomSample<T>                     rng{-1, 1};
    rdiff::flat_linear_allocator<T, 16> float_allocator{};
    std::vector<T>                      test_vector;

    for (int i = 0; i < 2 * buffer_size - 1; i++) {
        test_vector.push_back(rng.next());
        float_allocator.emplace_back(test_vector[i]);
    }

    auto it1 = float_allocator.template emplace_back_n<4>();
    for (int i = 0; i < 4; i++) {
        T literal = rng.next();
        test_vector.push_back(literal);
        *(it1 + i) = literal;
    }

    auto it2 = float_allocator.emplace_back_n(buffer_size);
    for (int i = 0; i < buffer_size; i++) {
        T literal = rng.next();
        test_vector.push_back(literal);
        *(it2 + i) = literal;
    }

    auto vit      = test_vector.begin();
    auto alloc_it = float_allocator.begin();
    for (; vit != test_vector.end(); vit++, alloc_it++) {
        BOOST_CHECK_EQUAL(
            *vit,
            *alloc_it); // should be ok to check floats like this since they are expected to be the same.
    }

    for (int i = 0; i < test_vector.size(); i++) {
        BOOST_CHECK_EQUAL(test_vector[i], float_allocator[i]); // check random access aswell;
    }

    BOOST_CHECK_EQUAL(test_vector.size(), float_allocator.size());
}
BOOST_AUTO_TEST_CASE_TEMPLATE(flat_linear_allocator_test_checkpointing, T, all_float_types)
{
    size_t                              buffer_size = 16;
    RandomSample<T>                     rng{-1, 1};
    rdiff::flat_linear_allocator<T, 16> float_allocator{};
    std::vector<T>                      test_vector;
    std::vector<size_t>                 checkpoint_indices{2, 11, 15, 16, 17, 28};
    std::vector<T>                      expected_value_at_checkpoint;

    size_t                              ckp_id = 0;
    for (int i = 0; i < 2 * buffer_size; i++) {
        T literal = rng.next();
        float_allocator.emplace_back(literal);
        if (i == checkpoint_indices[ckp_id]) {
            float_allocator.add_checkpoint();
            expected_value_at_checkpoint.push_back(literal);
            ++ckp_id;
        }
    }
    for (int i = 0; i < checkpoint_indices.size(); i++) {
        auto it = float_allocator.checkpoint_at(i);
        BOOST_CHECK_EQUAL(*it, expected_value_at_checkpoint[i]);
    }

    auto first_ckp = float_allocator.first_checkpoint();
    auto last_ckp  = float_allocator.last_checkpoint();

    BOOST_CHECK_EQUAL(*first_ckp, expected_value_at_checkpoint[0]);
    BOOST_CHECK_EQUAL(*last_ckp, expected_value_at_checkpoint.back());

    float_allocator.rewind_to_last_checkpoint();
    BOOST_CHECK_EQUAL(float_allocator.size(), checkpoint_indices.back());
    BOOST_CHECK_EQUAL(float_allocator.capacity(), 2 * buffer_size);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(vars_n_stuff, T, all_float_types)
{
    reverse_mode::gradient_tape<T>    tape{};        // depth = 0 for rvar<T>
    reverse_mode::scoped_tape_context context(tape); // Activate tape for current thread

    reverse_mode::rvar<T> x(3.0), y(2.0);
    auto                  z = x * x * x;
    std::cout << decltype(z)::num_literals << std::endl;
    // std::cout << z.evaluate() << std::endl;
    reverse_mode::rvar<T> a(z);
    std::cout << a.evaluate() << std::endl;
    a.print_der_info();
    a.backward();
    std::cout << x.adjoint() << std::endl;
}
BOOST_AUTO_TEST_SUITE_END()
