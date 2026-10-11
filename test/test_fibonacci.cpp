//  Copyright Madhur Chauhan 2020.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/policies/error_handling.hpp>
#include <boost/math/special_functions/fibonacci.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include "math_unit_test.hpp"
#include <cstdint>
#include <exception>

using boost::math::fibonacci;
using boost::math::unchecked_fibonacci;
using namespace boost::multiprecision;
typedef cpp_int BST;
typedef number<backends::cpp_int_backend<2048, 2048, unsigned_magnitude, unchecked>> BST_2048;

// Some sanity checks using OEIS A000045
void sanity_checks() {
    CHECK_TRUE(fibonacci<char>(0) == 0);    // Base case
    CHECK_TRUE(fibonacci<int>(1) == 1);     // Base case
    CHECK_TRUE(fibonacci<int16_t>(2) == 1); // Base Case
    CHECK_TRUE(fibonacci<int16_t>(3) == 2); // First computation
    CHECK_TRUE(fibonacci<int16_t>(10) == 55);
    CHECK_TRUE(fibonacci<int16_t>(15) == 610);
    CHECK_TRUE(fibonacci<int16_t>(18) == 2584);
    CHECK_TRUE(fibonacci<int32_t>(40) == 102334155);
}

// Tests unchecked_fibonacci by computing naively
void big_integer_check() {
    // GMP is used as type for fibonacci and cpp_int is for naive computation
    BST val = 0, a = 0, b = 1;
    for (int i = 0; i <= 1e4; ++i) {
        CHECK_TRUE(fibonacci<BST>(i) == a);
        val = b, b += a, a = val;
    }
}

// Check for overflow throw using magic constants
void overflow_check() {

    // 1. check for unsigned integer overflow
    CHECK_NO_THROW(fibonacci<uint64_t>(93));
    CHECK_THROW(fibonacci<uint64_t>(94), std::exception);
    CHECK_NO_THROW(unchecked_fibonacci<uint64_t>(94));

    // 2. check for signed integer overflow
    CHECK_NO_THROW(fibonacci<int64_t>(91));
    // BOOST_CHECK_THROW(fibonacci<int64_t>(92), std::exception); // this should be the correct value but imprecisions
    CHECK_THROW(fibonacci<int64_t>(93), std::exception);
    // In UBSAN this will error out, but it is expected
    CHECK_NO_THROW(unchecked_fibonacci<int64_t>(93));

    // 3. check for floating point (double)
    CHECK_NO_THROW(fibonacci<double>(78));
    CHECK_THROW(fibonacci<double>(79), std::exception);
    CHECK_NO_THROW(unchecked_fibonacci<double>(79));

    // 4. check using boost's multiprecision unchecked integer
    CHECK_NO_THROW(fibonacci<BST_2048>(2950));
    // BOOST_CHECK_THROW(fibonacci<T>(2951), std::exception); // this should be the correct value but imprecisions
    CHECK_THROW(fibonacci<BST_2048>(2952), std::exception);
    CHECK_NO_THROW(unchecked_fibonacci<BST_2048>(2952));
}

void generator_check() {
    // first 5 values
    boost::math::fibonacci_generator<BST_2048> gen;
    for (int i : {0, 1, 1, 2, 3, 5, 8, 13, 21}) {
        CHECK_TRUE(gen() == i);
    }

    // test whether the generator is set correctly to the given index --- next
    const int next = 1000; // next <=2950 (checked from test above)
    gen.set(next);
    BST_2048 a = fibonacci<BST_2048>(next), b = fibonacci<BST_2048>(next + 1);
    for (int i = next; i < next + 50; ++i) {
        CHECK_TRUE(gen() == a);
        swap(a, b);
        b += a;
    }

    // shift the generator back to next and check again
    a = fibonacci<BST_2048>(next), b = fibonacci<BST_2048>(next + 1);
    gen.set(next);
    for (int i = next; i < next + 50; ++i) {
        CHECK_TRUE(gen() == a);
        swap(a, b);
        b += a;
    }
}

#ifndef BOOST_NO_CXX14_CONSTEXPR

void constexpr_check() {
    constexpr int x = boost::math::unchecked_fibonacci<int>(32);
    static_assert(x == 2178309, "constexpr check 1");

    constexpr double y = boost::math::unchecked_fibonacci<double>(40);
    static_assert(y == 102334155.0, "constexpr check 2");

    constexpr int xx = boost::math::fibonacci<int>(32);
    static_assert(xx == 2178309, "constexpr check 3");

    constexpr double yy = boost::math::fibonacci<double>(40);
    static_assert(yy == 102334155.0, "constexpr check 4");

}

#endif

int main()
{
   sanity_checks();
   big_integer_check();
   overflow_check();
   generator_check();
#ifndef BOOST_NO_CXX14_CONSTEXPR
   constexpr_check();
#endif
   return boost::math::test::report_errors();
}
