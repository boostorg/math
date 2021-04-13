//  Copyright Matt Borland 2021.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Basic sanity check that header
// #includes all the files that it needs to.
//
#include <boost/math/cstdfloat/cstdfloat_limits.hpp>
//
// Note this header includes no other headers, this is
// important if this test is to be meaningful:
//
#include "test_compile_result.hpp"

void compile_and_link_test()
{
    #ifdef __float128
    check_result<bool>(numeric_limits<__float128>::is_specialized);
    check_result<__float128>(numeric_limts<__float128>::(min)());
    check_result<__float128>(numeric_limts<__float128>::(max)());
    check_result<__float128>(numeric_limts<__float128>::lowest());
    check_result<int>(numeric_limts<__float128>::digits);
    check_result<int>(numeric_limts<__float128>::digits10);
    check_result<int>(numeric_limts<__float128>::max_digits10);
    check_result<bool>(numeric_limts<__float128>::is_signed);
    check_result<bool>(numeric_limts<__float128>::is_integer);
    check_result<bool>(numeric_limts<__float128>::is_exact);
    check_result<int>(numeric_limts<__float128>::radix);
    check_result<__float128>(numeric_limits<__float128>::epsilon());
    check_result<int>(numeric_limts<__float128>::min_exponent);
    check_result<int>(numeric_limts<__float128>::min_exponent10);
    check_result<int>(numeric_limts<__float128>::max_exponent);
    check_result<int>(numeric_limts<__float128>::max_exponent10);
    check_result<bool>(numeric_limts<__float128>::has_infinity);
    check_result<bool>(numeric_limts<__float128>::has_quiet_NaN);
    check_result<bool>(numeric_limts<__float128>::has_signaling_NaN);
    check_result<float_denorm_style>(numeric_limts<__float128>::has_denorm);
    check_result<bool>(numeric_limts<__float128>::has_denorm_loss);
    check_result<__float128>(numeric_limits<__float128>::infinity());
    check_result<__float128>(numeric_limits<__float128>::quiet_NaN());
    check_result<__float128>(numeric_limits<__float128>::signaling_NaN());
    check_result<__float128>(numeric_limits<__float128>::denorm_min());
    check_result<bool>(numeric_limts<__float128>::is_iec559);
    check_result<bool>(numeric_limts<__float128>::is_bounded);
    check_result<bool>(numeric_limts<__float128>::is_modulo);
    check_result<bool>(numeric_limts<__float128>::traps);
    check_result<bool>(numeric_limts<__float128>::tinyness_before);
    check_result<float_round_style>(numeric_limts<__float128>::round_style);
    #endif
}
