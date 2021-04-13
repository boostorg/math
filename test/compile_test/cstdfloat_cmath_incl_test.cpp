//  Copyright Matt Borland 2021.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Basic sanity check that header
// #includes all the files that it needs to.
//
#include <boost/math/cstdfloat/cstdfloat_cmath.hpp>
//
// Note this header includes no other headers, this is
// important if this test is to be meaningful:
//
#include "test_compile_result.hpp"

void compile_and_link_test()
{
    #ifdef __float128
    __float128 f128 = 0;

    check_result<__float128>(ldexp(f128, 0));
    check_result<__float128>(frexp(f128, 0));
    check_result<__float128>(fabs(f128));
    check_result<__float128>(abs(f128));
    check_result<__float128>(floor(f128));
    check_result<__float128>(ceil(f128));
    check_result<__float128>(sqrt(f128));
    check_result<__float128>(trunc(f128));
    check_result<__float128>(exp(f128));
    check_result<__float128>(expm1(f128));
    check_result<__float128>(pow(f128, 0));
    check_result<__float128>(log(f128));
    check_result<__float128>(log10(f128));
    check_result<__float128>(sin(f128));
    check_result<__float128>(cos(f128));
    check_result<__float128>(tan(f128));
    check_result<__float128>(asin(f128));
    check_result<__float128>(acos(f128));
    check_result<__float128>(atan(f128));
    check_result<__float128>(sinh(f128));
    check_result<__float128>(cosh(f128));
    check_result<__float128>(tanh(f128));
    check_result<__float128>(asinh(f128));
    check_result<__float128>(acosh(f128));
    check_result<__float128>(atanh(f128));
    check_result<__float128>(fmod(f128, f128));
    check_result<__float128>(atan2(f128, f128));
    check_result<__float128>(lgamma(f128));
    check_result<__float128>(tgamma(f128));
    check_result<__float128>(remainder(f128, f128));
    check_result<__float128>(remquo(f128, f128, 0));
    check_result<__float128>(fma(f128, f128, f128));
    check_result<__float128>(fmax(f128, f128));
    check_result<__float128>(fmin(f128, f128));
    check_result<__float128>(fdim(f128, f128));
    check_result<__float128>(nanq(0));
    check_result<__float128>(exp2(f128));
    check_result<__float128>(log2(f128));
    check_result<__float128>(log1p(f128));
    check_result<__float128>(cbrt(f128));
    check_result<__float128>(hypot(f128, f128));
    check_result<__float128>(erf(f128));
    check_result<__float128>(erfc(f128));
    check_result<__float128>(llround(f128));
    check_result<__float128>(lround(f128));
    check_result<__float128>(round(f128));
    check_result<__float128>(nearbyint(f128));
    check_result<__float128>(llrint(f128));
    check_result<__float128>(lrint(f128));
    check_result<__float128>(rint(f128));
    check_result<__float128>(modf(f128, nullptr));
    check_result<__float128>(scalbln(f128, 0));
    check_result<__float128>(scalbn(f128, 0));
    check_result<__float128>(ilogb(f128));
    check_result<__float128>(logb(f128));
    check_result<__float128>(nextafter(f128, f128));
    check_result<__float128>(nexttoward(f128, f128));
    check_result<__float128>(copysign(f128, f128));
    check_result<__float128>(std::signbit(f128));
    check_result<__float128>(std::fpclassify(f128));
    check_result<__float128>(std::isfinite(f128));
    check_result<__float128>(std::isinf(f128));
    check_result<__float128>(std::isnan(f128));
    check_result<__float128>(std::isnormal(f128));
    check_result<__float128>(std::isgreater(f128, f128));
    check_result<__float128>(std::isgreaterequal(f128, f128));
    check_result<__float128>(std::isless(f128, f128));
    check_result<__float128>(std::islessequal(f128, f128));
    check_result<__float128>(std::islessgreater(f128, f128));
    check_result<__float128>(std::isunordered(f128, f128));
    #endif // __float128
}
