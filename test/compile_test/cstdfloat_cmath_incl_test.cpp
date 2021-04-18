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
    #ifdef BOOST_FLOAT128_C
    boost::float128_t f128 = 0;

    check_result<boost::float128_t>(ldexp(f128, 0));
    check_result<boost::float128_t>(frexp(f128, 0));
    check_result<boost::float128_t>(fabs(f128));
    check_result<boost::float128_t>(abs(f128));
    check_result<boost::float128_t>(floor(f128));
    check_result<boost::float128_t>(ceil(f128));
    check_result<boost::float128_t>(sqrt(f128));
    check_result<boost::float128_t>(trunc(f128));
    check_result<boost::float128_t>(exp(f128));
    check_result<boost::float128_t>(expm1(f128));
    check_result<boost::float128_t>(pow(f128, 0));
    check_result<boost::float128_t>(log(f128));
    check_result<boost::float128_t>(log10(f128));
    check_result<boost::float128_t>(sin(f128));
    check_result<boost::float128_t>(cos(f128));
    check_result<boost::float128_t>(tan(f128));
    check_result<boost::float128_t>(asin(f128));
    check_result<boost::float128_t>(acos(f128));
    check_result<boost::float128_t>(atan(f128));
    check_result<boost::float128_t>(sinh(f128));
    check_result<boost::float128_t>(cosh(f128));
    check_result<boost::float128_t>(tanh(f128));
    check_result<boost::float128_t>(asinh(f128));
    check_result<boost::float128_t>(acosh(f128));
    check_result<boost::float128_t>(atanh(f128));
    check_result<boost::float128_t>(fmod(f128, f128));
    check_result<boost::float128_t>(atan2(f128, f128));
    check_result<boost::float128_t>(lgamma(f128));
    check_result<boost::float128_t>(tgamma(f128));
    check_result<boost::float128_t>(remainder(f128, f128));
    check_result<boost::float128_t>(remquo(f128, f128, 0));
    check_result<boost::float128_t>(fma(f128, f128, f128));
    check_result<boost::float128_t>(fmax(f128, f128));
    check_result<boost::float128_t>(fmin(f128, f128));
    check_result<boost::float128_t>(fdim(f128, f128));
    check_result<boost::float128_t>(nanq(0));
    check_result<boost::float128_t>(exp2(f128));
    check_result<boost::float128_t>(log2(f128));
    check_result<boost::float128_t>(log1p(f128));
    check_result<boost::float128_t>(cbrt(f128));
    check_result<boost::float128_t>(hypot(f128, f128));
    check_result<boost::float128_t>(erf(f128));
    check_result<boost::float128_t>(erfc(f128));
    check_result<boost::float128_t>(llround(f128));
    check_result<boost::float128_t>(lround(f128));
    check_result<boost::float128_t>(round(f128));
    check_result<boost::float128_t>(nearbyint(f128));
    check_result<boost::float128_t>(llrint(f128));
    check_result<boost::float128_t>(lrint(f128));
    check_result<boost::float128_t>(rint(f128));
    check_result<boost::float128_t>(modf(f128, nullptr));
    check_result<boost::float128_t>(scalbln(f128, 0));
    check_result<boost::float128_t>(scalbn(f128, 0));
    check_result<boost::float128_t>(ilogb(f128));
    check_result<boost::float128_t>(logb(f128));
    check_result<boost::float128_t>(nextafter(f128, f128));
    check_result<boost::float128_t>(nexttoward(f128, f128));
    check_result<boost::float128_t>(copysign(f128, f128));
    check_result<boost::float128_t>(std::signbit(f128));
    check_result<boost::float128_t>(std::fpclassify(f128));
    check_result<boost::float128_t>(std::isfinite(f128));
    check_result<boost::float128_t>(std::isinf(f128));
    check_result<boost::float128_t>(std::isnan(f128));
    check_result<boost::float128_t>(std::isnormal(f128));
    check_result<boost::float128_t>(std::isgreater(f128, f128));
    check_result<boost::float128_t>(std::isgreaterequal(f128, f128));
    check_result<boost::float128_t>(std::isless(f128, f128));
    check_result<boost::float128_t>(std::islessequal(f128, f128));
    check_result<boost::float128_t>(std::islessgreater(f128, f128));
    check_result<boost::float128_t>(std::isunordered(f128, f128));
    #endif // boost::float128_t
}
