//  (C) Copyright Nick Thompson 2020.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#include <boost/multiprecision/mpfr.hpp>
#include <boost/math/tools/luroth.hpp>
#include <boost/math/constants/constants.hpp>

using boost::math::constants::pi;
using boost::math::tools::luroth_expansion;
using boost::multiprecision::mpfr_float;

int main() {
    using Real = mpfr_float;
    mpfr_float::default_precision(1024);
    auto lur = luroth(pi<Real>());
    std::cout << lur << "\n";
}
