//  (C) Copyright Nick Thompson 2020.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#include <iostream>
#include <boost/math/tools/agm.hpp>
#include <boost/math/constants/constants.hpp>

// This example computes the lemniscate constant to high precision using the agm:
using boost::math::tools::agm;
using boost::math::constants::pi;

int main() {
    using Real = long double;
    Real G = agm(sqrt(Real(2)), Real(1));
    std::cout << std::setprecision(std::numeric_limits<Real>::max_digits10);
    std::cout << "Gauss's lemniscate constant = " << pi<Real>()/G << "\n";
}