// Copyright Jacob Hass, 2026
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_QUADRATURE_MULTI_INTEGRALS_HPP
#define BOOST_MATH_QUADRATURE_MULTI_INTEGRALS_HPP

#include <iostream>
#include <vector>
#include <string>
#include <functional>
#include <cmath>
#include <stdexcept>
#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>
#include <boost/math/policies/policy.hpp>
#include <boost/math/policies/error_handling.hpp>

namespace boost { namespace math { namespace quadrature { namespace detail {

// Recursive engine tracking runtime depth and coordinate state
template <typename F, typename Integrator, typename RealType>
auto integrate_recursive(const F& f, 
                         const std::vector<RealType>& a, 
                         const std::vector<RealType>& b, 
                         std::vector<RealType>& working_coords, 
                         const size_t& depth,
                         const Integrator& integrate) -> decltype(f(a))
{
    RealType low = a[depth];
    RealType high = b[depth];
    
    // Base Case: Innermost dimension reached
    if (depth == a.size() - 1) {
        auto innermost_f = [&](const RealType& x) {
            working_coords[depth] = x; // Assign final coordinate
            return f(working_coords);  // Evaluate the arbitrary-D function
        };
        return integrate(innermost_f, low, high);
    }
    
    // Recursive Step: Outer dimensions
    auto next_dimension_f = [&](const RealType& x) {
        working_coords[depth] = x; // Assign current dimension's coordinate
        return integrate_recursive(f, a, b, working_coords, depth + 1, integrate);
    };
    
    return integrate(next_dimension_f, low, high);
}

} // namespace detail

template <typename F, typename Integrator, typename RealType, class Policy>
auto integrateND(const F& f, const std::vector<RealType>& a, const std::vector<RealType>& b, 
                                   const Integrator& integrate, const Policy& pol) -> decltype(f(a))
{
    static const char* function = "boost::math::quadrature::integrateND(f, a, b, integrate, pol)";

    // Error checking on integration bounds
    if (a.empty() || b.empty()){
        boost::math::policies::raise_evaluation_error(function, "Integration limits cannot be empty.\n", 0, pol);
    }
    if (a.size() != b.size()) {
        std::string msg = "Integration limits must be the same size but got (" + std::to_string(a.size()) + "!=" + std::to_string(b.size()) + ")";
        boost::math::policies::raise_evaluation_error(function, msg.c_str(), 0, pol);
    }
    
    // Allocate the coordinate state tracking vector dynamically based on runtime size
    std::vector<RealType> working_coords(a.size(), 0.0);
    return detail::integrate_recursive(f, a, b, working_coords, 0, integrate);
}

template <typename F, typename Integrator, typename RealType>
auto integrateND(const F& f, const std::vector<RealType>& a, const std::vector<RealType>& b, 
                 const Integrator& integrate) -> decltype(f(a))
{
    return integrateND(f, a, b, integrate, boost::math::policies::policy<>());
}

}}}

#endif // BOOST_MATH_QUADRATURE_MULTI_INTEGRALS_HPP