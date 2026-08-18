#ifndef BOOST_MATH_QUADRATURE_MULTI_INTEGRALS_HPP
#define BOOST_MATH_QUADRATURE_MULTI_INTEGRALS_HPP

#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include <stdexcept>
#include <boost/math/quadrature/gauss.hpp>

namespace boost { namespace math { namespace quadrature { namespace detail {

// Recursive engine tracking runtime depth and coordinate state
template <typename F, typename Integrator>
double integrate_recursive(const F& f, 
                           const std::vector<double>& a, 
                           const std::vector<double>& b, 
                           std::vector<double>& working_coords, 
                           const size_t& depth,
                           const Integrator& integrate) {
    double low = a[depth];
    double high = b[depth];
    
    // Base Case: Innermost dimension reached
    if (depth == a.size() - 1) {
        // Create a 1D lambda that satisfies Boost's double(double) requirement
        auto innermost_f = [&](double x) {
            working_coords[depth] = x; // Assign final coordinate
            return f(working_coords);  // Evaluate the arbitrary-D function
        };
        return integrate(innermost_f, low, high);
    } 
    
    // Recursive Step: Outer dimensions
    auto next_dimension_f = [&](double x) {
        working_coords[depth] = x; // Assign current dimension's coordinate
        return integrate_recursive(f, a, b, working_coords, depth + 1, integrate);
    };
    
    return integrate(next_dimension_f, low, high);
}

} // namespace detail

template <typename F, typename Integrator>
double integrateND(const F& f, const std::vector<double>& a, const std::vector<double>& b, const Integrator& integrate)
{
    if (a.size() != b.size() || a.empty()) {
        throw std::invalid_argument("Limits vectors must be of equal size and non-empty.");
    }
    
    // Allocate the coordinate state tracking vector dynamically based on runtime size
    std::vector<double> working_coords(a.size(), 0.0);
    return detail::integrate_recursive(f, a, b, working_coords, 0, integrate);
}

template <typename F>
double integrateND(const F& f, const std::vector<double>& a, const std::vector<double>& b)
{
    typedef std::function<double(double)> func;
    typedef std::function<double(func, double, double)> f2;

    auto integrate = [](func f, double a, double b)
    {
        return gauss<double, 20>::integrate(f, a, b);
    };
    
    // Allocate the coordinate state tracking vector dynamically based on runtime size
    std::vector<double> working_coords(a.size(), 0);
    return integrateND(f, a, b, integrate);
}

}}}

#endif // BOOST_MATH_QUADRATURE_MULTI_INTEGRALS_HPP