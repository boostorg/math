#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include <stdexcept>
#include <boost/math/quadrature/gauss.hpp>

using namespace boost::math::quadrature;

// Recursive engine tracking runtime depth and coordinate state
template <typename F, typename Integrator>
double integrate_recursive(F& f, 
                           const std::vector<double>& a, 
                           const std::vector<double>& b, 
                           std::vector<double>& working_coords, 
                           size_t depth,
                           Integrator integrate) {
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

// User-facing entry point supporting any size N
template <typename F>
double integrate_nd(F& f, const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size() || a.empty()) {
        throw std::invalid_argument("Limits vectors must be of equal size and non-empty.");
    }

    typedef std::function<double(double)> func;
    typedef std::function<double(func, double, double)> f2;

    auto integrate = [](func f, double a, double b)
    {
        return boost::math::quadrature::gauss<double, 20>::integrate(f, a, b);
    };
    
    // Allocate the coordinate state tracking vector dynamically based on runtime size
    std::vector<double> working_coords(a.size(), 0.0);
    return integrate_recursive(f, a, b, working_coords, 0, integrate);
}

int main() {
    // Example: 5-Dimensional Integrand f(x₀, x₁, x₂, x₃, x₄) = x₀ * x₁ * x₂ * x₃ * x₄
    auto integrand_5d = [](const std::vector<double>& coords) {
        return sin(coords[0] * coords[1]);
    };

    // Define 5D hypercube boundaries
    std::vector<double> lower_bounds = {0.0, 0.0};
    std::vector<double> upper_bounds = {1.0, 0.5}; // [0,1]^5

    // Exact result: (0.5)^5 = 0.03125
    double result = integrate_nd(integrand_5d, lower_bounds, upper_bounds);
    std::cout << "5D Integrated Result: " << result << std::endl;

    return 0;
}