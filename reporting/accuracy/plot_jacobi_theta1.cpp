
#include <iostream>
#include <boost/math/tools/ulps_plot.hpp>
#include <boost/core/demangle.hpp>
#include <boost/math/special_functions/jacobi_theta.hpp>

using boost::math::tools::ulps_plot;

int main() {
    using PreciseReal = long double; // boost::multiprecision::float128;
    using CoarseReal = double;

    auto jacobi_theta1_coarse = [](CoarseReal z) {
        return boost::math::jacobi_theta1<CoarseReal>(z, 0.5);
    };
    auto jacobi_theta1_precise = [](PreciseReal z) {
        return boost::math::jacobi_theta1<PreciseReal>(z, 0.5);
    };

    std::string filename = "jacobi_theta1_" + boost::core::demangle(typeid(CoarseReal).name()) + ".svg";
    int samples = 2500;
    int width = 1100;
    PreciseReal clip = 100;
    auto plot = ulps_plot<decltype(jacobi_theta1_precise), PreciseReal, CoarseReal>(jacobi_theta1_precise, -10.0, 10.0, samples);
    plot.clip(clip).width(width);
    std::string title = "jacobi_theta1 ULP plot at " + boost::core::demangle(typeid(CoarseReal).name()) + " precision";
    //plot.title(title);
    plot.vertical_lines(10);
    plot.add_fn(jacobi_theta1_coarse);
    plot.write(filename);
}
