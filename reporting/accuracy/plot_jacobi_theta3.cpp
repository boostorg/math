
#include <iostream>
#include <boost/math/tools/ulps_plot.hpp>
#include <boost/core/demangle.hpp>
#include <boost/math/special_functions/jacobi_theta.hpp>

using boost::math::tools::ulps_plot;

int main() {
    using PreciseReal = long double;
    using CoarseReal = float;

    auto jacobi_theta_coarse = [](CoarseReal z) {
        return boost::math::jacobi_theta3m1<CoarseReal>(z, 0.05);
    };
    auto jacobi_theta_precise = [](PreciseReal z) {
        return boost::math::jacobi_theta3m1<PreciseReal>(z, 0.05);
    };

    std::string filename = "jacobi_theta3_" + boost::core::demangle(typeid(CoarseReal).name()) + ".svg";
    int samples = 2500;
    int width = 1100;
    PreciseReal clip = 100;
    auto plot = ulps_plot<decltype(jacobi_theta_precise), PreciseReal, CoarseReal>(jacobi_theta_precise, CoarseReal(-10.0), CoarseReal(10.0), samples);
    plot.clip(clip).width(width);
    std::string title = "jacobi_theta3 ULP plot at " + boost::core::demangle(typeid(CoarseReal).name()) + " precision";
    plot.title(title);
    plot.vertical_lines(10);
    plot.add_fn(jacobi_theta_coarse);
    plot.write(filename);
}
