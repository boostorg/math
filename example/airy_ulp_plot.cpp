#include <iostream>
#include <boost/math/tools/ulp_plot.hpp>
#include <boost/core/demangle.hpp>
#include <boost/multiprecision/float128.hpp>
#include <boost/math/special_functions/airy.hpp>

int main() {
    using PreciseReal = long double;
    using CoarseReal = float;

    auto ai_coarse = boost::math::airy_ai<CoarseReal>;
    auto ai_precise = boost::math::airy_ai<PreciseReal>;

    std::string title = "Airy Ai ULP plot at " + boost::core::demangle(typeid(CoarseReal).name()) + " precision";
    std::string filename = "airy_ai_" + boost::core::demangle(typeid(CoarseReal).name()) + ".svg";
    int samples = 10000;
    int width = 700;
    int clip = 20;
    auto plot = boost::math::tools::ulp_plot<decltype(ai_precise), PreciseReal, CoarseReal>(ai_precise, CoarseReal(-3), CoarseReal(3), true, samples);
    plot.set_clip(clip);
    plot.set_width(width);
    plot.set_title(title);
    plot.set_vertical_lines(6);
    plot.add_fn(ai_coarse);
    std::cout << plot << "\n";
    plot.write(filename);
}