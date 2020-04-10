//  (C) Copyright Nick Thompson 2020.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#include <iostream>
#include <boost/math/tools/ulp_plot.hpp>
#include <boost/core/demangle.hpp>
#include <boost/math/special_functions/airy.hpp>

using boost::math::tools::ulp_plot;

int main() {
    using PreciseReal = long double;
    using CoarseReal = float;

    auto ai_coarse = boost::math::airy_ai<CoarseReal>;
    auto ai_precise = boost::math::airy_ai<PreciseReal>;

    std::string filename = "airy_ai_" + boost::core::demangle(typeid(CoarseReal).name()) + ".svg";
    int samples = 10000;
    // How many pixels wide do you want your .svg?
    int width = 700;
    // Near a root, we have unbounded relative error. So for functions with roots, we define a ULP clip:
    int clip = 20;
    // Should we perturb the abscissas? i.e., should we compute the high precision function f at x,
    // and the low precision function at the nearest representable x̂ to x?
    // Or should we compute both the high precision and low precision function at a low precision representable x̂?
    // If we are displaying an ULP envelope, it makes sense to perturb the abscissas.
    bool perturb_abscissas = true;
    auto plot = ulp_plot<decltype(ai_precise), PreciseReal, CoarseReal>(ai_precise, CoarseReal(-3), CoarseReal(3), perturb_abscissas, samples);
    plot.set_clip(clip);
    plot.set_width(width);
    // Sometimes it's useful to set a title, but in many cases it's more useful to just use a caption (then no title required)
    //std::string title = "Airy Ai ULP plot at " + boost::core::demangle(typeid(CoarseReal).name()) + " precision";
    //plot.set_title(title);
    plot.set_vertical_lines(6);
    plot.add_fn(ai_coarse);
    // You can write the plot to a stream:
    std::cout << plot;
    // Or to a file:
    plot.write(filename);

    // Don't like the default dark theme?
    plot.set_background_color("white");
    plot.set_font_color("black");
    filename =  "airy_ai_" + boost::core::demangle(typeid(CoarseReal).name()) + "_white.svg";
    plot.write(filename);

    // Don't like the envelope?
    plot.write_ulp_envelope(false);
    filename =  "airy_ai_" + boost::core::demangle(typeid(CoarseReal).name()) + "_white_no_envelope.svg";
    plot.write(filename);
}
