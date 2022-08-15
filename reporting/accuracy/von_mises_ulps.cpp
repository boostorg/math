//    Copyright Philipp C. J. Muenster, 2020.
//    Copyright Matt Borland, 2022.

//    Use, modification and distribution are subject to the
//    Boost Software License, Version 1.0. (See accompanying file
//    LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <string>

#include <boost/math/distributions/von_mises.hpp>
#include <boost/math/tools/ulps_plot.hpp>

using precise_real = long double;

template <typename CoarseReal>
void generate_ulps_plot_pdf(CoarseReal concentration)
{
    auto hi_conc = static_cast<precise_real>(concentration);
    auto high_precision = [hi_conc](precise_real x)
    {
        boost::math::von_mises_distribution<precise_real> dist(0, hi_conc);
        return boost::math::pdf(dist, x);
    };

    auto low_precision = [concentration](CoarseReal x)
    {
        boost::math::von_mises_distribution<CoarseReal> dist(0, concentration);
        return boost::math::pdf(dist, x);
    };

    using hp_func = decltype (high_precision);

    boost::math::tools::ulps_plot<hp_func, precise_real, CoarseReal> plot(
                high_precision, 0, boost::math::constants::pi<CoarseReal>());
    plot.add_fn(low_precision);
    std::string filename = "von_mises_ulps_pdf_"
                           + std::to_string(static_cast<int>(concentration))
                           + typeid(CoarseReal).name()
                           + ".svg";
    plot.write(filename);
}

template <typename CoarseReal>
void generate_ulps_plot_cdf(CoarseReal concentration)
{
    auto hi_conc = static_cast<precise_real>(concentration);
    auto high_precision = [hi_conc](precise_real x)
    {
        boost::math::von_mises_distribution<precise_real> dist(0, hi_conc);
        return boost::math::cdf(dist, x);
    };

    auto low_precision = [concentration](CoarseReal x)
    {
        boost::math::von_mises_distribution<CoarseReal> dist(0, concentration);
        return boost::math::cdf(dist, x);
    };

    using hp_func = decltype (high_precision);

    auto pi = boost::math::constants::pi<CoarseReal>();
    boost::math::tools::ulps_plot<hp_func, precise_real, CoarseReal> plot(
                high_precision, -pi, +pi);
    plot.add_fn(low_precision);
    std::string filename = "von_mises_ulps_cdf_"
                           + std::to_string(static_cast<int>(concentration))
                           + typeid(CoarseReal).name()
                           + ".svg";
    plot.write(filename);
}

void generate_ulps_plots(double concentration)
{
    auto fconc = static_cast<float>(concentration);
    generate_ulps_plot_pdf<float>(fconc);
    generate_ulps_plot_cdf<float>(fconc);

    generate_ulps_plot_pdf<double>(concentration);
    generate_ulps_plot_cdf<double>(concentration);
}

int main()
{
    generate_ulps_plots(4);
    generate_ulps_plots(64);

    return 0;
}