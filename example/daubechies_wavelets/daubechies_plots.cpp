/*
 * Copyright Nick Thompson, 2020
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */
#include <iostream>
#include <boost/core/demangle.hpp>
#include <boost/multiprecision/float128.hpp>
#include <boost/math/special_functions/daubechies_scaling.hpp>
#include <quicksvg/graph_fn.hpp>
#include <quicksvg/ulp_plot.hpp>


using boost::multiprecision::float128;
constexpr const int GRAPH_WIDTH = 700;

template<typename Real, int p>
void plot_phi(int grid_refinements = -1)
{
    auto phi = boost::math::daubechies_scaling<Real, p>();
    if (grid_refinements >= 0)
    {
        phi = boost::math::daubechies_scaling<Real, p>(grid_refinements);
    }
    Real a = 0;
    Real b = phi.support().second;
    std::string title = "Daubechies " + std::to_string(p) + " scaling function";
    std::string filename = "daubechies_" + std::to_string(p) + "_scaling.svg";
    int samples = 1024;
    quicksvg::graph_fn daub(a, b, title, filename, samples, GRAPH_WIDTH);
    daub.set_stroke_width(1);
    daub.add_fn(phi);
    daub.write_all();
}

template<typename Real, int p>
void plot_dphi(int grid_refinements = -1)
{
    auto phi = boost::math::daubechies_scaling<Real, p>();
    if (grid_refinements >= 0)
    {
        phi = boost::math::daubechies_scaling<Real, p>(grid_refinements);
    }
    Real a = 0;
    Real b = phi.support().second;
    std::string title = "Daubechies " + std::to_string(p) + " scaling function derivative";
    std::string filename = "daubechies_" + std::to_string(p) + "_scaling_prime.svg";
    int samples = 1024;
    quicksvg::graph_fn daub(a, b, title, filename, samples, GRAPH_WIDTH);
    daub.set_stroke_width(1);
    auto dphi = [phi](Real x)->Real { return phi.prime(x); };
    daub.add_fn(dphi);
    daub.write_all();
}

template<typename Real, int p>
void plot_convergence()
{
    auto phi0 = boost::math::daubechies_scaling<Real, p>(0);
    Real a = 0;
    Real b = phi0.support().second;
    std::string title = "Daubechies " + std::to_string(p) + " scaling at 0 (green), 1 (orange), 2 (red), and 24 (blue) grid refinements";
    std::string filename = "daubechies_" + std::to_string(p) + "_scaling_convergence.svg";

    quicksvg::graph_fn daub(a, b, title, filename, 1024, GRAPH_WIDTH);
    daub.set_stroke_width(1);

    daub.add_fn(phi0, "green");
    auto phi1 = boost::math::daubechies_scaling<Real, p>(1);
    daub.add_fn(phi1, "orange");
    auto phi2 = boost::math::daubechies_scaling<Real, p>(2);
    daub.add_fn(phi2, "red");

    auto phi21 = boost::math::daubechies_scaling<Real, p>(21);
    daub.add_fn(phi21);
    daub.write_all();
}

template<typename Real, int p>
void plot_condition_number()
{
    using std::abs;
    using std::log;
    static_assert(p >= 3, "p = 2 is not differentiable, so condition numbers cannot be effectively evaluated.");
    auto phi = boost::math::daubechies_scaling<Real, p>();
    Real a = std::sqrt(std::numeric_limits<Real>::epsilon());
    Real b = phi.support().second - 1000*std::sqrt(std::numeric_limits<Real>::epsilon());
    std::string title = "log10 of condition number of function evaluation for Daubechies " + std::to_string(p) + " scaling function.";
    std::string filename = "daubechies_" + std::to_string(p) + "_scaling_condition_number.svg";


    quicksvg::graph_fn daub(a, b, title, filename, 2048, GRAPH_WIDTH);
    daub.set_stroke_width(1);

    auto cond = [&phi](Real x)
    {
        Real y = phi(x);
        Real dydx = phi.prime(x);
        Real z = abs(x*dydx/y);
        using std::isnan;
        if (z==0)
        {
            return Real(-1);
        }
        if (isnan(z))
        {
            // Graphing libraries don't like nan's:
            return Real(1);
        }
        return log10(z);
    };
    daub.add_fn(cond);
    daub.write_all();
}

template<typename CoarseReal, typename PreciseReal, int p, class PhiPrecise>
void do_ulp(int coarse_refinements, PhiPrecise phi_precise)
{
    auto phi_coarse = boost::math::daubechies_scaling<CoarseReal, p>(coarse_refinements);

    std::string title = std::to_string(p) + " vanishing moment ULP plot at " + std::to_string(coarse_refinements) + " refinements and " + boost::core::demangle(typeid(CoarseReal).name()) + " precision";
    std::string filename = "daubechies_" + std::to_string(p) + "_" + boost::core::demangle(typeid(CoarseReal).name()) + "_" + std::to_string(coarse_refinements) + "_refinements.svg";
    int samples = 20000;
    int clip = 20;
    quicksvg::ulp_plot<decltype(phi_coarse), CoarseReal, decltype(phi_precise), PreciseReal>(phi_coarse, phi_precise, CoarseReal(0), phi_coarse.support().second, title, filename, samples, GRAPH_WIDTH, clip);
}


int main()
{
    plot_phi<double, 2>();
    plot_phi<double, 3>();
    plot_phi<double, 4>();
    plot_phi<double, 5>();
    plot_phi<double, 6>();
    plot_phi<double, 7>();
    plot_phi<double, 8>();
    plot_phi<double, 9>();
    plot_phi<double, 10>();
    plot_phi<double, 11>();
    plot_phi<double, 12>();
    plot_phi<double, 13>();
    plot_phi<double, 14>();
    plot_phi<double, 15>();

    plot_dphi<double, 3>();
    plot_dphi<double, 4>();
    plot_dphi<double, 5>();
    plot_dphi<double, 6>();
    plot_dphi<double, 7>();
    plot_dphi<double, 8>();
    plot_dphi<double, 9>();
    plot_dphi<double, 10>();
    plot_dphi<double, 11>();
    plot_dphi<double, 12>();
    plot_dphi<double, 13>();
    plot_dphi<double, 14>();
    plot_dphi<double, 15>();

    plot_condition_number<long double, 3>();
    plot_condition_number<long double, 4>();
    plot_condition_number<long double, 5>();
    plot_condition_number<long double, 6>();
    plot_condition_number<long double, 7>();
    plot_condition_number<long double, 8>();
    plot_condition_number<long double, 9>();
    plot_condition_number<long double, 10>();
    plot_condition_number<long double, 11>();
    plot_condition_number<long double, 12>();
    plot_condition_number<long double, 13>();

    plot_convergence<double, 2>();
    plot_convergence<double, 3>();
    plot_convergence<double, 4>();
    plot_convergence<double, 5>();
    plot_convergence<double, 6>();
    plot_convergence<double, 7>();
    plot_convergence<double, 8>();
    plot_convergence<double, 9>();
    plot_convergence<double, 10>();
    plot_convergence<double, 11>();
    plot_convergence<double, 12>();
    plot_convergence<double, 13>();
    plot_convergence<double, 14>();
    plot_convergence<double, 15>();


    using PreciseReal = float128;
    using CoarseReal = double;
    int precise_refinements = 22;
    constexpr const int p = 12;
    std::cout << "Computing precise scaling function in " << boost::core::demangle(typeid(PreciseReal).name()) << " precision.\n";
    auto phi_precise = boost::math::daubechies_scaling<PreciseReal, p>(precise_refinements);
    std::cout << "Beginning comparison with functions computed in " << boost::core::demangle(typeid(CoarseReal).name()) << " precision.\n";
    for (int i = 7; i <= precise_refinements-1; ++i)
    {
        std::cout << "\tCoarse refinement " << i << "\n";
        do_ulp<CoarseReal, PreciseReal, p>(i, phi_precise);
    }
}
