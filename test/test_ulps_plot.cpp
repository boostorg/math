//  (C) Copyright Nick Thompson 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Tests the condition numbers that ulps_plot and ulps_heatmap compute by finite differences on their grids,
// against analytic ones, and the files they write, with and without a PNG encoder.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#include <boost/math/tools/ulps_heatmap.hpp>
#include <boost/math/tools/ulps_plot.hpp>
#include "math_unit_test.hpp"

using boost::math::tools::ulps_heatmap;
using boost::math::tools::ulps_plot;
namespace detail = boost::math::tools::detail;

// A uniform grid of n cell centres on [a, b], and m random probes in [a, b] plus the given extra ones.
struct samples_1d
{
    std::vector<double> x, fx, px, fpx;
};

template<class F>
samples_1d make_samples(F f, double a, double b, std::size_t n, std::size_t m, std::vector<double> const & extra)
{
    samples_1d s;
    for (std::size_t i = 0; i < n; ++i)
    {
        s.x.push_back(a + (b - a)*(i + 0.5)/n);
        s.fx.push_back(f(s.x.back()));
    }
    std::mt19937_64 gen(12345);
    std::uniform_real_distribution<double> dis(a, b);
    for (std::size_t i = 0; i < m; ++i)
    {
        s.px.push_back(dis(gen));
    }
    s.px.insert(s.px.end(), extra.begin(), extra.end());
    for (double p : s.px)
    {
        s.fpx.push_back(f(p));
    }
    return s;
}

// Checks the condition numbers at the grid and the probes against the exact ones, to a relative tolerance,
// and returns the number of local fallbacks.
template<class F, class Exact>
std::size_t check_1d(F f, Exact exact, double a, double b, std::vector<double> const & extra, double tol)
{
    samples_1d s = make_samples(f, a, b, 500, 200, extra);
    std::vector<double> pcond;
    std::size_t fallbacks = 0;
    std::vector<double> cond = detail::grid_condition_numbers(f, s.x, s.fx, s.px, s.fpx, pcond, &fallbacks);
    for (std::size_t i = 0; i < s.x.size(); ++i)
    {
        double e = exact(s.x[i]);
        if (!CHECK_LE(std::abs(cond[i] - e), tol*(std::max)(1.0, e)))
        {
            std::cerr << "  at grid point x = " << s.x[i] << ": " << cond[i] << " vs " << e << "\n";
        }
    }
    for (std::size_t i = 0; i < s.px.size(); ++i)
    {
        double e = exact(s.px[i]);
        if (!CHECK_LE(std::abs(pcond[i] - e), tol*(std::max)(1.0, e)))
        {
            std::cerr << "  at probe x = " << s.px[i] << ": " << pcond[i] << " vs " << e << "\n";
        }
    }
    return fallbacks;
}

void test_condition_numbers_1d()
{
    // exp: |x f'/f| = |x|. Smooth, so the grid is trusted except within 8 grid steps of 0, where the stencils are
    // wide compared to |x| (detail::fd_max_span):
    auto exp_fallbacks = check_1d([](double x) { return std::exp(x); }, [](double x) { return std::abs(x); }, -5.0, 5.0, {}, 1e-6);
    CHECK_LE(exp_fallbacks, std::size_t(30));

    // log near 1: |1/log(x)|, a narrow spike at the zero x = 1. The probes next to it must see the spike, not a smeared version:
    std::vector<double> near_one = {1 + 1e-9, 1 - 1e-9, 1 + 1e-6, 1 - 1e-6, 1 + 1e-3, std::nextafter(1.0, 2.0)};
    auto log_fallbacks = check_1d([](double x) { return std::log(x); }, [](double x) { return std::abs(1/std::log(x)); }, 0.5, 2.0, near_one, 1e-6);
    CHECK_EQUAL(log_fallbacks, std::size_t(0));

    // sin near pi: |x cos(x)/sin(x)|, with a spike at the zero x = pi:
    double const pi = 3.141592653589793;
    std::vector<double> near_pi = {pi, std::nextafter(pi, 4.0), pi + 1e-8, pi - 1e-8, pi + 1e-5};
    auto sin_fallbacks = check_1d([](double x) { return std::sin(x); }, [](double x) { return std::abs(x*std::cos(x)/std::sin(x)); }, 2.0, 4.0, near_pi, 1e-6);
    CHECK_EQUAL(sin_fallbacks, std::size_t(0));

    // tan has a pole at pi/2, across which finite differences on the grid are meaningless. They must be detected,
    // and the condition number |x/(sin(x) cos(x))| computed locally. The envelope is meant to be good to a few percent
    // (detail::fd_tolerance), and next to a pole that is what it gets:
    std::vector<double> near_pole = {pi/2 + 1e-3, pi/2 - 1e-3, pi/2 + 1e-2};
    auto tan_fallbacks = check_1d([](double x) { return std::tan(x); }, [](double x) { return std::abs(x/(std::sin(x)*std::cos(x))); }, 1.0, 2.0, near_pole, 1e-2);
    CHECK_LE(std::size_t(1), tan_fallbacks);
    // Only the few points near the pole need them:
    CHECK_LE(tan_fallbacks, std::size_t(30));

    // log(x) - 1 on [0.5, 32768] with the 8 grid points of a 10-sample logarithmic ulps_plot, a factor of 4 apart.
    // The five-point stencils are as wide as x itself, and agree with the three-point ones although both are badly wrong,
    // so they must not be trusted:
    {
        auto f = [](double x) { return std::log(x) - 1; };
        std::vector<double> x, fx, px = {3.0, 100.0, 5000.0}, fpx, pcond;
        for (int i = 0; i < 8; ++i)
        {
            x.push_back(std::exp(std::log(0.5) + (std::log(32768.0) - std::log(0.5))*(i + 0.5)/8));
            fx.push_back(f(x.back()));
        }
        for (double p : px) { fpx.push_back(f(p)); }
        std::vector<double> cond = detail::grid_condition_numbers(f, x, fx, px, fpx, pcond);
        for (std::size_t i = 0; i < x.size(); ++i)
        {
            CHECK_LE(std::abs(cond[i] - std::abs(1/fx[i])), 1e-6*(std::max)(1.0, std::abs(1/fx[i])));
        }
        for (std::size_t i = 0; i < px.size(); ++i)
        {
            CHECK_LE(std::abs(pcond[i] - std::abs(1/fpx[i])), 1e-6*(std::max)(1.0, std::abs(1/fpx[i])));
        }
    }
}

void test_condition_numbers_2d()
{
    // f(x, y) = sin(x) exp(y): |x cot(x)| + |y|, with a line of zeros at x = pi.
    auto f = [](double x, double y) { return std::sin(x)*std::exp(y); };
    auto exact = [](double x, double y) { return std::abs(x*std::cos(x)/std::sin(x)) + std::abs(y); };
    std::size_t const W = 120;
    std::size_t const H = 80;
    std::vector<double> x(W), y(H), z(W*H);
    for (std::size_t i = 0; i < W; ++i) { x[i] = 2 + 2*(i + 0.5)/W; }
    // Rows from the top, as ulps_heatmap has them:
    for (std::size_t j = 0; j < H; ++j) { y[j] = 1 - 2*(j + 0.5)/H; }
    for (std::size_t j = 0; j < H; ++j)
    {
        for (std::size_t i = 0; i < W; ++i)
        {
            z[j*W + i] = f(x[i], y[j]);
        }
    }
    double const pi = 3.141592653589793;
    std::vector<double> qx = {pi + 1e-9, pi - 1e-6, 2.5, 3.9, 2.01};
    std::vector<double> qy = {0.3, -0.7, 0.0001, 0.99, -0.99};
    std::mt19937_64 gen(5);
    std::uniform_real_distribution<double> dx(2, 4), dy(-1, 1);
    for (int k = 0; k < 500; ++k)
    {
        qx.push_back(dx(gen));
        qy.push_back(dy(gen));
    }
    std::vector<double> qz(qx.size());
    for (std::size_t p = 0; p < qx.size(); ++p) { qz[p] = f(qx[p], qy[p]); }
    std::vector<double> cond, qcond;
    detail::grid_condition_numbers_2d(f, x, y, z, qx, qy, qz, cond, qcond, 2);
    for (std::size_t j = 0; j < H; ++j)
    {
        for (std::size_t i = 0; i < W; ++i)
        {
            double e = exact(x[i], y[j]);
            CHECK_LE(std::abs(cond[j*W + i] - e), 1e-5*(std::max)(1.0, e));
        }
    }
    for (std::size_t p = 0; p < qx.size(); ++p)
    {
        double e = exact(qx[p], qy[p]);
        if (!CHECK_LE(std::abs(qcond[p] - e), 1e-5*(std::max)(1.0, e)))
        {
            std::cerr << "  at probe (" << qx[p] << ", " << qy[p] << "): " << qcond[p] << " vs " << e << "\n";
        }
    }
}

std::string read_file(std::string const & name)
{
    std::ifstream fs(name);
    std::stringstream ss;
    ss << fs.rdbuf();
    return ss.str();
}

std::size_t count(std::string const & s, std::string const & what)
{
    std::size_t n = 0;
    for (std::size_t k = s.find(what); k != std::string::npos; k = s.find(what, k + 1))
    {
        ++n;
    }
    return n;
}

// The parts of an SVG file that must be there for it to be one:
void check_svg(std::string const & s)
{
    CHECK_EQUAL(s.rfind("<?xml", 0), std::size_t(0));
    CHECK_EQUAL(count(s, "<svg "), std::size_t(1));
    CHECK_TRUE(s.size() >= 7 && s.substr(s.size() - 7) == "</svg>\n");
    CHECK_EQUAL(count(s, "<g"), count(s, "</g>"));
    CHECK_EQUAL(count(s, "<mask "), count(s, "</mask>"));
    CHECK_EQUAL(count(s, "<defs>"), count(s, "</defs>"));
}

// Stands in for a PNG encoder: records the images, and returns their pixels as the "PNG".
struct recording_encoder
{
    std::vector<std::size_t> * widths;
    std::vector<std::size_t> * heights;
    bool * gray;
    std::vector<unsigned char> operator()(std::vector<unsigned char> const & rgb, std::size_t w, std::size_t h) const
    {
        widths->push_back(w);
        heights->push_back(h);
        for (std::size_t k = 0; k + 2 < rgb.size(); k += 3)
        {
            *gray = *gray && rgb[k] == rgb[k + 1] && rgb[k] == rgb[k + 2];
        }
        return rgb;
    }
};

// The masks of a raster canvas, and the order of their colors:
struct masks
{
    std::vector<std::vector<unsigned char>> images;
    std::vector<std::string> colors;
};

masks finish(detail::raster_canvas const & canvas)
{
    masks m;
    std::ostringstream os;
    canvas.finish(os, [&](std::vector<unsigned char> const & rgb, std::size_t, std::size_t) { m.images.push_back(rgb); return rgb; }, "t-");
    std::string s = os.str();
    for (std::size_t k = s.find("<rect "); k != std::string::npos; k = s.find("<rect ", k + 1))
    {
        std::size_t a = s.find("fill='", k) + 6;
        m.colors.push_back(s.substr(a, s.find('\'', a) - a));
    }
    return m;
}

void test_raster_canvas()
{
    // Overlapping shapes of one color cover their union: two unit disks centred on the corners of pixel (0, 0) each cover
    // 13 of its 16 subsamples, and together all of them:
    {
        detail::raster_canvas canvas(0, 0, 2, 2);
        canvas.disk("red", 0, 0, 1);
        canvas.disk("red", 1, 0, 1);
        masks m = finish(canvas);
        CHECK_EQUAL(m.images.size(), std::size_t(1));
        CHECK_EQUAL(int(m.images[0][0]), 255);
    }
    // A polyline's joints aren't painted twice: two collinear segments make the same mask as one:
    {
        detail::raster_canvas one(0, 0, 20, 10);
        detail::raster_canvas two(0, 0, 20, 10);
        one.segment("red", 1.3, 2.2, 18.7, 7.9, 1);
        two.polyline("red", {{1.3, 2.2}, {10.0, 5.05}, {18.7, 7.9}}, 1);
        CHECK_TRUE(finish(one).images == finish(two).images);
    }
    // Shapes are painted in the order they are drawn: a red disk drawn over a blue one is on top of it...
    {
        detail::raster_canvas canvas(0, 0, 8, 8);
        canvas.disk("red", 2, 2, 1.5);
        canvas.disk("blue", 4, 4, 2);
        canvas.disk("red", 4, 4, 1);
        masks m = finish(canvas);
        CHECK_TRUE(m.colors == std::vector<std::string>({"red", "blue", "red"}));
    }
    // ...but one that doesn't touch the blue one shares the first red mask:
    {
        detail::raster_canvas canvas(0, 0, 16, 8);
        canvas.disk("red", 2, 2, 1.5);
        canvas.disk("blue", 8, 4, 2);
        canvas.disk("red", 14, 4, 1);
        masks m = finish(canvas);
        CHECK_TRUE(m.colors == std::vector<std::string>({"red", "blue"}));
    }
}

std::string first_mask_id(std::string const & svg)
{
    std::size_t a = svg.find("<mask id='") + 10;
    return svg.substr(a, svg.find('\'', a) - a);
}

void test_files()
{
    // lodepng is the default encoder only on request, so that this program, which isn't linked with lodepng, links even if
    // lodepng.h is on the include path:
#if !defined(BOOST_MATH_USE_LODEPNG) && !defined(LODEPNG_H)
    CHECK_TRUE(!detail::default_png_encoder());
#endif

    CHECK_EQUAL(detail::base64({'M', 'a', 'n'}), std::string("TWFu"));
    CHECK_EQUAL(detail::base64({'M', 'a'}), std::string("TWE="));
    CHECK_EQUAL(detail::base64({'M'}), std::string("TQ=="));

    auto f = [](double x) { return std::sin(x); };
    auto g = [](float x) { return std::sin(x); };
    auto h = [](float x) { return static_cast<float>(std::sin(static_cast<double>(x))); };
    ulps_plot<decltype(f), double, float> plot(f, 0.5f, 6.0f, 500, false, 3);
    plot.add_fn(g).add_fn(h, "orange").clip(5).width(500).show_function(true);

    // Without an encoder, everything is vector, as before: a circle for each point within the clip, and none of them for
    // the cropped points, which are crosses.
    plot.png_encoder(nullptr);
    std::string const vector_name = "test_ulps_plot_vector.svg";
    plot.write(vector_name);
    std::string v = read_file(vector_name);
    check_svg(v);
    CHECK_EQUAL(count(v, "<image"), std::size_t(0));
    CHECK_LE(std::size_t(500), count(v, "<circle"));
    std::remove(vector_name.c_str());

    // With one, there is a mask for each color in the data layer (two functions and the envelope), and one for the
    // function's curve, each a grayscale image of the plot area:
    std::vector<std::size_t> widths, heights;
    bool gray = true;
    plot.png_encoder(recording_encoder{&widths, &heights, &gray});
    std::string const raster_name = "test_ulps_plot_raster.svg";
    plot.write(raster_name);
    std::string r = read_file(raster_name);
    check_svg(r);
    CHECK_EQUAL(count(r, "<circle"), std::size_t(0));
    CHECK_EQUAL(count(r, "<mask "), std::size_t(4));
    CHECK_EQUAL(count(r, "href='data:image/png;base64,"), std::size_t(4));
    CHECK_EQUAL(widths.size(), std::size_t(4));
    CHECK_TRUE(gray);
    for (std::size_t w : widths)
    {
        CHECK_EQUAL(w, widths[0]);
    }
    CHECK_LE(std::size_t(400), widths[0]);
    std::remove(raster_name.c_str());

    // Plots that differ only in their settings get different mask ids, so they can be inlined in one page:
    {
        std::ostringstream a, b;
        plot.clip(0.5);
        a << plot;
        plot.clip(5);
        b << plot;
        CHECK_TRUE(a.str() != b.str());
        CHECK_TRUE(first_mask_id(a.str()) != first_mask_id(b.str()));
        std::ostringstream c;
        plot.envelope_color("red");
        c << plot;
        CHECK_TRUE(first_mask_id(b.str()) != first_mask_id(c.str()));
    }

    // The summary counts the grid and the probes:
    CHECK_LE(std::size_t(1), count(plot.summary(), "samples=500 "));

    // A logarithmic axis:
    ulps_plot<decltype(f), double, float> log_plot(f, 0.01f, 100.0f, 200, true, 3, true);
    log_plot.add_fn(g).png_encoder(nullptr);
    std::ostringstream os;
    os << log_plot;
    check_svg(os.str());

    // The heatmap, with probes:
    auto f2 = [](double x, double y) { return std::sin(x)*std::exp(y); };
    auto g2 = [](float x, float y) { return std::sin(x)*std::exp(y); };
    widths.clear();
    heights.clear();
    ulps_heatmap<decltype(f2), double, float> heat(f2, 0.5f, 6.0f, -1.0f, 1.0f, 40, 30, false, false, 2);
    heat.add_fn(g2, "g2").png_encoder(recording_encoder{&widths, &heights, &gray});
    std::string const heat_name = "test_ulps_heatmap.svg";
    heat.write_ulps_envelope(heat_name);
    check_svg(read_file(heat_name));
    CHECK_EQUAL(widths.size(), std::size_t(1));
    CHECK_EQUAL(widths[0], std::size_t(40));
    CHECK_EQUAL(heights[0], std::size_t(30));
    std::remove(heat_name.c_str());
    CHECK_LE(std::size_t(1), count(heat.summary(), "1200 pixels (300 probes)"));
}

int main()
{
    test_condition_numbers_1d();
    test_condition_numbers_2d();
    test_raster_canvas();
    test_files();
    return boost::math::test::report_errors();
}
