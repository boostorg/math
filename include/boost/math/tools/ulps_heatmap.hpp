//  (C) Copyright Nick Thompson 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Two-dimensional ulps plots: heatmaps of the error of a function f(x, y) of two
// variables over a rectangle, the 2D analogue of boost/math/tools/ulps_plot.hpp.
//
// Two plots can be written from the same data:
//   write_ulps:          max(1/2, |ulps|). A correctly rounded implementation stays at
//                        1/2, whatever the conditioning of f.
//   write_ulps_envelope: max(1, |ulps|/max(1/2, cond/2)), where cond = |x f_x/f| + |y f_y/f|
//                        is the condition number of evaluating f with both arguments
//                        rounded. Values above 1 are errors the arguments' rounding can't
//                        explain.
// Both use a sequential colormap (inferno by default), so the floor is dark and errors
// glow. Each can use a linear or a logarithmic color scale.
//
// The heatmap is a PNG embedded in an SVG, which carries the axes, colorbar and labels.
// The PNG encoder can be set with png_encoder. If lodepng.h is on the include path, the
// default encoder is lodepng (https://github.com/lvandeve/lodepng; link lodepng.cpp).

#ifndef BOOST_MATH_TOOLS_ULPS_HEATMAP_HPP
#define BOOST_MATH_TOOLS_ULPS_HEATMAP_HPP
#if (__cplusplus < 201703L) && (!defined(_MSVC_LANG) || (_MSVC_LANG < 201703L))
#error "The header <boost/math/tools/ulps_heatmap.hpp> can only be used in C++17 and later."
#endif
#ifndef BOOST_MATH_BUILD_MODULE
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>
#endif
#include <boost/math/tools/color_maps.hpp>
#include <boost/math/tools/ulps_plot.hpp>

#if __has_include("lodepng.h")
#include "lodepng.h"
#define BOOST_MATH_ULPS_HEATMAP_HAS_LODEPNG
#endif

BOOST_MATH_NAMESPACE_BEGIN namespace tools {

namespace detail {


// The color_maps.hpp tables are already sRGB, so they need no gamma correction:
inline std::array<std::uint8_t, 3> to_8bit_rgb(std::array<double, 3> const & c)
{
    std::array<std::uint8_t, 3> rgb;
    for (std::size_t k = 0; k < 3; ++k)
    {
        rgb[k] = static_cast<std::uint8_t>(std::lround(255*std::clamp(c[k], 0.0, 1.0)));
    }
    return rgb;
}

inline std::string base64(std::vector<unsigned char> const & bytes)
{
    static constexpr char table[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    std::string out;
    out.reserve(4*((bytes.size() + 2)/3));
    for (std::size_t i = 0; i < bytes.size(); i += 3)
    {
        std::uint32_t n = std::uint32_t(bytes[i]) << 16;
        if (i + 1 < bytes.size()) { n |= std::uint32_t(bytes[i + 1]) << 8; }
        if (i + 2 < bytes.size()) { n |= bytes[i + 2]; }
        out += table[(n >> 18) & 63];
        out += table[(n >> 12) & 63];
        out += i + 1 < bytes.size() ? table[(n >> 6) & 63] : '=';
        out += i + 2 < bytes.size() ? table[n & 63] : '=';
    }
    return out;
}

inline std::function<std::vector<unsigned char>(std::vector<unsigned char> const &, std::size_t, std::size_t)> default_png_encoder()
{
#ifdef BOOST_MATH_ULPS_HEATMAP_HAS_LODEPNG
    return [](std::vector<unsigned char> const & rgb, std::size_t width, std::size_t height) {
        std::vector<unsigned char> png;
        unsigned error = lodepng::encode(png, rgb, static_cast<unsigned>(width), static_cast<unsigned>(height), LCT_RGB);
        if (error)
        {
            throw std::runtime_error(std::string("lodepng: ") + lodepng_error_text(error));
        }
        return png;
    };
#else
    return nullptr;
#endif
}

// A diverging map for signed data on a dark background: light blue -> dark blue -> near black (zero) -> dark red -> light orange.
// Nothing is glaringly bright, which the cool-warm maps are. Similar in spirit to Crameri's "berlin".
inline std::array<double, 3> dark_diverging(double t)
{
    static constexpr double stops[5][3] = {{0.45, 0.65, 0.95}, {0.10, 0.30, 0.60}, {0.06, 0.05, 0.07}, {0.55, 0.18, 0.10}, {0.95, 0.62, 0.38}};
    t = std::clamp(t, 0.0, 1.0)*4;
    int i = (std::min)(3, static_cast<int>(t));
    double u = t - i;
    return {stops[i][0]*(1 - u) + stops[i + 1][0]*u, stops[i][1]*(1 - u) + stops[i + 1][1]*u, stops[i][2]*(1 - u) + stops[i + 1][2]*u};
}

// Renders a tick label; whole powers of ten on a logarithmic axis are written as 10^k.
inline std::string tick_label(double v, bool log_axis)
{
    std::ostringstream os;
    if (log_axis && v > 0)
    {
        double e = std::round(std::log10(v));
        if (std::abs(v - std::pow(10.0, e)) <= 1e-9*v)
        {
            os << "10<tspan dy='-6' font-size='10'>" << static_cast<int>(e) << "</tspan>";
            return os.str();
        }
    }
    os << std::setprecision(4) << v;
    return os.str();
}

// Approximate width in pixels of a tick label at font-size 13, ignoring its SVG markup:
inline int label_width(std::string const & label)
{
    int chars = 0;
    bool in_tag = false;
    for (char c : label)
    {
        if (c == '<') { in_tag = true; }
        else if (c == '>') { in_tag = false; }
        else if (!in_tag) { ++chars; }
    }
    return 7*chars;
}

// Labels v > 0 on a logarithmic colorbar as a power of ten, 10^2 or 10^2.45, so that every label on it reads alike.
inline std::string power_label(double v)
{
    double e = std::log10(v);
    std::ostringstream os;
    if (std::abs(e - std::round(e)) < 0.005)
    {
        os << static_cast<int>(std::round(e));
    }
    else
    {
        os << std::fixed << std::setprecision(2) << e;
    }
    return "10<tspan dy='-6' font-size='10'>" + os.str() + "</tspan>";
}

// Coarse abscissas at the centres of n pixels spanning [a, b], evenly spaced in x or log(x).
template<class CoarseReal>
std::vector<CoarseReal> heatmap_abscissas(CoarseReal a, CoarseReal b, std::size_t n, bool log_axis)
{
    std::vector<CoarseReal> xs(n);
    long double lo = log_axis ? std::log(static_cast<long double>(a)) : static_cast<long double>(a);
    long double hi = log_axis ? std::log(static_cast<long double>(b)) : static_cast<long double>(b);
    for (std::size_t i = 0; i < n; ++i)
    {
        long double u = lo + (hi - lo)*(i + 0.5L)/n;
        xs[i] = static_cast<CoarseReal>(log_axis ? std::exp(u) : u);
    }
    return xs;
}

// |x f'(x)/f(x)|, given fx = f(x), by a one-sided difference with relative step sqrt(eps). The envelope
// needs only a few digits, which this gives with one more evaluation of f; evaluation_condition_number
// takes six. Steps backward if the forward step leaves the domain.
template<class F, class Real>
Real heatmap_condition_number(F f, Real const & x, Real const & fx)
{
    using std::abs;
    using std::isfinite;
    using std::sqrt;
    if (x == 0)
    {
        return 0;
    }
    Real h = sqrt(std::numeric_limits<Real>::epsilon());
    for (Real step : {h, -h})
    {
        try
        {
            Real xh = x + x*step;
            Real slope = (f(xh) - fx)/(xh - x);
            if (isfinite(slope))
            {
                return abs(x*slope/fx);
            }
        }
        catch (...)
        {
        }
    }
    return std::numeric_limits<Real>::quiet_NaN();
}

// Runs body(j) for j in [0, n) on several threads.
template<class Body>
void parallel_for(std::size_t n, unsigned threads, Body body)
{
    if (threads == 0)
    {
        threads = (std::max)(1u, std::thread::hardware_concurrency());
    }
    std::atomic<std::size_t> next{0};
    std::vector<std::thread> pool;
    for (unsigned t = 0; t < threads; ++t)
    {
        pool.emplace_back([&]() {
            for (std::size_t j = next++; j < n; j = next++)
            {
                body(j);
            }
        });
    }
    for (auto & th : pool)
    {
        th.join();
    }
}

} // namespace detail

template<class F, typename PreciseReal, typename CoarseReal>
class ulps_heatmap {
public:
    // Samples f at the pixel centres of a width x height grid on [x_min, x_max] x [y_min, y_max].
    // A logarithmic axis is sampled evenly in log and requires a positive minimum.
    // hi_acc_impl(x, y) takes PreciseReal arguments and must be safe to call from several threads.
    ulps_heatmap(F hi_acc_impl, CoarseReal x_min, CoarseReal x_max, CoarseReal y_min, CoarseReal y_max,
                 std::size_t width = 800, std::size_t height = 500,
                 bool log_x = false, bool log_y = false, unsigned threads = 0);

    // Computes the error of the implementation g(x, y), which takes CoarseReal arguments. Each
    // function added is drawn as its own panel, side by side with the others and on the same scale.
    template<class G>
    ulps_heatmap& add_fn(G g, std::string const & label = "");

    // Saturates the color scale at clip. By default the clip is chosen from the values in all
    // panels: the largest on a logarithmic scale, and the 99th percentile on a linear one.
    ulps_heatmap& clip(double clip);

    // The color scale runs from the floor (1/2 for ulps, 1 for ulps/envelope) to the clip,
    // linearly or logarithmically.
    ulps_heatmap& log_scale(bool log_scale);

    // Maps [0, 1] to sRGB in [0, 1]^3; inferno by default. Any of the color_maps.hpp maps fit.
    ulps_heatmap& color_map(std::function<std::array<double, 3>(double)> color_map);

    // Encodes 8-bit RGB pixels, row-major from the top row, as a PNG; lodepng by default, if lodepng.h is found.
    using png_encoder_type = std::function<std::vector<unsigned char>(std::vector<unsigned char> const & rgb, std::size_t width, std::size_t height)>;
    ulps_heatmap& png_encoder(png_encoder_type encoder);

    ulps_heatmap& title(std::string const & title);

    ulps_heatmap& x_label(std::string const & label);

    ulps_heatmap& y_label(std::string const & label);

    ulps_heatmap& background_color(std::string const & color);

    ulps_heatmap& font_color(std::string const & color);

    // Draws the function itself (the reference values) as a heatmap below each error panel, on the same axes
    // and with its own colorbar. Off by default. The color mapping is chosen from the values; see write_svg.
    ulps_heatmap& show_function(bool show_function);

    // The caption of the function panels; by default the title, or f(x, y) if there is none.
    ulps_heatmap& function_label(std::string const & label);

    // The function panels are colored linearly, from min(f) to max(f) (symmetric about 0, with a dark center, if f has both signs).
    // If f > 0 everywhere and spans three or more decades, a logarithmic scale is used automatically; log_function(true) forces it
    // when f > 0 everywhere, and log_function(false) turns it off.
    ulps_heatmap& log_function(bool log_function);

    void write_ulps(std::string const & filename) const;

    void write_ulps_envelope(std::string const & filename) const;

    // For each panel, counts of pixels without a value and outside 1/2 ulp and the envelope,
    // and the worst errors and where they occur.
    std::string summary() const;

private:
    struct panel
    {
        std::string label;
        // Row-major, with row 0 at y_max (the top of the image):
        std::vector<double> ulps;
        std::vector<double> ulps_envelope;
        std::vector<detail::ulps_class> state;
    };

    void write_svg(std::string const & filename, std::vector<double> panel::* values, double floor,
                   std::string const & quantity) const;

    std::vector<CoarseReal> xs_;
    std::vector<CoarseReal> ys_;
    // Row-major, with row 0 at y_max (the top of the image):
    std::vector<PreciseReal> precise_;
    std::vector<double> envelope_;
    std::vector<panel> panels_;
    CoarseReal x_min_;
    CoarseReal x_max_;
    CoarseReal y_min_;
    CoarseReal y_max_;
    bool log_x_;
    bool log_y_;
    unsigned threads_;
    double clip_ = -1;
    bool log_scale_ = false;
    std::function<std::array<double, 3>(double)> color_map_ = inferno<double>;
    png_encoder_type png_encoder_ = detail::default_png_encoder();
    std::string title_;
    std::string x_label_ = "x";
    std::string y_label_ = "y";
    std::string background_color_ = "black";
    std::string font_color_ = "white";
    bool show_function_ = false;
    std::string function_label_;
    int log_function_ = -1;
};

template<class F, typename PreciseReal, typename CoarseReal>
ulps_heatmap<F, PreciseReal, CoarseReal>::ulps_heatmap(F hi_acc_impl, CoarseReal x_min, CoarseReal x_max,
             CoarseReal y_min, CoarseReal y_max, std::size_t width, std::size_t height,
             bool log_x, bool log_y, unsigned threads)
    : x_min_(x_min), x_max_(x_max), y_min_(y_min), y_max_(y_max), log_x_(log_x), log_y_(log_y), threads_(threads)
{
    static_assert(std::numeric_limits<PreciseReal>::digits10 >= std::numeric_limits<CoarseReal>::digits10, "PreciseReal must have higher precision that CoarseReal");
    if (!(x_min < x_max) || !(y_min < y_max))
    {
        throw std::domain_error("Require x_min < x_max and y_min < y_max.");
    }
    if ((log_x && !(x_min > 0)) || (log_y && !(y_min > 0)))
    {
        throw std::domain_error("A logarithmic axis requires a positive minimum.");
    }
    if (width < 2 || height < 2)
    {
        throw std::domain_error("The heatmap must be at least 2 x 2 pixels.");
    }
    xs_ = detail::heatmap_abscissas(x_min, x_max, width, log_x);
    ys_ = detail::heatmap_abscissas(y_min, y_max, height, log_y);
    std::reverse(ys_.begin(), ys_.end());

    precise_.resize(width*height);
    envelope_.resize(width*height, std::numeric_limits<double>::quiet_NaN());
    detail::parallel_for(height, threads_, [&](std::size_t j) {
        using std::abs;
        PreciseReal y = ys_[j];
        for (std::size_t i = 0; i < width; ++i)
        {
            PreciseReal x = xs_[i];
            PreciseReal & z = precise_[j*width + i];
            try
            {
                z = hi_acc_impl(x, y);
                // Relative perturbations of x and y each contribute their own condition number:
                auto f_x = [&](PreciseReal const & t) { return hi_acc_impl(t, y); };
                auto f_y = [&](PreciseReal const & t) { return hi_acc_impl(x, t); };
                PreciseReal cond = detail::heatmap_condition_number(f_x, x, z) + detail::heatmap_condition_number(f_y, y, z);
                // As in ulps_plot, the envelope never drops below the 1/2 ulp of correct rounding:
                envelope_[j*width + i] = (std::max)(0.5, static_cast<double>(cond/2));
            }
            catch (...)
            {
                z = std::numeric_limits<PreciseReal>::quiet_NaN();
            }
        }
    });
}

template<class F, typename PreciseReal, typename CoarseReal>
template<class G>
ulps_heatmap<F, PreciseReal, CoarseReal>& ulps_heatmap<F, PreciseReal, CoarseReal>::add_fn(G g, std::string const & label)
{
    std::size_t width = xs_.size();
    std::size_t n = precise_.size();
    panel p;
    p.label = label;
    p.ulps.assign(n, std::numeric_limits<double>::quiet_NaN());
    p.ulps_envelope.assign(n, std::numeric_limits<double>::quiet_NaN());
    p.state.assign(n, detail::ulps_class::ok);
    detail::parallel_for(ys_.size(), threads_, [&](std::size_t j) {
        using std::abs;
        using std::isfinite;
        using std::nextafter;
        for (std::size_t i = 0; i < width; ++i)
        {
            std::size_t k = j*width + i;
            PreciseReal const & z = precise_[k];
            // The implementation is evaluated everywhere, even where the true value is out of range, to check
            // that overflow and underflow are handled correctly:
            CoarseReal w = 0;
            detail::ulps_class c = detail::classify_result<CoarseReal>(z, [&]() { return g(xs_[i], ys_[j]); }, w);
            p.state[k] = c;
            if (c != detail::ulps_class::ok)
            {
                continue;
            }
            // The same ulp distance as ulps_plot:
            CoarseReal absz = static_cast<CoarseReal>(abs(z));
            PreciseReal dist = static_cast<PreciseReal>(nextafter(absz, (std::numeric_limits<CoarseReal>::max)()) - absz);
            p.ulps[k] = static_cast<double>((static_cast<PreciseReal>(w) - z)/dist);
            p.ulps_envelope[k] = p.ulps[k]/envelope_[k];
        }
    });
    panels_.push_back(std::move(p));
    return *this;
}

template<class F, typename PreciseReal, typename CoarseReal>
ulps_heatmap<F, PreciseReal, CoarseReal>& ulps_heatmap<F, PreciseReal, CoarseReal>::clip(double clip)
{
    clip_ = clip;
    return *this;
}

template<class F, typename PreciseReal, typename CoarseReal>
ulps_heatmap<F, PreciseReal, CoarseReal>& ulps_heatmap<F, PreciseReal, CoarseReal>::log_scale(bool log_scale)
{
    log_scale_ = log_scale;
    return *this;
}

template<class F, typename PreciseReal, typename CoarseReal>
ulps_heatmap<F, PreciseReal, CoarseReal>& ulps_heatmap<F, PreciseReal, CoarseReal>::color_map(std::function<std::array<double, 3>(double)> color_map)
{
    color_map_ = color_map;
    return *this;
}

template<class F, typename PreciseReal, typename CoarseReal>
ulps_heatmap<F, PreciseReal, CoarseReal>& ulps_heatmap<F, PreciseReal, CoarseReal>::png_encoder(png_encoder_type encoder)
{
    png_encoder_ = encoder;
    return *this;
}

template<class F, typename PreciseReal, typename CoarseReal>
ulps_heatmap<F, PreciseReal, CoarseReal>& ulps_heatmap<F, PreciseReal, CoarseReal>::title(std::string const & title)
{
    title_ = title;
    return *this;
}

template<class F, typename PreciseReal, typename CoarseReal>
ulps_heatmap<F, PreciseReal, CoarseReal>& ulps_heatmap<F, PreciseReal, CoarseReal>::x_label(std::string const & label)
{
    x_label_ = label;
    return *this;
}

template<class F, typename PreciseReal, typename CoarseReal>
ulps_heatmap<F, PreciseReal, CoarseReal>& ulps_heatmap<F, PreciseReal, CoarseReal>::y_label(std::string const & label)
{
    y_label_ = label;
    return *this;
}

template<class F, typename PreciseReal, typename CoarseReal>
ulps_heatmap<F, PreciseReal, CoarseReal>& ulps_heatmap<F, PreciseReal, CoarseReal>::background_color(std::string const & color)
{
    background_color_ = color;
    return *this;
}

template<class F, typename PreciseReal, typename CoarseReal>
ulps_heatmap<F, PreciseReal, CoarseReal>& ulps_heatmap<F, PreciseReal, CoarseReal>::font_color(std::string const & color)
{
    font_color_ = color;
    return *this;
}

template<class F, typename PreciseReal, typename CoarseReal>
ulps_heatmap<F, PreciseReal, CoarseReal>& ulps_heatmap<F, PreciseReal, CoarseReal>::show_function(bool show_function)
{
    show_function_ = show_function;
    return *this;
}

template<class F, typename PreciseReal, typename CoarseReal>
ulps_heatmap<F, PreciseReal, CoarseReal>& ulps_heatmap<F, PreciseReal, CoarseReal>::log_function(bool log_function)
{
    log_function_ = log_function ? 1 : 0;
    return *this;
}

template<class F, typename PreciseReal, typename CoarseReal>
ulps_heatmap<F, PreciseReal, CoarseReal>& ulps_heatmap<F, PreciseReal, CoarseReal>::function_label(std::string const & label)
{
    function_label_ = label;
    return *this;
}

template<class F, typename PreciseReal, typename CoarseReal>
void ulps_heatmap<F, PreciseReal, CoarseReal>::write_ulps(std::string const & filename) const
{
    write_svg(filename, &panel::ulps, 0.5, "|ulps|");
}

template<class F, typename PreciseReal, typename CoarseReal>
void ulps_heatmap<F, PreciseReal, CoarseReal>::write_ulps_envelope(std::string const & filename) const
{
    write_svg(filename, &panel::ulps_envelope, 1, "|ulps| / max(1/2, cond/2)");
}

template<class F, typename PreciseReal, typename CoarseReal>
std::string ulps_heatmap<F, PreciseReal, CoarseReal>::summary() const
{
    std::ostringstream os;
    std::size_t width = xs_.size();
    for (panel const & p : panels_)
    {
        using detail::ulps_class;
        std::size_t counts[7] = {0, 0, 0, 0, 0, 0, 0};
        std::size_t first[7] = {0, 0, 0, 0, 0, 0, 0};
        std::size_t bx0[7], bx1[7], by0[7], by1[7];
        for (int c = 0; c < 7; ++c) { bx0[c] = width; bx1[c] = 0; by0[c] = ys_.size(); by1[c] = 0; }
        std::size_t over_half = 0;
        std::size_t over_envelope = 0;
        std::size_t worst_ulps = 0;
        std::size_t worst_ratio = 0;
        for (std::size_t k = 0; k < p.state.size(); ++k)
        {
            std::size_t c = static_cast<std::size_t>(p.state[k]);
            if (counts[c]++ == 0) { first[c] = k; }
            bx0[c] = (std::min)(bx0[c], k % width); bx1[c] = (std::max)(bx1[c], k % width);
            by0[c] = (std::min)(by0[c], k / width); by1[c] = (std::max)(by1[c], k / width);
            if (p.state[k] != ulps_class::ok) { continue; }
            if (std::abs(p.ulps[k]) > 0.5) { ++over_half; }
            if (std::abs(p.ulps_envelope[k]) > 1) { ++over_envelope; }
            if (!(std::abs(p.ulps[k]) <= std::abs(p.ulps[worst_ulps]))) { worst_ulps = k; }
            if (!(std::abs(p.ulps_envelope[k]) <= std::abs(p.ulps_envelope[worst_ratio]))) { worst_ratio = k; }
        }
        if (!p.label.empty())
        {
            os << p.label << ": ";
        }
        os << p.state.size() << " pixels: ok=" << counts[0] << " undefined=" << counts[1] << " overflow_handled=" << counts[2]
           << " underflow_handled=" << counts[3] << " green_nan_or_domain_error=" << counts[4] << " cyan_spurious_overflow=" << counts[5]
           << " mishandled=" << counts[6] << " over_half_ulp=" << over_half << " over_envelope=" << over_envelope << "\n";
        char const * names[7] = {"", "", "", "", "Green (NaN or domain error)", "Cyan (spurious overflow)", "Mishandled out-of-range result (maximal error)"};
        for (std::size_t c = 4; c < 7; ++c)
        {
            if (counts[c] > 0)
            {
                os << "  " << names[c] << ": first at (" << std::setprecision(std::numeric_limits<CoarseReal>::max_digits10)
                   << xs_[first[c] % width] << ", " << ys_[first[c] / width] << "); x in [" << xs_[bx0[c]] << ", " << xs_[bx1[c]]
                   << "], y in [" << ys_[by1[c]] << ", " << ys_[by0[c]] << "]\n";
            }
        }
        os << "  Worst ulps:          " << std::setprecision(4) << p.ulps[worst_ulps]
           << std::setprecision(std::numeric_limits<CoarseReal>::max_digits10)
           << " at (" << xs_[worst_ulps % width] << ", " << ys_[worst_ulps / width] << ")\n";
        os << "  Worst ulps/envelope: " << std::setprecision(4) << p.ulps_envelope[worst_ratio]
           << std::setprecision(std::numeric_limits<CoarseReal>::max_digits10)
           << " at (" << xs_[worst_ratio % width] << ", " << ys_[worst_ratio / width] << ")\n";
    }
    return os.str();
}

template<class F, typename PreciseReal, typename CoarseReal>
void ulps_heatmap<F, PreciseReal, CoarseReal>::write_svg(std::string const & filename, std::vector<double> panel::* values,
                                                         double floor, std::string const & quantity) const
{
    if (panels_.empty())
    {
        throw std::domain_error("No function added; call add_fn first.");
    }
    if (!detail::ends_with(filename, ".svg"))
    {
        throw std::logic_error("Only svg files are supported at this time.");
    }
    if (!png_encoder_)
    {
        throw std::logic_error("No PNG encoder: set one with png_encoder, or put lodepng.h on the include path.");
    }
    // Gray: correctly handled: the true value is undefined, or out of range and the implementation overflowed (or underflowed to zero
    //       or a subnormal) as it should. Green: NaN or a domain error where the true value is a finite double.
    static constexpr std::array<std::uint8_t, 3> no_reference_rgb = {70, 70, 70};
    static constexpr std::array<std::uint8_t, 3> failed_rgb = {127, 255, 0};
    // Cyan: the implementation threw std::overflow_error or returned an infinity, where the true value is a finite double.
    static constexpr std::array<std::uint8_t, 3> spurious_overflow_rgb = {0, 229, 255};

    // Plot max(floor, |value|), so that log scales work:
    auto magnitude = [&](double v) { return (std::max)(floor, std::abs(v)); };
    double clip = clip_;
    if (clip <= 0)
    {
        // Choose the clip from the values of every panel, so that they share a scale:
        std::vector<double> magnitudes;
        for (panel const & p : panels_)
        {
            for (std::size_t k = 0; k < p.state.size(); ++k)
            {
                if (p.state[k] == detail::ulps_class::ok && std::isfinite((p.*values)[k]))
                {
                    magnitudes.push_back(magnitude((p.*values)[k]));
                }
            }
        }
        // A few huge errors would wash out a linear scale, so clip it at the 99th percentile:
        if (!magnitudes.empty())
        {
            std::size_t rank = log_scale_ ? magnitudes.size() - 1 : (99*magnitudes.size())/100;
            std::nth_element(magnitudes.begin(), magnitudes.begin() + rank, magnitudes.end());
            clip = magnitudes[rank];
        }
    }
    // The scale needs some room above the floor:
    if (!(clip > floor*1.01))
    {
        clip = log_scale_ ? 10*floor : 2*floor;
    }
    // Maps a value to [0, 1]:
    auto scale = [&](double v) {
        double m = magnitude(v);
        double t = log_scale_ ? std::log(m/floor)/std::log(clip/floor) : (m - floor)/(clip - floor);
        return (std::min)(1.0, t);
    };

    std::size_t width = xs_.size();
    std::size_t height = ys_.size();
    bool labelled = std::any_of(panels_.begin(), panels_.end(), [](panel const & p) { return !p.label.empty(); });
    int const margin_left = 80;
    int const margin_top = (title_.empty() ? 20 : 50) + (labelled ? 24 : 0);
    int const margin_bottom = 60;
    int const panel_gap = 30;
    int const bar_gap = 25;
    int const bar_width = 22;
    int const margin_right = 215;
    int const W = static_cast<int>(width);
    int const H = static_cast<int>(height);
    int const panels = static_cast<int>(panels_.size());
    int const plots_width = panels*W + (panels - 1)*panel_gap;
    int const total_width = margin_left + plots_width + margin_right;
    // As many axis ticks as fit without their labels overlapping:
    int const x_ticks = (std::max)(2, (std::min)(8, W/60));
    int const y_ticks = (std::max)(2, (std::min)(6, H/40));
    // The function heatmaps go below the error heatmaps, separated by room for the latter's x tick labels:
    int const function_gap = 66;
    int const function_top = margin_top + H + function_gap;
    int const total_height = show_function_ ? function_top + H + margin_bottom : margin_top + H + margin_bottom;
    auto x_pos = [&](double x) {
        double u = log_x_ ? std::log(x/x_min_)/std::log(double(x_max_)/x_min_) : (x - x_min_)/(double(x_max_) - x_min_);
        return u*W;
    };
    auto y_pos = [&](double y) {
        double u = log_y_ ? std::log(y/y_min_)/std::log(double(y_max_)/y_min_) : (y - y_min_)/(double(y_max_) - y_min_);
        return margin_top + (1 - u)*H;
    };
    auto rgb_string = [](std::array<std::uint8_t, 3> const & c) {
        return "rgb(" + std::to_string(c[0]) + "," + std::to_string(c[1]) + "," + std::to_string(c[2]) + ")";
    };
    // Axis ticks; reuse ulps_plot's gridline abscissas, plus the minimum on a linear axis:
    auto ticks = [](CoarseReal a, CoarseReal b, int lines, bool log_axis) {
        std::vector<CoarseReal> t = detail::vertical_gridlines(a, b, lines, log_axis);
        if (!log_axis)
        {
            t.insert(t.begin(), a);
        }
        return t;
    };

    // The function's color mapping: f itself, linearly from min(f) to max(f) with a sequential map (viridis); if f has both signs,
    // symmetric about zero with a dark-centered diverging map. If f > 0 everywhere and spans three or more decades, color by log10(f)
    // (see log_function). If f has both signs and its largest values dwarf the bulk of the plot (0F1(b, z) reaches 1e9 in a corner
    // and is of order 1 elsewhere), a linear map would show only black, so color by sgn(f) ln(1 + |f|/s), with s chosen so that the
    // 90th percentile of |f| falls halfway between s and max |f| on a log scale. Pixels with no finite reference in CoarseReal's range are
    // drawn gray.
    enum class fn_mode { log_seq, linear_seq, linear_div, symlog_div };
    fn_mode mode = fn_mode::linear_seq;
    double fmin = std::numeric_limits<double>::infinity();
    double fmax = -std::numeric_limits<double>::infinity();
    double fsym = 0;
    double fscale = 1;
    std::vector<unsigned char> fn_png;
    std::vector<unsigned char> const * fn_png_ptr = nullptr;
    auto fn_ok = [&](std::size_t k) {
        using std::abs;
        using std::isfinite;
        return isfinite(precise_[k]) && abs(precise_[k]) <= (std::numeric_limits<CoarseReal>::max)();
    };
    auto fn_t = [&](double v) {
        double t = 0;
        switch (mode)
        {
        case fn_mode::log_seq:     t = (std::log10(v) - std::log10(fmin))/(std::log10(fmax) - std::log10(fmin)); break;
        case fn_mode::linear_seq:  t = (v - fmin)/(fmax - fmin); break;
        case fn_mode::linear_div:  t = 0.5 + 0.5*v/fsym; break;
        case fn_mode::symlog_div:  t = 0.5 + 0.5*std::copysign(std::log1p(std::abs(v)/fscale), v)/std::log1p(fsym/fscale); break;
        }
        return std::clamp(t, 0.0, 1.0);
    };
    if (show_function_)
    {
        for (std::size_t k = 0; k < width*height; ++k)
        {
            if (fn_ok(k))
            {
                double v = static_cast<double>(precise_[k]);
                fmin = (std::min)(fmin, v);
                fmax = (std::max)(fmax, v);
            }
        }
        if (!(fmin <= fmax))
        {
            fmin = 0;
            fmax = 1;
        }
        if (fmin > 0 && (log_function_ == 1 || (log_function_ == -1 && fmax >= 1000*fmin)) && fmax > fmin)
        {
            mode = fn_mode::log_seq;
        }
        else if (fmin < 0 && fmax > 0)
        {
            mode = fn_mode::linear_div;
            fsym = (std::max)(-fmin, fmax);
            std::vector<double> fabs_values;
            for (std::size_t k = 0; k < width*height; ++k)
            {
                if (fn_ok(k)) { fabs_values.push_back(std::abs(static_cast<double>(precise_[k]))); }
            }
            auto p90 = fabs_values.begin() + static_cast<std::ptrdiff_t>(0.9*(fabs_values.size() - 1));
            std::nth_element(fabs_values.begin(), p90, fabs_values.end());
            if (*p90 > 0 && fsym > 1000 * *p90)
            {
                mode = fn_mode::symlog_div;
                fscale = (std::max)(*p90 * (*p90 / fsym), (std::numeric_limits<double>::min)());
            }
        }
        else
        {
            mode = fn_mode::linear_seq;
            if (!(fmax > fmin))
            {
                // A constant function:
                fmin -= 0.5*std::abs(fmin) + 0.5;
                fmax += 0.5*std::abs(fmax) + 0.5;
            }
        }
        std::vector<unsigned char> rgb(3*width*height);
        for (std::size_t k = 0; k < width*height; ++k)
        {
            std::array<std::uint8_t, 3> c = no_reference_rgb;
            if (fn_ok(k))
            {
                double t = fn_t(static_cast<double>(precise_[k]));
                c = detail::to_8bit_rgb(mode == fn_mode::linear_div || mode == fn_mode::symlog_div ? detail::dark_diverging(t) : viridis<double>(t));
            }
            std::copy(c.begin(), c.end(), rgb.begin() + 3*k);
        }
        fn_png = png_encoder_(rgb, width, height);
        fn_png_ptr = &fn_png;
    }

    std::ofstream fs(filename);
    fs << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n"
       << "<svg xmlns='http://www.w3.org/2000/svg' width='" << total_width << "' height='" << total_height << "'>\n"
       << "<style>svg { font-family: times; fill: " << font_color_ << "; }</style>\n"
       << "<rect width='100%' height='100%' fill='" << background_color_ << "'/>\n";
    if (!title_.empty())
    {
        fs << "<text x='" << margin_left + plots_width/2 << "' y='31' text-anchor='middle' font-size='18'>" << title_ << "</text>\n";
    }
    bool any_failed = false;
    bool any_spurious = false;
    bool any_no_reference = false;
    for (int q = 0; q < panels; ++q)
    {
        panel const & p = panels_[q];
        int const left = margin_left + q*(W + panel_gap);
        std::vector<unsigned char> rgb(3*width*height);
        for (std::size_t k = 0; k < width*height; ++k)
        {
            double v = (p.*values)[k];
            std::array<std::uint8_t, 3> c = no_reference_rgb;
            switch (p.state[k])
            {
            case detail::ulps_class::ok:
                // An envelope can be missing (NaN) where the condition number could not be estimated; draw it at the floor.
                c = detail::to_8bit_rgb(color_map_(std::isfinite(v) ? scale(v) : 0.0));
                break;
            case detail::ulps_class::failed:
                c = failed_rgb;
                any_failed = true;
                break;
            case detail::ulps_class::spurious_overflow:
                c = spurious_overflow_rgb;
                any_spurious = true;
                break;
            case detail::ulps_class::mishandled:
                // An out-of-range result handled wrongly is the maximal error: the top of the color scale.
                c = detail::to_8bit_rgb(color_map_(1.0));
                break;
            default:
                any_no_reference = true;
                break;
            }
            std::copy(c.begin(), c.end(), rgb.begin() + 3*k);
        }
        std::vector<unsigned char> png = png_encoder_(rgb, width, height);
        if (!p.label.empty())
        {
            fs << "<text x='" << left + W/2 << "' y='" << margin_top - 8 << "' text-anchor='middle' font-size='15'>" << p.label << "</text>\n";
        }
        fs << "<image x='" << left << "' y='" << margin_top << "' width='" << W << "' height='" << H
           << "' style='image-rendering: pixelated' href='data:image/png;base64," << detail::base64(png) << "'/>\n";
        fs << "<rect x='" << left << "' y='" << margin_top << "' width='" << W << "' height='" << H
           << "' fill='none' stroke='gray'/>\n";
        for (CoarseReal x : ticks(x_min_, x_max_, x_ticks, log_x_))
        {
            double px = left + x_pos(x);
            fs << "<line x1='" << px << "' y1='" << margin_top + H << "' x2='" << px << "' y2='" << margin_top + H + 5
               << "' stroke='gray'/>\n"
               << "<text x='" << px << "' y='" << margin_top + H + 20 << "' text-anchor='middle' font-size='13'>"
               << (log_x_ ? detail::tick_label(x, true) : detail::axis_label(static_cast<double>(x), static_cast<double>(x_min_), static_cast<double>(x_max_), x_ticks, false)) << "</text>\n";
        }
        if (!show_function_)
        {
            fs << "<text x='" << left + W/2 << "' y='" << margin_top + H + 45 << "' text-anchor='middle' font-size='15'>"
               << x_label_ << "</text>\n";
        }
        else
        {
            fs << "<image x='" << left << "' y='" << function_top << "' width='" << W << "' height='" << H
               << "' style='image-rendering: pixelated' href='data:image/png;base64," << detail::base64(*fn_png_ptr) << "'/>\n";
            fs << "<rect x='" << left << "' y='" << function_top << "' width='" << W << "' height='" << H
               << "' fill='none' stroke='gray'/>\n";
            for (CoarseReal x : ticks(x_min_, x_max_, x_ticks, log_x_))
            {
                double px = left + x_pos(x);
                fs << "<line x1='" << px << "' y1='" << function_top + H << "' x2='" << px << "' y2='" << function_top + H + 5
                   << "' stroke='gray'/>\n"
                   << "<text x='" << px << "' y='" << function_top + H + 20 << "' text-anchor='middle' font-size='13'>"
                   << (log_x_ ? detail::tick_label(x, true) : detail::axis_label(static_cast<double>(x), static_cast<double>(x_min_), static_cast<double>(x_max_), x_ticks, false)) << "</text>\n";
            }
            fs << "<text x='" << left + W/2 << "' y='" << function_top + H + 45 << "' text-anchor='middle' font-size='15'>"
               << x_label_ << "</text>\n";
            fs << "<text x='" << left + W/2 << "' y='" << function_top - 8 << "' text-anchor='middle' font-size='13'>" << (!function_label_.empty() ? function_label_ : (!title_.empty() ? title_ : std::string("f(x, y)"))) << "</text>\n";
        }
    }
    // The panels share the y axis, which is labelled on the first:
    for (CoarseReal y : ticks(y_min_, y_max_, y_ticks, log_y_))
    {
        double py = y_pos(y);
        fs << "<line x1='" << margin_left - 5 << "' y1='" << py << "' x2='" << margin_left << "' y2='" << py
           << "' stroke='gray'/>\n"
           << "<text x='" << margin_left - 8 << "' y='" << py + 4 << "' text-anchor='end' font-size='13'>"
           << (log_y_ ? detail::tick_label(y, true) : detail::axis_label(static_cast<double>(y), static_cast<double>(y_min_), static_cast<double>(y_max_), y_ticks, false)) << "</text>\n";
    }
    fs << "<text x='" << 18 << "' y='" << margin_top + H/2 << "' text-anchor='middle' font-size='15' transform='rotate(-90 18 "
       << margin_top + H/2 << ")'>" << y_label_ << "</text>\n";
    if (show_function_)
    {
        for (CoarseReal y : ticks(y_min_, y_max_, y_ticks, log_y_))
        {
            double py = function_top + (y_pos(y) - margin_top);
            fs << "<line x1='" << margin_left - 5 << "' y1='" << py << "' x2='" << margin_left << "' y2='" << py
               << "' stroke='gray'/>\n"
               << "<text x='" << margin_left - 8 << "' y='" << py + 4 << "' text-anchor='end' font-size='13'>"
               << (log_y_ ? detail::tick_label(y, true) : detail::axis_label(static_cast<double>(y), static_cast<double>(y_min_), static_cast<double>(y_max_), y_ticks, false)) << "</text>\n";
        }
        fs << "<text x='" << 18 << "' y='" << function_top + H/2 << "' text-anchor='middle' font-size='15' transform='rotate(-90 18 "
           << function_top + H/2 << ")'>" << y_label_ << "</text>\n";
    }

    // Colorbar, with the floor at the bottom and the clip at the top:
    int const bar_x = margin_left + plots_width + bar_gap;
    int const steps = 256;
    for (int s = 0; s < steps; ++s)
    {
        double t = 1 - (s + 0.5)/steps;
        fs << "<rect x='" << bar_x << "' y='" << margin_top + double(H)*s/steps << "' width='" << bar_width
           << "' height='" << double(H)/steps + 0.5 << "' fill='" << rgb_string(detail::to_8bit_rgb(color_map_(t))) << "'/>\n";
    }
    fs << "<rect x='" << bar_x << "' y='" << margin_top << "' width='" << bar_width << "' height='" << H
       << "' fill='none' stroke='gray'/>\n";
    int bar_label_width = 0;
    auto bar_tick = [&](double t, std::string const & label) {
        bar_label_width = (std::max)(bar_label_width, detail::label_width(label));
        double py = margin_top + (1 - t)*H;
        fs << "<line x1='" << bar_x + bar_width << "' y1='" << py << "' x2='" << bar_x + bar_width + 5 << "' y2='" << py
           << "' stroke='gray'/>\n"
           << "<text x='" << bar_x + bar_width + 8 << "' y='" << py + 4 << "' font-size='13'>" << label << "</text>\n";
    };
    auto number = [](double v) {
        std::ostringstream os;
        os << std::setprecision(3) << v;
        return os.str();
    };
    if (log_scale_)
    {
        double decades = std::log10(clip/floor);
        // Label every power of ten, or every step-th one when the scale spans many decades, so that the labels don't overlap:
        int step = (std::max)(1, static_cast<int>(std::ceil(decades/(std::max)(1, (std::min)(8, H/24)))));
        for (double e = step*std::ceil(std::log10(floor)/step); e < std::log10(clip); e += step)
        {
            double t = (e - std::log10(floor))/decades;
            // Keep clear of the end labels:
            if (t*H > 16 && (1 - t)*H > 16)
            {
                bar_tick(t, detail::power_label(std::pow(10.0, e)));
            }
        }
    }
    else
    {
        for (int i = 1; i < 4; ++i)
        {
            bar_tick(i/4.0, number(floor + (clip - floor)*i/4));
        }
    }
    bar_tick(0, log_scale_ ? detail::power_label(floor) : number(floor));
    // Say "at least" only when some value exceeds the clip, i.e. is hidden by the saturated scale:
    bool clipped = false;
    for (panel const & p : panels_)
    {
        for (std::size_t k = 0; k < p.state.size(); ++k)
        {
            if (p.state[k] == detail::ulps_class::ok && std::isfinite((p.*values)[k]) && magnitude((p.*values)[k]) > clip*(1 + 1e-9))
            {
                clipped = true;
            }
        }
    }
    {
        std::ostringstream top;
        top << std::setprecision(4) << clip;
        bar_tick(1, (clipped ? "&#8805; " : "") + (log_scale_ ? detail::power_label(clip) : top.str()));
    }
    // Legend for the pixels without a color on the scale:
    int legend_y = margin_top + H + 18;
    auto swatch = [&](std::array<std::uint8_t, 3> const & c, std::string const & text) {
        fs << "<rect x='" << bar_x << "' y='" << legend_y - 10 << "' width='12' height='12' fill='" << rgb_string(c) << "'/>\n"
           << "<text x='" << bar_x + 18 << "' y='" << legend_y << "' font-size='11'>" << text << "</text>\n";
        legend_y += 14;
    };
    if (any_failed) { swatch(failed_rgb, "NaN or domain error"); }
    if (any_spurious) { swatch(spurious_overflow_rgb, "spurious overflow"); }
    if (any_no_reference) { swatch(no_reference_rgb, "out of range, handled correctly"); }
    // The colorbar's label runs down its right side, clear of the title:
    int const label_x = (std::min)(bar_x + bar_width + 8 + bar_label_width + 16, total_width - 10);
    fs << "<text x='" << label_x << "' y='" << margin_top + H/2 << "' text-anchor='middle' font-size='13' transform='rotate(90 "
       << label_x << " " << margin_top + H/2 << ")'>" << quantity << "</text>\n";
    if (show_function_)
    {
        for (int st = 0; st < steps; ++st)
        {
            double t = 1 - (st + 0.5)/steps;
            auto c = detail::to_8bit_rgb(mode == fn_mode::linear_div || mode == fn_mode::symlog_div ? detail::dark_diverging(t) : viridis<double>(t));
            fs << "<rect x='" << bar_x << "' y='" << function_top + double(H)*st/steps << "' width='" << bar_width
               << "' height='" << double(H)/steps + 0.5 << "' fill='" << rgb_string(c) << "'/>\n";
        }
        fs << "<rect x='" << bar_x << "' y='" << function_top << "' width='" << bar_width << "' height='" << H
           << "' fill='none' stroke='gray'/>\n";
        int fn_label_width = 0;
        auto fn_tick = [&](double t, std::string const & label) {
            fn_label_width = (std::max)(fn_label_width, detail::label_width(label));
            double py = function_top + (1 - t)*H;
            fs << "<line x1='" << bar_x + bar_width << "' y1='" << py << "' x2='" << bar_x + bar_width + 5 << "' y2='" << py
               << "' stroke='gray'/>\n"
               << "<text x='" << bar_x + bar_width + 8 << "' y='" << py + 4 << "' font-size='13'>" << label << "</text>\n";
        };
        std::string fn_quantity = "f(x, y)";
        if (!function_label_.empty()) { fn_quantity = function_label_; }
        else if (!title_.empty()) { fn_quantity = title_; }
        // Enough digits that the ticks, a quarter of the colorbar apart, don't all read alike:
        int const fn_precision = detail::tick_precision(mode == fn_mode::linear_div ? -fsym : fmin, mode == fn_mode::linear_div ? fsym : fmax, 4);
        auto fn_number = [&](double v) { std::ostringstream os; os << std::setprecision(fn_precision) << v; return os.str(); };
        switch (mode)
        {
        case fn_mode::log_seq:
        {
            double lo = std::log10(fmin);
            double hi = std::log10(fmax);
            int step = (std::max)(1, static_cast<int>(std::ceil((hi - lo)/(std::max)(1, (std::min)(6, H/24)))));
            for (int e = static_cast<int>(std::ceil(lo)); e <= hi; e += step)
            {
                double t = fn_t(std::pow(10.0, e));
                if (t*H > 16 && (1 - t)*H > 16)
                {
                    fn_tick(t, detail::power_label(std::pow(10.0, e)));
                }
            }
            break;
        }
        case fn_mode::linear_seq:
            for (int i = 1; i < 4; ++i)
            {
                fn_tick(i/4.0, fn_number(fmin + (fmax - fmin)*i/4));
            }
            break;
        case fn_mode::linear_div:
            fn_tick(0.25, fn_number(-fsym/2));
            fn_tick(0.5, "0");
            fn_tick(0.75, fn_number(fsym/2));
            break;
        case fn_mode::symlog_div:
        {
            fn_tick(0.5, "0");
            double lo = std::log10(fscale);
            double hi = std::log10(fsym);
            int step = (std::max)(1, static_cast<int>(std::ceil((hi - lo)/(std::max)(1, (std::min)(3, H/48)))));
            for (int e = static_cast<int>(std::ceil(lo)); e <= hi; e += step)
            {
                double t = fn_t(std::pow(10.0, e));
                if ((t - 0.5)*H > 16 && (1 - t)*H > 16)
                {
                    fn_tick(t, detail::power_label(std::pow(10.0, e)));
                    fn_tick(1 - t, "-" + detail::power_label(std::pow(10.0, e)));
                }
            }
            break;
        }
        }
        // The ends of the colorbar are the extreme values, not clipped:
        if (mode == fn_mode::symlog_div)
        {
            fn_tick(0, "-" + detail::power_label(fsym));
            fn_tick(1, detail::power_label(fsym));
        }
        else
        {
            fn_tick(0, mode == fn_mode::log_seq ? detail::power_label(fmin) : fn_number(mode == fn_mode::linear_div ? -fsym : fmin));
            fn_tick(1, mode == fn_mode::log_seq ? detail::power_label(fmax) : fn_number(mode == fn_mode::linear_div ? fsym : fmax));
        }
        int const fn_label_x = (std::min)(bar_x + bar_width + 8 + fn_label_width + 16, total_width - 10);
        fs << "<text x='" << fn_label_x << "' y='" << function_top + H/2 << "' text-anchor='middle' font-size='13' transform='rotate(90 "
           << fn_label_x << " " << function_top + H/2 << ")'>" << fn_quantity << "</text>\n";
    }
    fs << "</svg>\n";
}

} BOOST_MATH_NAMESPACE_END
#endif // BOOST_MATH_TOOLS_ULPS_HEATMAP_HPP
