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

enum class heatmap_pixel : std::uint8_t
{
    ok,
    // The reference is undefined (NaN or threw) or outside CoarseReal's range:
    masked,
    // The reference is finite but the implementation returned NaN or inf, or threw:
    failed
};

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
        std::vector<detail::heatmap_pixel> state;
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
    p.state.assign(n, detail::heatmap_pixel::ok);
    detail::parallel_for(ys_.size(), threads_, [&](std::size_t j) {
        using std::abs;
        using std::isfinite;
        using std::nextafter;
        for (std::size_t i = 0; i < width; ++i)
        {
            std::size_t k = j*width + i;
            PreciseReal const & z = precise_[k];
            if (!isfinite(z) || abs(z) > (std::numeric_limits<CoarseReal>::max)())
            {
                p.state[k] = detail::heatmap_pixel::masked;
                continue;
            }
            CoarseReal w;
            try
            {
                w = g(xs_[i], ys_[j]);
            }
            catch (...)
            {
                p.state[k] = detail::heatmap_pixel::failed;
                continue;
            }
            if (!isfinite(w))
            {
                p.state[k] = detail::heatmap_pixel::failed;
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
        std::size_t no_value = 0;
        std::size_t over_half = 0;
        std::size_t over_envelope = 0;
        std::size_t worst_ulps = 0;
        std::size_t worst_ratio = 0;
        for (std::size_t k = 0; k < p.state.size(); ++k)
        {
            if (p.state[k] != detail::heatmap_pixel::ok) { ++no_value; continue; }
            if (std::abs(p.ulps[k]) > 0.5) { ++over_half; }
            if (std::abs(p.ulps_envelope[k]) > 1) { ++over_envelope; }
            if (!(std::abs(p.ulps[k]) <= std::abs(p.ulps[worst_ulps]))) { worst_ulps = k; }
            if (!(std::abs(p.ulps_envelope[k]) <= std::abs(p.ulps_envelope[worst_ratio]))) { worst_ratio = k; }
        }
        if (!p.label.empty())
        {
            os << p.label << ": ";
        }
        os << p.state.size() << " pixels: " << no_value << " without a value (NaN, inf, throw or no reference), "
           << over_half << " over 1/2 ulp, " << over_envelope << " over the envelope.\n";
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
    // Pixels without a value: the implementation returned NaN or inf or threw, or there is no reference.
    static constexpr std::array<std::uint8_t, 3> no_value_rgb = {127, 255, 0};

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
                if (p.state[k] == detail::heatmap_pixel::ok && std::isfinite((p.*values)[k]))
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
    int const margin_right = 150;
    int const W = static_cast<int>(width);
    int const H = static_cast<int>(height);
    int const panels = static_cast<int>(panels_.size());
    int const plots_width = panels*W + (panels - 1)*panel_gap;
    int const total_width = margin_left + plots_width + margin_right;
    int const total_height = margin_top + H + margin_bottom;
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

    std::ofstream fs(filename);
    fs << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n"
       << "<svg xmlns='http://www.w3.org/2000/svg' width='" << total_width << "' height='" << total_height << "'>\n"
       << "<style>svg { font-family: times; fill: " << font_color_ << "; }</style>\n"
       << "<rect width='100%' height='100%' fill='" << background_color_ << "'/>\n";
    if (!title_.empty())
    {
        fs << "<text x='" << margin_left + plots_width/2 << "' y='31' text-anchor='middle' font-size='18'>" << title_ << "</text>\n";
    }
    for (int q = 0; q < panels; ++q)
    {
        panel const & p = panels_[q];
        int const left = margin_left + q*(W + panel_gap);
        std::vector<unsigned char> rgb(3*width*height);
        for (std::size_t k = 0; k < width*height; ++k)
        {
            double v = (p.*values)[k];
            std::array<std::uint8_t, 3> c = (p.state[k] == detail::heatmap_pixel::ok && std::isfinite(v))
                ? detail::to_8bit_rgb(color_map_(scale(v))) : no_value_rgb;
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
        for (CoarseReal x : ticks(x_min_, x_max_, 8, log_x_))
        {
            double px = left + x_pos(x);
            fs << "<line x1='" << px << "' y1='" << margin_top + H << "' x2='" << px << "' y2='" << margin_top + H + 5
               << "' stroke='gray'/>\n"
               << "<text x='" << px << "' y='" << margin_top + H + 20 << "' text-anchor='middle' font-size='13'>"
               << detail::tick_label(x, log_x_) << "</text>\n";
        }
        fs << "<text x='" << left + W/2 << "' y='" << margin_top + H + 45 << "' text-anchor='middle' font-size='15'>"
           << x_label_ << "</text>\n";
    }
    // The panels share the y axis, which is labelled on the first:
    for (CoarseReal y : ticks(y_min_, y_max_, 6, log_y_))
    {
        double py = y_pos(y);
        fs << "<line x1='" << margin_left - 5 << "' y1='" << py << "' x2='" << margin_left << "' y2='" << py
           << "' stroke='gray'/>\n"
           << "<text x='" << margin_left - 8 << "' y='" << py + 4 << "' text-anchor='end' font-size='13'>"
           << detail::tick_label(y, log_y_) << "</text>\n";
    }
    fs << "<text x='" << 18 << "' y='" << margin_top + H/2 << "' text-anchor='middle' font-size='15' transform='rotate(-90 18 "
       << margin_top + H/2 << ")'>" << y_label_ << "</text>\n";

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
    auto bar_tick = [&](double t, std::string const & label) {
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
        // Label 2 and 5 times each power of ten too when the scale spans less than two decades:
        std::vector<double> mantissas = decades < 2 ? std::vector<double>{1, 2, 5} : std::vector<double>{1};
        for (double e = std::floor(std::log10(floor)); e < std::log10(clip); ++e)
        {
            for (double m : mantissas)
            {
                double v = m*std::pow(10.0, e);
                double t = std::log10(v/floor)/decades;
                // Keep clear of the end labels:
                if (t > 0.05 && t < 0.95)
                {
                    bar_tick(t, m == 1 ? detail::tick_label(v, true) : number(v));
                }
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
    bar_tick(0, number(floor));
    bar_tick(1, "&#8805; " + number(clip));
    // The colorbar's label runs down its right side, clear of the title:
    int const label_x = bar_x + bar_width + 90;
    fs << "<text x='" << label_x << "' y='" << margin_top + H/2 << "' text-anchor='middle' font-size='13' transform='rotate(90 "
       << label_x << " " << margin_top + H/2 << ")'>" << quantity << (log_scale_ ? " (log scale)" : "") << "</text>\n";
    fs << "</svg>\n";
}

} BOOST_MATH_NAMESPACE_END
#endif // BOOST_MATH_TOOLS_ULPS_HEATMAP_HPP
