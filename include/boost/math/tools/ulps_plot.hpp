//  (C) Copyright Nick Thompson 2020.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#ifndef BOOST_MATH_TOOLS_ULP_PLOT_HPP
#define BOOST_MATH_TOOLS_ULP_PLOT_HPP
#ifndef BOOST_MATH_BUILD_MODULE
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <vector>
#include <utility>
#include <fstream>
#include <string>
#include <list>
#include <random>
#include <limits>
#include <stdexcept>
#include <cstdint>
#include <cmath>
#include <sstream>
#include <functional>
#endif
#include <boost/math/tools/is_standalone.hpp>
#include <boost/math/tools/config.hpp>

#ifndef BOOST_MATH_STANDALONE
#include <boost/random/uniform_real_distribution.hpp>
#endif

// The data layer (points, markers, envelope, curves) is drawn into PNGs embedded in the SVG, while the axes and text
// stay vector. If lodepng.h is on the include path, lodepng encodes them (link lodepng.cpp); see png_encoder.
#if defined(__has_include)
#if __has_include("lodepng.h")
#include "lodepng.h"
#define BOOST_MATH_ULPS_HAS_LODEPNG
#define BOOST_MATH_ULPS_HEATMAP_HAS_LODEPNG
#endif
#endif

// Design of this function comes from:
// https://blogs.mathworks.com/cleve/2017/01/23/ulps-plots-reveal-math-function-accurary/

// The envelope is the maximum of 1/2 and half the condition number of function evaluation.
//
// Sampling: the abscissas are a uniform grid in the plotted coordinate (x, or log(x) on a logarithmic axis), plus random
// probes between the grid points. The condition number |x f'(x)/f(x)| at the grid points comes from finite differences of
// the reference values on the grid (see detail::grid_condition_numbers), so it costs no extra evaluations of f.

BOOST_MATH_NAMESPACE_BEGIN namespace tools {

namespace detail {

// Encodes 8-bit RGB pixels, row-major from the top row, as a PNG.
using png_encoder_type = std::function<std::vector<unsigned char>(std::vector<unsigned char> const &, std::size_t, std::size_t)>;

inline png_encoder_type default_png_encoder()
{
#ifdef BOOST_MATH_ULPS_HAS_LODEPNG
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

// ---- Condition numbers by finite differences on a grid ----

// How well a finite-difference derivative must agree with a lower-order one to be trusted, as a fraction of
// max(|f|/|x|, |f'|): that is, of the scale on which an error in f' changes the condition number |x f'/f|.
constexpr double fd_tolerance = 0.05;

// f'(t[j]) from the polynomial interpolating (t[k], f[k]), k < n: sum over k != j of L_k'(t[j]) (f[k] - f[j]).
template<class Real>
Real derivative_at_node(Real const * t, Real const * f, int n, int j)
{
    Real d = 0;
    for (int k = 0; k < n; ++k)
    {
        if (k == j)
        {
            continue;
        }
        Real num = 1;
        Real den = t[k] - t[j];
        for (int m = 0; m < n; ++m)
        {
            if (m != k && m != j)
            {
                num *= t[j] - t[m];
                den *= t[k] - t[m];
            }
        }
        d += (num/den)*(f[k] - f[j]);
    }
    return d;
}

// The Lagrange basis polynomials of the nodes t[0..n) at s.
template<class Real>
void lagrange_weights(Real const * t, int n, Real const & s, Real * w)
{
    for (int k = 0; k < n; ++k)
    {
        Real num = 1;
        Real den = 1;
        for (int m = 0; m < n; ++m)
        {
            if (m != k)
            {
                num *= s - t[m];
                den *= t[k] - t[m];
            }
        }
        w[k] = num/den;
    }
}

// The start of the window of four consecutive nodes of t[0..n) (sorted, ascending or descending) around s.
template<class Real>
std::size_t interpolation_window(Real const * t, std::size_t n, Real const & s)
{
    bool const ascending = t[0] < t[n - 1];
    std::size_t lo = 0;
    std::size_t hi = n;
    // The number of nodes on the near side of s:
    while (lo < hi)
    {
        std::size_t mid = (lo + hi)/2;
        if (ascending ? (t[mid] <= s) : (t[mid] >= s)) { lo = mid + 1; } else { hi = mid; }
    }
    std::size_t start = lo >= 2 ? lo - 2 : 0;
    return (std::min)(start, n - 4);
}

// The derivative of f at each node x[j] of a grid of n nodes (sorted, ascending or descending), from the values f[j*stride].
// Each comes from a five-point stencil: centered (fourth order) where it can be, shifted one-sided near the edges and next to
// values that are not finite. trusted[j*stride] says whether it agrees with the three-point (second order) estimate from the
// same values to within fd_tolerance*max(|f|/|x|, |f'|): if not, f is not smooth on the scale of the grid there (a pole,
// a kink, a narrow feature, or too coarse a grid), and the condition number should be computed some other way.
// d is NaN where f is not finite or no stencil fits.
template<class Real>
void grid_derivatives(Real const * x, Real const * f, std::size_t stride, std::size_t n, Real * d, unsigned char * trusted)
{
    using std::abs;
    using std::isfinite;
    static constexpr int offsets[5] = {-2, -1, -3, 0, -4};
    for (std::size_t j = 0; j < n; ++j)
    {
        d[j*stride] = std::numeric_limits<Real>::quiet_NaN();
        trusted[j*stride] = 0;
        if (!isfinite(f[j*stride]) || n < 5)
        {
            continue;
        }
        long start = -1;
        for (int o : offsets)
        {
            long s = static_cast<long>(j) + o;
            if (s < 0 || s + 4 >= static_cast<long>(n))
            {
                continue;
            }
            bool finite = true;
            for (long k = s; k < s + 5; ++k)
            {
                finite = finite && isfinite(f[static_cast<std::size_t>(k)*stride]);
            }
            if (finite)
            {
                start = s;
                break;
            }
        }
        if (start < 0)
        {
            continue;
        }
        Real ts[5];
        Real fs[5];
        for (int k = 0; k < 5; ++k)
        {
            ts[k] = x[start + k];
            fs[k] = f[static_cast<std::size_t>(start + k)*stride];
        }
        int const jj = static_cast<int>(static_cast<long>(j) - start);
        Real const d5 = derivative_at_node(ts, fs, 5, jj);
        int const s3 = (std::min)(2, (std::max)(0, jj - 1));
        Real const d3 = derivative_at_node(ts + s3, fs + s3, 3, jj - s3);
        d[j*stride] = d5;
        // |x| |d5 - d3| <= tol max(|f|, |x d5|), which holds trivially at x = 0, where the condition number is 0:
        Real const ax = abs(x[j]);
        trusted[j*stride] = isfinite(d5) && ax*abs(d5 - d3) <= Real(fd_tolerance)*(std::max)(Real(abs(fs[jj])), Real(ax*abs(d5)));
    }
}

// |x f'(x)/fx| at a single point, with fx = f(x), from a central difference with relative step cbrt(eps), or a one-sided one
// with relative step sqrt(eps) if a side is undefined. Two or three extra evaluations of f; used where the grid can't be trusted.
template<class F, class Real>
Real local_condition_number(F const & f, Real const & x, Real const & fx)
{
    using std::abs;
    using std::isfinite;
    using std::pow;
    using std::sqrt;
    if (x == 0)
    {
        return 0;
    }
    static Real const cbrt_eps = pow(std::numeric_limits<Real>::epsilon(), Real(1)/3);
    Real const sqrt_eps = sqrt(std::numeric_limits<Real>::epsilon());
#ifndef BOOST_MATH_NO_EXCEPTIONS
    try
    {
#endif
        Real xp = x + x*cbrt_eps;
        Real xm = x - x*cbrt_eps;
        Real slope = (f(xp) - f(xm))/(xp - xm);
        if (isfinite(slope))
        {
            return abs(x*slope/fx);
        }
#ifndef BOOST_MATH_NO_EXCEPTIONS
    }
    catch (...)
    {
    }
#endif
    for (int sign : {1, -1})
    {
#ifndef BOOST_MATH_NO_EXCEPTIONS
        try
        {
#endif
            Real xh = x + sign*x*sqrt_eps;
            Real slope = (f(xh) - fx)/(xh - x);
            if (isfinite(slope))
            {
                return abs(x*slope/fx);
            }
#ifndef BOOST_MATH_NO_EXCEPTIONS
        }
        catch (...)
        {
        }
#endif
    }
    return std::numeric_limits<Real>::quiet_NaN();
}

// Given the grid derivatives d (from grid_derivatives) of f on the grid x[0..n), estimates f'(s) at a point s between the
// nodes from the cubic through the four nearest grid derivatives. f' is smooth where f has a simple zero, so dividing by the
// exact f(s) = fs gives the narrow spike of |s f'/f| there exactly, which interpolating the condition number would smear out.
// The estimate is used only if the four nodes are trusted and the cubic through the grid values of f predicts fs well enough
// that the condition number is good to about fd_tolerance; otherwise f has a feature between the nodes, and this returns false.
template<class Real>
bool interpolate_derivative(Real const * x, Real const * f, Real const * d, unsigned char const * trusted, std::size_t n,
                            Real const & s, Real const & fs, Real & ds)
{
    using std::abs;
    using std::isfinite;
    if (n < 5)
    {
        return false;
    }
    std::size_t const start = interpolation_window(x, n, s);
    for (std::size_t k = start; k < start + 4; ++k)
    {
        if (!trusted[k])
        {
            return false;
        }
    }
    Real w[4];
    lagrange_weights(x + start, 4, s, w);
    Real fi = 0;
    ds = 0;
    for (std::size_t k = 0; k < 4; ++k)
    {
        fi += w[k]*f[start + k];
        ds += w[k]*d[start + k];
    }
    // An error e in the interpolated f means an error of about e/h in f', so |s| e/(h |f|) in the condition number:
    Real const h = abs(x[start + 2] - x[start + 1]);
    Real const as = abs(s);
    return isfinite(ds) && as*abs(fi - fs) <= Real(fd_tolerance)*h*(std::max)(Real(abs(fs)), Real(as*abs(ds)));
}

// The condition numbers |x f'(x)/f(x)| at the nodes of the grid x (sorted) with reference values fx = f(x), and at the probes
// px (anywhere in the grid's range) with values fpx, returned in pcond. Inf where f = 0, NaN where f is not finite.
// The grid's come from finite differences of fx (grid_derivatives), and the probes' from interpolate_derivative; where these
// can't be trusted, local_condition_number evaluates f at nearby points. If fallbacks isn't null, it receives the number of those.
template<class F, class Real>
std::vector<Real> grid_condition_numbers(F const & f, std::vector<Real> const & x, std::vector<Real> const & fx,
                                         std::vector<Real> const & px, std::vector<Real> const & fpx, std::vector<Real> & pcond,
                                         std::size_t * fallbacks = nullptr)
{
    using std::abs;
    using std::isfinite;
    std::size_t const n = x.size();
    std::size_t local = 0;
    std::vector<Real> d(n);
    std::vector<unsigned char> trusted(n);
    grid_derivatives(x.data(), fx.data(), 1, n, d.data(), trusted.data());
    std::vector<Real> cond(n, std::numeric_limits<Real>::quiet_NaN());
    for (std::size_t j = 0; j < n; ++j)
    {
        if (!isfinite(fx[j]))
        {
            continue;
        }
        if (trusted[j])
        {
            cond[j] = abs(x[j]*d[j]/fx[j]);
        }
        else
        {
            cond[j] = local_condition_number(f, x[j], fx[j]);
            ++local;
        }
    }
    pcond.assign(px.size(), std::numeric_limits<Real>::quiet_NaN());
    for (std::size_t j = 0; j < px.size(); ++j)
    {
        if (!isfinite(fpx[j]))
        {
            continue;
        }
        Real dj;
        if (interpolate_derivative(x.data(), fx.data(), d.data(), trusted.data(), n, px[j], fpx[j], dj))
        {
            pcond[j] = abs(px[j]*dj/fpx[j]);
        }
        else
        {
            pcond[j] = local_condition_number(f, px[j], fpx[j]);
            ++local;
        }
    }
    if (fallbacks)
    {
        *fallbacks = local;
    }
    return cond;
}

// ---- Drawing the data layer ----

// Clips the segment (x0, y0)-(x1, y1) to the rectangle [xmin, xmax] x [ymin, ymax] (Liang-Barsky); false if nothing is left.
inline bool clip_segment(double & x0, double & y0, double & x1, double & y1, double xmin, double ymin, double xmax, double ymax)
{
    if (!(std::isfinite(x0) && std::isfinite(y0) && std::isfinite(x1) && std::isfinite(y1)))
    {
        return false;
    }
    double t0 = 0;
    double t1 = 1;
    double const dx = x1 - x0;
    double const dy = y1 - y0;
    double const p[4] = {-dx, dx, -dy, dy};
    double const q[4] = {x0 - xmin, xmax - x0, y0 - ymin, ymax - y0};
    for (int i = 0; i < 4; ++i)
    {
        if (p[i] == 0)
        {
            if (q[i] < 0) { return false; }
            continue;
        }
        double r = q[i]/p[i];
        if (p[i] < 0) { t0 = (std::max)(t0, r); } else { t1 = (std::min)(t1, r); }
        if (t0 > t1) { return false; }
    }
    double const ox = x0;
    double const oy = y0;
    x0 = ox + t0*dx;
    y0 = oy + t0*dy;
    x1 = ox + t1*dx;
    y1 = oy + t1*dy;
    return true;
}

struct plot_point { double x; double y; };

// Draws the primitives of the data layer as SVG elements, clipped to a rectangle. Used when there is no PNG encoder.
class vector_canvas
{
public:
    vector_canvas(std::ostream & os, double x0, double y0, double width, double height, std::string const & id)
        : os_(os), x0_(x0), y0_(y0), x1_(x0 + width), y1_(y0 + height)
    {
        os_ << "<defs><clipPath id='" << id << "'><rect x='" << x0 << "' y='" << y0 << "' width='" << width << "' height='" << height
            << "'/></clipPath></defs>\n<g clip-path='url(#" << id << ")'>\n";
    }
    void disk(std::string const & color, double cx, double cy, double r)
    {
        os_ << "<circle cx='" << cx << "' cy='" << cy << "' r='" << r << "' fill='" << color << "'/>\n";
    }
    void segment(std::string const & color, double ax, double ay, double bx, double by, double width)
    {
        if (clip_segment(ax, ay, bx, by, x0_ - 1, y0_ - 1, x1_ + 1, y1_ + 1))
        {
            os_ << "<line x1='" << ax << "' y1='" << ay << "' x2='" << bx << "' y2='" << by << "' stroke='" << color
                << "' stroke-width='" << width << "'/>\n";
        }
    }
    void polygon(std::string const & color, std::vector<plot_point> const & p)
    {
        os_ << "<path d='M";
        for (std::size_t i = 0; i < p.size(); ++i)
        {
            os_ << (i ? " L" : "") << p[i].x << " " << p[i].y;
        }
        os_ << " Z' fill='" << color << "'/>\n";
    }
    void polyline(std::string const & color, std::vector<plot_point> const & p, double width)
    {
        // One path per run of segments that stay (after clipping) inside the rectangle:
        bool pen_down = false;
        for (std::size_t i = 0; i + 1 < p.size(); ++i)
        {
            double ax = p[i].x, ay = p[i].y, bx = p[i + 1].x, by = p[i + 1].y;
            if (!clip_segment(ax, ay, bx, by, x0_ - 1, y0_ - 1, x1_ + 1, y1_ + 1))
            {
                end_path(pen_down, color, width);
                continue;
            }
            if (!pen_down || ax != p[i].x || ay != p[i].y)
            {
                end_path(pen_down, color, width);
                os_ << "<path d='M" << ax << " " << ay;
                pen_down = true;
            }
            os_ << " L" << bx << " " << by;
            if (bx != p[i + 1].x || by != p[i + 1].y)
            {
                end_path(pen_down, color, width);
            }
        }
        end_path(pen_down, color, width);
    }
    void finish()
    {
        os_ << "</g>\n";
    }
private:
    void end_path(bool & pen_down, std::string const & color, double width)
    {
        if (pen_down)
        {
            os_ << "' stroke='" << color << "' stroke-width='" << width << "' fill='none'/>\n";
            pen_down = false;
        }
    }
    std::ostream & os_;
    double x0_, y0_, x1_, y1_;
};

// Rasterizes the primitives of the data layer, with 4x4 supersampling for antialiasing, into one coverage mask per color.
// finish() writes each mask as a PNG in the SVG and paints its color through it, so colors stay CSS strings and the
// background, gridlines and text show through. A pixel of the mask covers one unit of the SVG's user space.
class raster_canvas
{
public:
    raster_canvas(double x0, double y0, int width, int height)
        : x0_(x0), y0_(y0), width_(width), height_(height)
    {
    }
    void disk(std::string const & color, double cx, double cy, double r)
    {
        fill(layer(color), cx - r, cy - r, cx + r, cy + r, [&](double x, double y) {
            return (x - cx)*(x - cx) + (y - cy)*(y - cy) <= r*r;
        });
    }
    void segment(std::string const & color, double ax, double ay, double bx, double by, double width)
    {
        double const h = width/2;
        if (!clip_segment(ax, ay, bx, by, x0_ - h, y0_ - h, x0_ + width_ + h, y0_ + height_ + h))
        {
            return;
        }
        double const dx = bx - ax;
        double const dy = by - ay;
        double const len2 = dx*dx + dy*dy;
        fill(layer(color), (std::min)(ax, bx) - h, (std::min)(ay, by) - h, (std::max)(ax, bx) + h, (std::max)(ay, by) + h,
             [&](double x, double y) {
            double t = len2 > 0 ? ((x - ax)*dx + (y - ay)*dy)/len2 : 0;
            t = (std::min)(1.0, (std::max)(0.0, t));
            double ex = ax + t*dx - x;
            double ey = ay + t*dy - y;
            return ex*ex + ey*ey <= h*h;
        });
    }
    // A convex polygon:
    void polygon(std::string const & color, std::vector<plot_point> const & p)
    {
        double xmin = p[0].x, xmax = p[0].x, ymin = p[0].y, ymax = p[0].y;
        double area = 0;
        for (std::size_t i = 0; i < p.size(); ++i)
        {
            plot_point const & a = p[i];
            plot_point const & b = p[(i + 1) % p.size()];
            xmin = (std::min)(xmin, a.x); xmax = (std::max)(xmax, a.x);
            ymin = (std::min)(ymin, a.y); ymax = (std::max)(ymax, a.y);
            area += a.x*b.y - b.x*a.y;
        }
        fill(layer(color), xmin, ymin, xmax, ymax, [&](double x, double y) {
            for (std::size_t i = 0; i < p.size(); ++i)
            {
                plot_point const & a = p[i];
                plot_point const & b = p[(i + 1) % p.size()];
                if (((b.x - a.x)*(y - a.y) - (b.y - a.y)*(x - a.x))*area < 0)
                {
                    return false;
                }
            }
            return true;
        });
    }
    void polyline(std::string const & color, std::vector<plot_point> const & p, double width)
    {
        for (std::size_t i = 0; i + 1 < p.size(); ++i)
        {
            segment(color, p[i].x, p[i].y, p[i + 1].x, p[i + 1].y, width);
        }
    }
    // Writes each color's mask, in the order the colors were first used. id must be unique in the document.
    void finish(std::ostream & os, png_encoder_type const & encoder, std::string const & id) const
    {
        for (std::size_t l = 0; l < layers_.size(); ++l)
        {
            std::vector<unsigned char> rgb(3*layers_[l].coverage.size());
            for (std::size_t k = 0; k < layers_[l].coverage.size(); ++k)
            {
                // Coverage counts subsamples, 0 to 16:
                unsigned char v = static_cast<unsigned char>((255*layers_[l].coverage[k] + 8)/16);
                rgb[3*k] = rgb[3*k + 1] = rgb[3*k + 2] = v;
            }
            std::string mask_id = id + "-" + std::to_string(l);
            os << "<mask id='" << mask_id << "' maskUnits='userSpaceOnUse' x='" << x0_ << "' y='" << y0_ << "' width='" << width_
               << "' height='" << height_ << "'><image x='" << x0_ << "' y='" << y0_ << "' width='" << width_ << "' height='" << height_
               << "' href='data:image/png;base64," << base64(encoder(rgb, static_cast<std::size_t>(width_), static_cast<std::size_t>(height_)))
               << "'/></mask>\n"
               << "<rect x='" << x0_ << "' y='" << y0_ << "' width='" << width_ << "' height='" << height_ << "' fill='"
               << layers_[l].color << "' mask='url(#" << mask_id << ")'/>\n";
        }
    }
private:
    struct color_layer
    {
        std::string color;
        std::vector<std::uint8_t> coverage;
    };
    color_layer & layer(std::string const & color)
    {
        for (color_layer & l : layers_)
        {
            if (l.color == color)
            {
                return l;
            }
        }
        layers_.push_back(color_layer{color, std::vector<std::uint8_t>(static_cast<std::size_t>(width_)*static_cast<std::size_t>(height_), 0)});
        return layers_.back();
    }
    // Sets the coverage of each pixel meeting the box [xmin, xmax] x [ymin, ymax] (in user space) to the larger of its
    // coverage and the number of its 4x4 subsamples inside the shape. The larger, not the sum, so that overlapping
    // shapes of one color (the segments of a polyline at their joints) don't darken.
    template<class Inside>
    void fill(color_layer & l, double xmin, double ymin, double xmax, double ymax, Inside inside)
    {
        if (!(xmin <= xmax && ymin <= ymax))
        {
            return;
        }
        int const i0 = (std::max)(0, static_cast<int>(std::floor((std::max)(xmin - x0_, -1.0))));
        int const i1 = (std::min)(width_ - 1, static_cast<int>(std::floor((std::min)(xmax - x0_, double(width_)))));
        int const j0 = (std::max)(0, static_cast<int>(std::floor((std::max)(ymin - y0_, -1.0))));
        int const j1 = (std::min)(height_ - 1, static_cast<int>(std::floor((std::min)(ymax - y0_, double(height_)))));
        for (int j = j0; j <= j1; ++j)
        {
            for (int i = i0; i <= i1; ++i)
            {
                std::uint8_t count = 0;
                for (int sj = 0; sj < 4; ++sj)
                {
                    for (int si = 0; si < 4; ++si)
                    {
                        count += inside(x0_ + i + (si + 0.5)/4, y0_ + j + (sj + 0.5)/4) ? 1 : 0;
                    }
                }
                std::uint8_t & c = l.coverage[static_cast<std::size_t>(j)*static_cast<std::size_t>(width_) + static_cast<std::size_t>(i)];
                c = (std::max)(c, count);
            }
        }
    }
    double x0_;
    double y0_;
    int width_;
    int height_;
    std::vector<color_layer> layers_;
};

// A short name for a document's SVG ids, from its content, so that several plots inlined in one page don't collide.
inline std::string svg_id(std::string const & prefix, std::string const & content)
{
    std::uint64_t h = 1469598103934665603ULL;
    for (char c : content)
    {
        h = (h ^ static_cast<unsigned char>(c))*1099511628211ULL;
    }
    std::ostringstream os;
    os << prefix << std::hex << (h & 0xffffffffffULL);
    return os.str();
}
// Significant digits that tell apart tick labels spaced (hi - lo)/intervals apart, and at least 4:
// a function that varies only in its sixth digit over the plot must not get identical labels.
inline int tick_precision(double lo, double hi, int intervals)
{
    double spacing = std::abs(hi - lo)/(std::max)(1, intervals);
    double magnitude = (std::max)(std::abs(lo), std::abs(hi));
    if (!(spacing > 0) || !(magnitude > 0))
    {
        return 4;
    }
    int digits = static_cast<int>(std::ceil(std::log10(magnitude/spacing))) + 1;
    return (std::min)(17, (std::max)(4, digits));
}

// Labels a tick on an axis spanning [lo, hi] with the given number of intervals between ticks. On a linear axis
// whose values are tiny or huge, use compact scientific notation (2e-4, 1.5e6) rather than a string of zeros.
inline std::string axis_label(double v, double lo, double hi, int intervals, bool log_axis)
{
    std::ostringstream os;
    double magnitude = (std::max)(std::abs(lo), std::abs(hi));
    int precision = tick_precision(lo, hi, intervals);
    if (log_axis || v == 0 || !(magnitude < 1e-2 || magnitude >= 1e5))
    {
        os << std::setprecision(log_axis ? 4 : precision) << v;
        return os.str();
    }
    os << std::scientific << std::setprecision(precision - 1) << v;
    std::string s = os.str();
    std::string::size_type e = s.find('e');
    std::string mantissa = s.substr(0, e);
    if (mantissa.find('.') != std::string::npos)
    {
        mantissa.erase(mantissa.find_last_not_of('0') + 1);
        if (mantissa.back() == '.')
        {
            mantissa.pop_back();
        }
    }
    return mantissa + "e" + std::to_string(std::stoi(s.substr(e + 1)));
}


// How an implementation's result relates to the true value, which is used to decide whether the ulps are meaningful.
enum class ulps_class : std::uint8_t
{
    // Both are in range: the ulps are meaningful.
    ok,
    // The true value is NaN or infinite (or its evaluation threw): f is undefined there. Not a defect.
    undefined,
    // |true value| > max: the implementation threw std::overflow_error or returned an infinity of the right sign. Correct.
    overflow_handled,
    // The true value is nonzero but below the smallest normal: the implementation returned zero or a subnormal of the right sign. Correct.
    underflow_handled,
    // The true value is a finite double, but the implementation returned NaN or threw something other than an overflow_error. A bug.
    failed,
    // The true value is a finite double, but the implementation threw std::overflow_error or returned an infinity. A bug.
    spurious_overflow,
    // The true value is out of range, but the implementation returned an ordinary number, NaN, the wrong sign, or threw something else.
    // Counted as the maximal error.
    mishandled
};

// Evaluates call() (which returns the implementation's value) and classifies the outcome against the true value z.
// When the result is ulps_class::ok, w holds the implementation's value.
template<class CoarseReal, class PreciseReal, class Call>
ulps_class classify_result(PreciseReal const & z, Call call, CoarseReal & w, std::string * message = nullptr)
{
    using std::abs;
    using std::isfinite;
    using std::isinf;
    using std::isnan;
    if (!isfinite(z))
    {
        return ulps_class::undefined;
    }
    bool const over = abs(z) > static_cast<PreciseReal>((std::numeric_limits<CoarseReal>::max)());
    bool const under = z != 0 && abs(z) < static_cast<PreciseReal>((std::numeric_limits<CoarseReal>::min)());
    try
    {
        w = static_cast<CoarseReal>(call());
    }
    catch (std::overflow_error const & e)
    {
        if (message) { *message = std::string("std::overflow_error: ") + e.what(); }
        return over ? ulps_class::overflow_handled : ulps_class::spurious_overflow;
    }
    catch (std::underflow_error const & e)
    {
        if (message) { *message = std::string("std::underflow_error: ") + e.what(); }
        return under ? ulps_class::underflow_handled : ulps_class::failed;
    }
    catch (std::exception const & e)
    {
        if (message) { *message = e.what(); }
        return (over || under) ? ulps_class::mishandled : ulps_class::failed;
    }
    catch (...)
    {
        if (message) { *message = "unknown exception"; }
        return (over || under) ? ulps_class::mishandled : ulps_class::failed;
    }
    if (isnan(w))
    {
        if (message) { *message = "returned NaN"; }
        return (over || under) ? ulps_class::mishandled : ulps_class::failed;
    }
    if (isinf(w))
    {
        if (message) { *message = "returned an infinity"; }
        if (over)
        {
            return ((w > 0) == (z > 0)) ? ulps_class::overflow_handled : ulps_class::mishandled;
        }
        return ulps_class::spurious_overflow;
    }
    if (over)
    {
        return ulps_class::mishandled;
    }
    if (under)
    {
        bool const right_sign = (w == 0) || ((w > 0) == (z > 0));
        bool const tiny = abs(w) < (std::numeric_limits<CoarseReal>::min)();
        return (right_sign && tiny) ? ulps_class::underflow_handled : ulps_class::mishandled;
    }
    return ulps_class::ok;
}

// Abscissas of the vertical gridlines: evenly spaced, or at whole powers of ten on a logarithmic axis.
template<class CoarseReal>
std::vector<CoarseReal> vertical_gridlines(CoarseReal min_x, CoarseReal max_x, int vertical_lines, bool log_x)
{
    using std::ceil;
    using std::floor;
    using std::log10;
    using std::pow;
    std::vector<CoarseReal> xs;
    if (!log_x)
    {
        for (int i = 1; i <= vertical_lines; ++i)
        {
            xs.push_back(min_x + ((max_x - min_x)*i)/vertical_lines);
        }
        return xs;
    }
    const int lo = static_cast<int>(ceil(log10(min_x)));
    const int hi = static_cast<int>(floor(log10(max_x)));
    const int step = (std::max)(1, static_cast<int>(ceil(double(hi - lo)/(std::max)(1, vertical_lines))));
    for (int e = hi; e >= lo; e -= step)
    {
        xs.push_back(static_cast<CoarseReal>(pow(10.0L, e)));
    }
    std::reverse(xs.begin(), xs.end());
    return xs;
}

template<class F1, class F2, class CoarseReal, class PreciseReal>
void write_gridlines(std::ostream& fs, int horizontal_lines, int vertical_lines,
                     F1 x_scale, F2 y_scale, CoarseReal min_x, CoarseReal max_x, PreciseReal min_y, PreciseReal max_y,
                     int graph_width, int graph_height, int margin_left, std::string const & font_color, bool log_x = false)
{
  // Make a grid:
  for (int i = 1; i <= horizontal_lines; ++i) {
      PreciseReal y_cord_dataspace = min_y +  ((max_y - min_y)*i)/horizontal_lines;
      auto y = y_scale(y_cord_dataspace);
      fs << "<line x1='0' y1='" << y << "' x2='" << graph_width
         << "' y2='" << y
         << "' stroke='gray' stroke-width='1' opacity='0.5' stroke-dasharray='4' />\n";

      fs << "<text x='" <<  -margin_left/4 + 5 << "' y='" << y - 3
         << "' font-family='times' font-size='10' fill='" << font_color << "' transform='rotate(-90 "
         << -margin_left/4 + 8 << " " << y + 5 << ")'>"
         << std::setprecision(4) << y_cord_dataspace << "</text>\n";
   }

    for (CoarseReal x_cord_dataspace : vertical_gridlines(min_x, max_x, vertical_lines, log_x)) {
        CoarseReal x = x_scale(x_cord_dataspace);
        fs << "<line x1='" << x << "' y1='0' x2='" << x
           << "' y2='" << graph_height
           << "' stroke='gray' stroke-width='1' opacity='0.5' stroke-dasharray='4' />\n";

        fs << "<text x='" <<  x - 10  << "' y='" << graph_height + 10
           << "' font-family='times' font-size='10' fill='" << font_color << "'>"
           << axis_label(static_cast<double>(x_cord_dataspace), static_cast<double>(min_x), static_cast<double>(max_x), vertical_lines, log_x) << "</text>\n";
    }
}
}

template<class F, typename PreciseReal, typename CoarseReal>
class ulps_plot {
public:
    // Of the samples, 4/5 are a uniform grid in the plotted coordinate and 1/5 random probes, uniform in the same coordinate
    // (seeded with random_seed, or randomly if it is -1). With perturb_abscissas, f is evaluated at the sample and the
    // implementation at the nearest CoarseReal; otherwise both at the CoarseReal.
    // With log_abscissas, the abscissas are sampled uniformly in log(x), and the x-axis is logarithmic; requires 0 < a.
    ulps_plot(F hi_acc_impl, CoarseReal a, CoarseReal b,
             size_t samples = 1000, bool perturb_abscissas = false, int random_seed = -1, bool log_abscissas = false);

    ulps_plot& clip(PreciseReal clip);

    ulps_plot& width(int width);

    ulps_plot& envelope_color(std::string const & color);

    ulps_plot& title(std::string const & title);

    ulps_plot& background_color(std::string const & background_color);

    ulps_plot& font_color(std::string const & font_color);

    ulps_plot& crop_color(std::string const & color);

    ulps_plot& nan_color(std::string const & color);

    ulps_plot& ulp_envelope(bool write_ulp);

    template<class G>
    ulps_plot& add_fn(G g, std::string const & color = "steelblue");

    ulps_plot& horizontal_lines(int horizontal_lines);

    ulps_plot& vertical_lines(int vertical_lines);

    // Draws the function itself (the high-accuracy ordinates) in a second plot below the ulps plot,
    // with the same width and x axis. Off by default.
    ulps_plot& show_function(bool show_function);

    // With show_function, asks for a logarithmic y axis, which is useful when f spans many decades. It is used
    // only if every finite value of f is positive; otherwise (f changes sign or is zero somewhere) the axis is linear.
    ulps_plot& log_function(bool log_function);

    ulps_plot& function_color(std::string const & color);

    // Encodes 8-bit RGB pixels, row-major from the top row, as a PNG; lodepng by default, if lodepng.h is found.
    // The points, markers and curves are drawn into PNGs embedded in the SVG, and the axes and text are vector.
    // With no encoder (nullptr, or no lodepng.h), everything is drawn as SVG elements, as before.
    using png_encoder_type = detail::png_encoder_type;
    ulps_plot& png_encoder(png_encoder_type encoder);

    void write(std::string const & filename) const;

    // For each function added: counts of samples by outcome (see detail::ulps_class), the worst ulps and the worst
    // ulps/envelope with the abscissas where they occur, and the first message of a failure.
    std::string summary() const;

    friend std::ostream& operator<<(std::ostream& fs, ulps_plot const & plot)
    {
        using std::abs;
        using std::floor;
        using std::isnan;
        using std::isfinite;
        if (plot.ulp_list_.size() == 0)
        {
            throw std::domain_error("No functions added for comparison.");
        }
        if (plot.width_ <= 1)
        {
            throw std::domain_error("Width = " + std::to_string(plot.width_) + ", which is too small.");
        }

        // The y range comes from the ulps alone (and clip). The condition number envelope never sets it: where the envelope
        // leaves the plot (e.g. near poles) it is clipped to the plot area.
        PreciseReal worst_ulp_distance = 0;
        PreciseReal min_y = (std::numeric_limits<PreciseReal>::max)();
        PreciseReal max_y = std::numeric_limits<PreciseReal>::lowest();
        for (auto const & ulp_vec : plot.ulp_list_)
        {
            for (auto const & ulp : ulp_vec)
            {
                if (static_cast<PreciseReal>(abs(ulp)) > worst_ulp_distance)
                {
                    worst_ulp_distance = static_cast<PreciseReal>(abs(ulp));
                }
                if (static_cast<PreciseReal>(ulp) < min_y)
                {
                    min_y = static_cast<PreciseReal>(ulp);
                }
                if (static_cast<PreciseReal>(ulp) > max_y)
                {
                    max_y = static_cast<PreciseReal>(ulp);
                }
            }
        }

        // half-ulp accuracy is the best that can be expected; sometimes we can get less, but barely less.
        // then the axes don't show up; painful!
        if (max_y < 0.5) {
            max_y = 0.5;
        }
        if (min_y > -0.5) {
            min_y = -0.5;
        }

        if (plot.clip_ > 0)
        {
            // With a clip, the range is symmetric, so that the ±envelope and the ticks are too: ±min(clip, worst |ulps|).
            PreciseReal extent = (std::max)(max_y, PreciseReal(-min_y));
            if (extent > plot.clip_)
            {
                extent = plot.clip_;
            }
            max_y = extent;
            min_y = -extent;
        }

        int height = static_cast<int>(floor(static_cast<double>(plot.width_)/1.61803));
        int margin_top = 40;
        int margin_left = 25;
        if (plot.title_.size() == 0)
        {
            margin_top = 10;
            margin_left = 15;
        }
        int margin_bottom = 20;
        int margin_right = 20;
        int graph_height = height - margin_bottom - margin_top;
        int graph_width = plot.width_ - margin_left - margin_right;

        // Maps [a,b] to [0, graph_width]
        auto x_scale = [&](CoarseReal x)->CoarseReal
        {
            using std::log;
            if (plot.log_abscissas_)
            {
                return (log(x/plot.a_)/log(plot.b_/plot.a_))*static_cast<CoarseReal>(graph_width);
            }
            return ((x-plot.a_)/(plot.b_ - plot.a_))*static_cast<CoarseReal>(graph_width);
        };

        auto y_scale = [&](PreciseReal y)->PreciseReal
        {
            return ((max_y - y)/(max_y - min_y) )*static_cast<PreciseReal>(graph_height);
        };

        // The function plot goes below the ulps plot, separated by room for the ulps plot's x tick labels:
        int const function_gap = 40;
        int const function_height = plot.show_function_ ? graph_height/2 : 0;
        int const total_height = plot.show_function_ ? height + function_gap + function_height + margin_bottom : height;
        fs << "<?xml version=\"1.0\" encoding='UTF-8' ?>\n"
           << "<svg xmlns='http://www.w3.org/2000/svg' width='"
           << plot.width_ << "' height='"
           << total_height << "'>\n"
           << "<style>\nsvg { background-color:" << plot.background_color_ << "; }\n"
           << "</style>\n";
        if (plot.title_.size() > 0)
        {
            fs << "<text x='" << floor(plot.width_/2)
               << "' y='" << floor(margin_top/2)
               << "' font-family='Palatino' font-size='25' fill='"
               << plot.font_color_  << "'  alignment-baseline='middle' text-anchor='middle'>"
               << plot.title_
               << "</text>\n";
        }

        // Construct SVG group to simplify the calculations slightly:
        fs << "<g transform='translate(" << margin_left << ", " << margin_top << ")'>\n";
            // y-axis:
        fs  << "<line x1='0' y1='0' x2='0' y2='" << graph_height
            << "' stroke='gray' stroke-width='1'/>\n";
        PreciseReal x_axis_loc = y_scale(static_cast<PreciseReal>(0));
        fs << "<line x1='0' y1='" << x_axis_loc
            << "' x2='" << graph_width << "' y2='" << x_axis_loc
            << "' stroke='gray' stroke-width='1'/>\n";

        if (worst_ulp_distance > 3)
        {
            detail::write_gridlines(fs, plot.horizontal_lines_, plot.vertical_lines_, x_scale, y_scale, plot.a_, plot.b_,
                                    min_y, max_y, graph_width, graph_height, margin_left, plot.font_color_, plot.log_abscissas_);
        }
        else
        {
            std::vector<double> ys{-3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0};
            for (double & i : ys)
            {
                if (min_y <= i && i <= max_y)
                {
                    PreciseReal y_cord_dataspace = i;
                    PreciseReal y = y_scale(y_cord_dataspace);
                    fs << "<line x1='0' y1='" << y << "' x2='" << graph_width
                       << "' y2='" << y
                       << "' stroke='gray' stroke-width='1' opacity='0.5' stroke-dasharray='4' />\n";

                    fs << "<text x='" <<  -margin_left/2 << "' y='" << y - 3
                       << "' font-family='times' font-size='10' fill='" << plot.font_color_ << "' transform='rotate(-90 "
                       << -margin_left/2 + 7 << " " << y << ")'>"
                       <<  std::setprecision(4) << y_cord_dataspace << "</text>\n";
                }
            }
            for (CoarseReal x_cord_dataspace : detail::vertical_gridlines(plot.a_, plot.b_, plot.vertical_lines_, plot.log_abscissas_))
            {
                CoarseReal x = x_scale(x_cord_dataspace);
                fs << "<line x1='" << x << "' y1='0' x2='" << x
                   << "' y2='" << graph_height
                   << "' stroke='gray' stroke-width='1' opacity='0.5' stroke-dasharray='4' />\n";

                fs << "<text x='" <<  x - 10  << "' y='" << graph_height + 10
                   << "' font-family='times' font-size='10' fill='" << plot.font_color_ << "'>"
                   << detail::axis_label(static_cast<double>(x_cord_dataspace), static_cast<double>(plot.a_), static_cast<double>(plot.b_), plot.vertical_lines_, plot.log_abscissas_) << "</text>\n";
            }
        }

        // Everything data-dependent (points, markers, the envelope) is clipped to the plot area, with room for the markers on its edges,
        // so that nothing runs into the neighboring panel. With a PNG encoder it is rasterized, and drawn as SVG elements otherwise.
        //
        // Results that are not comparisons are marked on the edges, and counted in a legend:
        //   green triangle (top):   NaN or a domain error, though the true value is a finite double (a bug);
        //   cyan diamond (top):     spurious overflow: std::overflow_error or an infinity, though the true value is a finite double (a bug);
        //   yellow circle (top):    out of range, but mishandled (counted as the maximal error);
        //   gray tick (bottom):     correctly handled: the true value is undefined, overflows (and std::overflow_error or an infinity
        //                           was returned), or underflows (and zero or a subnormal was returned).
        std::string const nan_color = plot.nan_color_.empty() ? "#7fff00" : plot.nan_color_;
        auto draw_data = [&](auto & canvas) {
            size_t color_idx = 0;
            for (auto const & ulp : plot.ulp_list_)
            {
                std::string const & color = plot.colors_[color_idx++];
                // Cropped points are crosses, in the color of their function unless one was set,
                // so that they stay distinguishable when several functions are plotted:
                std::string const & crop_color = plot.crop_color_set_ ? plot.crop_color_ : color;
                for (size_t j = 0; j < ulp.size(); ++j)
                {
                    if (isnan(ulp[j]))
                    {
                        // Not a comparison: see the markers below.
                        continue;
                    }
                    double x = static_cast<double>(x_scale(plot.coarse_abscissas_[j]));
                    if (plot.clip_ > 0 && static_cast<PreciseReal>(abs(ulp[j])) > plot.clip_)
                    {
                        if (crop_color.empty())
                        {
                            continue;
                        }
                        double y = static_cast<double>(y_scale(static_cast<PreciseReal>(ulp[j] < 0 ? -plot.clip_ : plot.clip_)));
                        canvas.segment(crop_color, x - 2.5, y - 2.5, x + 2.5, y + 2.5, 1);
                        canvas.segment(crop_color, x + 2.5, y - 2.5, x - 2.5, y + 2.5, 1);
                    }
                    else
                    {
                        canvas.disk(color, x, static_cast<double>(y_scale(static_cast<PreciseReal>(ulp[j]))), 1);
                    }
                }
            }
            for (auto const & classes : plot.class_list_)
            {
                for (size_t j = 0; j < classes.size(); ++j)
                {
                    double x = static_cast<double>(x_scale(plot.coarse_abscissas_[j]));
                    switch (classes[j])
                    {
                    case detail::ulps_class::ok:
                        break;
                    case detail::ulps_class::failed:
                        canvas.polygon(nan_color, {{x - 3, 1}, {x + 3, 1}, {x, 8}});
                        break;
                    case detail::ulps_class::spurious_overflow:
                        canvas.polygon("#00e5ff", {{x, 1}, {x + 4, 5}, {x, 9}, {x - 4, 5}});
                        break;
                    case detail::ulps_class::mishandled:
                        canvas.disk("black", x, 5, 4);
                        canvas.disk("#ffe600", x, 5, 3.5);
                        break;
                    default:
                        canvas.segment("#888888", x, graph_height - 7, x, graph_height - 1, 1);
                        break;
                    }
                }
            }
            if (plot.ulp_envelope_)
            {
                // ±cond, broken where it is undefined. Where it leaves the plot (near zeros and poles), it is clipped.
                for (int sign : {1, -1})
                {
                    std::vector<detail::plot_point> run;
                    for (size_t j = 0; j <= plot.cond_.size(); ++j)
                    {
                        if (j == plot.cond_.size() || isnan(plot.cond_[j]))
                        {
                            canvas.polyline(plot.envelope_color_, run, 1);
                            run.clear();
                            continue;
                        }
                        run.push_back({static_cast<double>(x_scale(plot.coarse_abscissas_[j])), static_cast<double>(y_scale(sign*plot.cond_[j]))});
                    }
                }
            }
        };
        std::string const id = detail::svg_id("ulps-", plot.svg_key());
        if (plot.png_encoder_)
        {
            detail::raster_canvas canvas(-4, -4, graph_width + 8, graph_height + 8);
            draw_data(canvas);
            canvas.finish(fs, plot.png_encoder_, id);
        }
        else
        {
            detail::vector_canvas canvas(fs, -4, -4, graph_width + 8, graph_height + 8, id);
            draw_data(canvas);
            canvas.finish();
        }
        {
            size_t counts[7] = {0, 0, 0, 0, 0, 0, 0};
            for (auto const & classes : plot.class_list_)
            {
                for (detail::ulps_class c : classes)
                {
                    ++counts[static_cast<size_t>(c)];
                }
            }
            size_t gray = counts[1] + counts[2] + counts[3];
            struct entry { std::string text; size_t n; int kind; };
            std::vector<entry> legend = {
                {"NaN or domain error", counts[4], 4},
                {"spurious overflow", counts[5], 5},
                {"out of range, mishandled", counts[6], 6},
                {"out of range or undefined, handled correctly", gray, 1}};
            int row = 0;
            for (auto const & e : legend)
            {
                if (e.n == 0)
                {
                    continue;
                }
                int ly = 16 + 15*row++;
                int lx = graph_width - 12;
                std::string fill = e.kind == 4 ? nan_color : e.kind == 5 ? "#00e5ff" : e.kind == 6 ? "#ffe600" : "#888888";
                fs << "<rect x='" << lx << "' y='" << ly - 9 << "' width='9' height='9' fill='" << fill << "'/>\n";
                fs << "<text x='" << lx - 5 << "' y='" << ly - 1 << "' font-family='times' font-size='12' text-anchor='end' fill='"
                   << plot.font_color_ << "'>" << e.text << " (" << e.n << ")</text>\n";
            }
        }
        fs << "</g>\n";
        if (plot.show_function_)
        {
            // The function, from the high-accuracy ordinates, on the same x axis as the ulps plot.
            // f itself is plotted, never |f|. A logarithmic axis needs every finite f to be positive.
            using std::log10;
            using std::pow;
            bool all_positive = true;
            for (auto const & y : plot.precise_ordinates_)
            {
                if (isfinite(y) && !(y > 0))
                {
                    all_positive = false;
                }
            }
            bool const log_f = plot.log_function_ && all_positive;
            auto transform = [&](PreciseReal y)->PreciseReal {
                return log_f ? log10(y) : y;
            };
            auto usable = [&](PreciseReal y)->bool {
                return isfinite(y);
            };
            PreciseReal lo = (std::numeric_limits<PreciseReal>::max)();
            PreciseReal hi = std::numeric_limits<PreciseReal>::lowest();
            for (auto const & y : plot.precise_ordinates_)
            {
                if (usable(y))
                {
                    PreciseReal t = transform(y);
                    lo = (std::min)(lo, t);
                    hi = (std::max)(hi, t);
                }
            }
            if (!(lo <= hi))
            {
                lo = 0;
                hi = 1;
            }
            if (log_f)
            {
                // Whole decades at the ends of the axis:
                lo = floor(lo);
                hi = -floor(-hi);
            }
            if (!(hi > lo))
            {
                // A constant function:
                lo -= 1;
                hi += 1;
            }
            auto f_scale = [&](PreciseReal t)->PreciseReal {
                return ((hi - t)/(hi - lo))*static_cast<PreciseReal>(function_height);
            };
            fs << "<g transform='translate(" << margin_left << ", " << margin_top + graph_height + function_gap << ")'>\n";
            fs << "<line x1='0' y1='0' x2='0' y2='" << function_height << "' stroke='gray' stroke-width='1'/>\n";
            fs << "<line x1='0' y1='" << function_height << "' x2='" << graph_width << "' y2='" << function_height
               << "' stroke='gray' stroke-width='1'/>\n";
            // y gridlines and labels:
            if (log_f)
            {
                int decades = static_cast<int>(hi - lo);
                int step = (std::max)(1, static_cast<int>(std::ceil(double(decades)/plot.horizontal_lines_)));
                for (int e = static_cast<int>(hi); e >= static_cast<int>(lo); e -= step)
                {
                    PreciseReal y = f_scale(static_cast<PreciseReal>(e));
                    fs << "<line x1='0' y1='" << y << "' x2='" << graph_width
                       << "' y2='" << y << "' stroke='gray' stroke-width='1' opacity='0.5' stroke-dasharray='4' />\n";
                    fs << "<text x='" << -margin_left/2 << "' y='" << y - 3
                       << "' font-family='times' font-size='10' fill='" << plot.font_color_ << "' transform='rotate(-90 "
                       << -margin_left/2 + 7 << " " << y << ")'>1e" << e << "</text>\n";
                }
            }
            else
            {
                int lines = (std::max)(2, plot.horizontal_lines_/2);
                int precision = detail::tick_precision(static_cast<double>(lo), static_cast<double>(hi), lines);
                for (int i = 0; i <= lines; ++i)
                {
                    PreciseReal v = lo + ((hi - lo)*i)/lines;
                    PreciseReal y = f_scale(v);
                    fs << "<line x1='0' y1='" << y << "' x2='" << graph_width
                       << "' y2='" << y << "' stroke='gray' stroke-width='1' opacity='0.5' stroke-dasharray='4' />\n";
                    fs << "<text x='" << -margin_left/4 + 5 << "' y='" << y - 3
                       << "' font-family='times' font-size='10' fill='" << plot.font_color_ << "' transform='rotate(-90 "
                       << -margin_left/4 + 8 << " " << y + 5 << ")'>" << std::setprecision(precision) << v << "</text>\n";
                }
            }
            for (CoarseReal x_cord_dataspace : detail::vertical_gridlines(plot.a_, plot.b_, plot.vertical_lines_, plot.log_abscissas_))
            {
                CoarseReal x = x_scale(x_cord_dataspace);
                fs << "<line x1='" << x << "' y1='0' x2='" << x << "' y2='" << function_height
                   << "' stroke='gray' stroke-width='1' opacity='0.5' stroke-dasharray='4' />\n";
                fs << "<text x='" << x - 10 << "' y='" << function_height + 10
                   << "' font-family='times' font-size='10' fill='" << plot.font_color_ << "'>"
                   << detail::axis_label(static_cast<double>(x_cord_dataspace), static_cast<double>(plot.a_), static_cast<double>(plot.b_), plot.vertical_lines_, plot.log_abscissas_) << "</text>\n";
            }
            fs << "<text x='" << graph_width - 5 << "' y='12' font-family='times' font-size='12' text-anchor='end' fill='"
               << plot.font_color_ << "'>" << "f(x)" << "</text>\n";
            // The curve, broken wherever the function is not finite:
            auto draw_curve = [&](auto & canvas) {
                std::vector<detail::plot_point> run;
                for (size_t j = 0; j <= plot.precise_ordinates_.size(); ++j)
                {
                    if (j == plot.precise_ordinates_.size() || !usable(plot.precise_ordinates_[j]))
                    {
                        canvas.polyline(plot.function_color_, run, 1);
                        run.clear();
                        continue;
                    }
                    run.push_back({static_cast<double>(x_scale(plot.coarse_abscissas_[j])),
                                   static_cast<double>(f_scale(transform(plot.precise_ordinates_[j])))});
                }
            };
            if (plot.png_encoder_)
            {
                detail::raster_canvas canvas(-4, -4, graph_width + 8, function_height + 8);
                draw_curve(canvas);
                canvas.finish(fs, plot.png_encoder_, id + "-f");
            }
            else
            {
                detail::vector_canvas canvas(fs, -4, -4, graph_width + 8, function_height + 8, id + "-f");
                draw_curve(canvas);
                canvas.finish();
            }
            fs << "</g>\n";
        }
        fs << "</svg>\n";
        return fs;
    }

private:
    std::vector<PreciseReal> precise_abscissas_;
    std::vector<CoarseReal> coarse_abscissas_;
    std::vector<PreciseReal> precise_ordinates_;
    std::vector<PreciseReal> cond_;
    std::list<std::vector<CoarseReal>> ulp_list_;
    std::list<std::vector<detail::ulps_class>> class_list_;
    std::list<std::string> first_message_;
    std::vector<std::string> colors_;
    CoarseReal a_;
    CoarseReal b_;
    bool log_abscissas_;
    PreciseReal clip_;
    int width_;
    std::string envelope_color_;
    bool ulp_envelope_;
    int horizontal_lines_;
    int vertical_lines_;
    std::string title_;
    std::string background_color_;
    std::string font_color_;
    std::string crop_color_;
    bool crop_color_set_ = false;
    std::string nan_color_;
    bool show_function_ = false;
    bool log_function_ = false;
    std::string function_color_ = "orange";
    png_encoder_type png_encoder_ = detail::default_png_encoder();

    // Distinguishes the plot's data, for svg_id:
    std::string svg_key() const
    {
        std::ostringstream os;
        os << std::setprecision(17) << title_ << ' ' << width_ << ' ' << a_ << ' ' << b_ << ' ' << coarse_abscissas_.size();
        for (auto const & ulps : ulp_list_)
        {
            for (size_t j = 0; j < ulps.size(); j += 1 + ulps.size()/64)
            {
                os << ' ' << ulps[j];
            }
        }
        return os.str();
    }
};

template<class F, typename PreciseReal, typename CoarseReal>
ulps_plot<F, PreciseReal, CoarseReal>& ulps_plot<F, PreciseReal, CoarseReal>::png_encoder(png_encoder_type encoder)
{
    png_encoder_ = encoder;
    return *this;
}

template<class F, typename PreciseReal, typename CoarseReal>
ulps_plot<F, PreciseReal, CoarseReal>& ulps_plot<F, PreciseReal, CoarseReal>::envelope_color(std::string const & color)
{
    envelope_color_ = color;
    return *this;
}

template<class F, typename PreciseReal, typename CoarseReal>
ulps_plot<F, PreciseReal, CoarseReal>& ulps_plot<F, PreciseReal, CoarseReal>::clip(PreciseReal clip)
{
    clip_ = clip;
    return *this;
}

template<class F, typename PreciseReal, typename CoarseReal>
ulps_plot<F, PreciseReal, CoarseReal>& ulps_plot<F, PreciseReal, CoarseReal>::width(int width)
{
    width_ = width;
    return *this;
}

template<class F, typename PreciseReal, typename CoarseReal>
ulps_plot<F, PreciseReal, CoarseReal>& ulps_plot<F, PreciseReal, CoarseReal>::horizontal_lines(int horizontal_lines)
{
    horizontal_lines_ = horizontal_lines;
    return *this;
}

template<class F, typename PreciseReal, typename CoarseReal>
ulps_plot<F, PreciseReal, CoarseReal>& ulps_plot<F, PreciseReal, CoarseReal>::vertical_lines(int vertical_lines)
{
    vertical_lines_ = vertical_lines;
    return *this;
}

template<class F, typename PreciseReal, typename CoarseReal>
ulps_plot<F, PreciseReal, CoarseReal>& ulps_plot<F, PreciseReal, CoarseReal>::show_function(bool show_function)
{
    show_function_ = show_function;
    return *this;
}

template<class F, typename PreciseReal, typename CoarseReal>
ulps_plot<F, PreciseReal, CoarseReal>& ulps_plot<F, PreciseReal, CoarseReal>::log_function(bool log_function)
{
    log_function_ = log_function;
    return *this;
}

template<class F, typename PreciseReal, typename CoarseReal>
ulps_plot<F, PreciseReal, CoarseReal>& ulps_plot<F, PreciseReal, CoarseReal>::function_color(std::string const & color)
{
    function_color_ = color;
    return *this;
}

template<class F, typename PreciseReal, typename CoarseReal>
ulps_plot<F, PreciseReal, CoarseReal>& ulps_plot<F, PreciseReal, CoarseReal>::title(std::string const & title)
{
    title_ = title;
    return *this;
}

template<class F, typename PreciseReal, typename CoarseReal>
ulps_plot<F, PreciseReal, CoarseReal>& ulps_plot<F, PreciseReal, CoarseReal>::background_color(std::string const & background_color)
{
    background_color_ = background_color;
    return *this;
}

template<class F, typename PreciseReal, typename CoarseReal>
ulps_plot<F, PreciseReal, CoarseReal>& ulps_plot<F, PreciseReal, CoarseReal>::font_color(std::string const & font_color)
{
    font_color_ = font_color;
    return *this;
}

template<class F, typename PreciseReal, typename CoarseReal>
ulps_plot<F, PreciseReal, CoarseReal>& ulps_plot<F, PreciseReal, CoarseReal>::crop_color(std::string const & color)
{
    crop_color_ = color;
    crop_color_set_ = true;
    return *this;
}

template<class F, typename PreciseReal, typename CoarseReal>
ulps_plot<F, PreciseReal, CoarseReal>& ulps_plot<F, PreciseReal, CoarseReal>::nan_color(std::string const & color)
{
    nan_color_ = color;
    return *this;
}

template<class F, typename PreciseReal, typename CoarseReal>
ulps_plot<F, PreciseReal, CoarseReal>& ulps_plot<F, PreciseReal, CoarseReal>::ulp_envelope(bool write_ulp_envelope)
{
    ulp_envelope_ = write_ulp_envelope;
    return *this;
}

namespace detail{
inline bool ends_with(std::string const& filename, std::string const& suffix)
{
    if(filename.size() < suffix.size())
    {
        return false;
    }

    return std::equal(std::begin(suffix), std::end(suffix), std::end(filename) - suffix.size());
}
}

template<class F, typename PreciseReal, typename CoarseReal>
std::string ulps_plot<F, PreciseReal, CoarseReal>::summary() const
{
    using std::abs;
    using std::isnan;
    std::ostringstream os;
    os << std::setprecision(std::numeric_limits<CoarseReal>::max_digits10);
    auto uit = ulp_list_.begin();
    auto cit = class_list_.begin();
    auto mit = first_message_.begin();
    for (size_t f = 0; f < ulp_list_.size(); ++f, ++uit, ++cit, ++mit)
    {
        size_t counts[7] = {0, 0, 0, 0, 0, 0, 0};
        size_t first[7];
        for (size_t& v : first) { v = static_cast<size_t>(-1); }
        size_t worst_ulps = static_cast<size_t>(-1);
        size_t worst_ratio = static_cast<size_t>(-1);
        for (size_t j = 0; j < cit->size(); ++j)
        {
            size_t c = static_cast<size_t>((*cit)[j]);
            if (counts[c]++ == 0) { first[c] = j; }
            if ((*cit)[j] != detail::ulps_class::ok) { continue; }
            if (worst_ulps == static_cast<size_t>(-1) || abs((*uit)[j]) > abs((*uit)[worst_ulps])) { worst_ulps = j; }
            if (!isnan(cond_[j]))
            {
                PreciseReal r = abs(static_cast<PreciseReal>((*uit)[j]))/cond_[j];
                if (worst_ratio == static_cast<size_t>(-1) || r > abs(static_cast<PreciseReal>((*uit)[worst_ratio]))/cond_[worst_ratio]) { worst_ratio = j; }
            }
        }
        os << "samples=" << cit->size() << " ok=" << counts[0] << " undefined=" << counts[1] << " overflow_handled=" << counts[2]
           << " underflow_handled=" << counts[3] << " green_nan_or_domain_error=" << counts[4] << " cyan_spurious_overflow=" << counts[5]
           << " mishandled=" << counts[6] << "\n";
        if (worst_ulps != static_cast<size_t>(-1))
        {
            os << "  worst_ulps=" << (*uit)[worst_ulps] << " at x=" << coarse_abscissas_[worst_ulps] << "\n";
        }
        if (worst_ratio != static_cast<size_t>(-1))
        {
            os << "  worst_ulps_over_envelope=" << (static_cast<PreciseReal>((*uit)[worst_ratio])/cond_[worst_ratio]) << " at x=" << coarse_abscissas_[worst_ratio] << "\n";
        }
        char const * names[7] = {"", "", "", "", "green", "cyan", "mishandled"};
        for (size_t c = 4; c < 7; ++c)
        {
            if (counts[c] > 0)
            {
                os << "  " << names[c] << ": first at x=" << coarse_abscissas_[first[c]] << (mit->empty() ? "" : ("; message: " + *mit)) << "\n";
            }
        }
    }
    return os.str();
}

template<class F, typename PreciseReal, typename CoarseReal>
void ulps_plot<F, PreciseReal, CoarseReal>::write(std::string const & filename) const
{
    if(!BOOST_MATH_NAMESPACE::tools::detail::ends_with(filename, ".svg"))
    {
        throw std::logic_error("Only svg files are supported at this time.");
    }
    std::ofstream fs(filename);
    fs << *this;
    fs.close();
}


template<class F, typename PreciseReal, typename CoarseReal>
ulps_plot<F, PreciseReal, CoarseReal>::ulps_plot(F hi_acc_impl, CoarseReal a, CoarseReal b,
             size_t samples, bool perturb_abscissas, int random_seed, bool log_abscissas) : crop_color_("red")
{
    // Use digits10 for this comparison in case the two types have differeing radixes:
    static_assert(std::numeric_limits<PreciseReal>::digits10 >= std::numeric_limits<CoarseReal>::digits10, "PreciseReal must have higher precision that CoarseReal");
    if (samples < 10)
    {
        throw std::domain_error("Must have at least 10 samples, samples = " + std::to_string(samples));
    }
    if (b <= a)
    {
        throw std::domain_error("On interval [a,b], b > a is required.");
    }
    if (log_abscissas && !(a > 0))
    {
        throw std::domain_error("Logarithmic abscissas require a > 0.");
    }
    a_ = a;
    b_ = b;
    log_abscissas_ = log_abscissas;

    using std::exp;
    using std::isfinite;
    using std::log;
    // A uniform grid in the plotted coordinate t (x, or log(x)), at the centres of grid_size cells spanning [a, b], gives the
    // envelope by finite differences and most of the points. A quarter as many random probes (uniform in t) find the errors
    // between the grid points: at integers and half-integers, which often take their own code paths and which a grid can
    // hit or miss, and in narrow features such as the spikes of the condition number at zeros.
    size_t const probes = samples/5;
    size_t const grid_size = samples - probes;
    PreciseReal const t_lo = log_abscissas ? static_cast<PreciseReal>(log(static_cast<PreciseReal>(a))) : static_cast<PreciseReal>(a);
    PreciseReal const t_hi = log_abscissas ? static_cast<PreciseReal>(log(static_cast<PreciseReal>(b))) : static_cast<PreciseReal>(b);
    // With perturb_abscissas, f is evaluated at x and the implementation at the nearest CoarseReal; otherwise both at the CoarseReal.
    auto make_abscissa = [&](PreciseReal const & t, PreciseReal & precise, CoarseReal & coarse) {
        PreciseReal x = log_abscissas ? static_cast<PreciseReal>(exp(t)) : t;
        x = (std::min)((std::max)(x, static_cast<PreciseReal>(a)), static_cast<PreciseReal>(b));
        coarse = static_cast<CoarseReal>(x);
        precise = perturb_abscissas ? x : static_cast<PreciseReal>(coarse);
    };
    std::vector<PreciseReal> grid_x(grid_size);
    std::vector<CoarseReal> grid_coarse(grid_size);
    for (size_t i = 0; i < grid_size; ++i)
    {
        make_abscissa(t_lo + (t_hi - t_lo)*(static_cast<PreciseReal>(i) + PreciseReal(0.5))/static_cast<PreciseReal>(grid_size), grid_x[i], grid_coarse[i]);
    }

    std::mt19937_64 gen;
    if (random_seed == -1)
    {
        std::random_device rd;
        gen.seed(rd());
    }
    else
    {
        gen.seed(static_cast<std::uint64_t>(random_seed));
    }
    // Boost's uniform_real_distribution can generate quad and multiprecision random numbers; std's cannot:
    #ifndef BOOST_MATH_STANDALONE
    boost::random::uniform_real_distribution<PreciseReal> dis(t_lo, t_hi);
    #else
    // Use std::random in standalone mode if it is a type that the standard library can support (float, double, or long double)
    static_assert(std::numeric_limits<PreciseReal>::digits10 <= std::numeric_limits<long double>::digits10, "Standalone mode does not support types with precision that exceeds long double");
    std::uniform_real_distribution<PreciseReal> dis(t_lo, t_hi);
    #endif
    std::vector<PreciseReal> probe_x(probes);
    std::vector<CoarseReal> probe_coarse(probes);
    for (size_t i = 0; i < probes; ++i)
    {
        make_abscissa(dis(gen), probe_x[i], probe_coarse[i]);
    }

    auto reference = [&](PreciseReal const & x) {
#ifndef BOOST_MATH_NO_EXCEPTIONS
        try
        {
#endif
            return static_cast<PreciseReal>(hi_acc_impl(x));
#ifndef BOOST_MATH_NO_EXCEPTIONS
        }
        catch (...)
        {
            return std::numeric_limits<PreciseReal>::quiet_NaN();
        }
#endif
    };
    std::vector<PreciseReal> grid_y(grid_size);
    for (size_t i = 0; i < grid_size; ++i)
    {
        grid_y[i] = reference(grid_x[i]);
    }
    std::vector<PreciseReal> probe_y(probes);
    for (size_t i = 0; i < probes; ++i)
    {
        probe_y[i] = reference(probe_x[i]);
    }
    std::vector<PreciseReal> probe_cond;
    std::vector<PreciseReal> grid_cond = detail::grid_condition_numbers(reference, grid_x, grid_y, probe_x, probe_y, probe_cond);

    // Merge the grid and the probes, in order:
    std::vector<size_t> order(samples);
    for (size_t i = 0; i < samples; ++i)
    {
        order[i] = i;
    }
    auto precise_x = [&](size_t i) -> PreciseReal const & { return i < grid_size ? grid_x[i] : probe_x[i - grid_size]; };
    std::stable_sort(order.begin(), order.end(), [&](size_t i, size_t j) { return precise_x(i) < precise_x(j); });
    precise_abscissas_.resize(samples);
    coarse_abscissas_.resize(samples);
    precise_ordinates_.resize(samples);
    cond_.resize(samples);
    for (size_t k = 0; k < samples; ++k)
    {
        size_t i = order[k];
        bool const grid = i < grid_size;
        size_t const m = grid ? i : i - grid_size;
        precise_abscissas_[k] = grid ? grid_x[m] : probe_x[m];
        coarse_abscissas_[k] = grid ? grid_coarse[m] : probe_coarse[m];
        precise_ordinates_[k] = grid ? grid_y[m] : probe_y[m];
        PreciseReal const c = grid ? grid_cond[m] : probe_cond[m];
        // cond_ is half the condition number, and never below the half-ulp of the correctly rounded result.
        // It is NaN where f is zero or undefined:
        cond_[k] = (isfinite(c) && precise_ordinates_[k] != 0) ? (std::max)(PreciseReal(c/2), PreciseReal(0.5)) : std::numeric_limits<PreciseReal>::quiet_NaN();
    }
    clip_ = -1;
    width_ = 1100;
    envelope_color_ = "chartreuse";
    ulp_envelope_ = true;
    horizontal_lines_ = 8;
    vertical_lines_ = 10;
    title_ = "";
    background_color_ = "black";
    font_color_ = "white";
}

template<class F, typename PreciseReal, typename CoarseReal>
template<class G>
ulps_plot<F, PreciseReal, CoarseReal>& ulps_plot<F, PreciseReal, CoarseReal>::add_fn(G g, std::string const & color)
{
    using std::abs;
    using std::isfinite;
    using std::nextafter;
    size_t samples = precise_abscissas_.size();
    std::vector<CoarseReal> ulps(samples, std::numeric_limits<CoarseReal>::quiet_NaN());
    std::vector<detail::ulps_class> classes(samples, detail::ulps_class::ok);
    std::string first_message;
    for (size_t i = 0; i < samples; ++i)
    {
        PreciseReal y_hi_acc = precise_ordinates_[i];
        CoarseReal w = 0;
        std::string message;
        // The implementation is evaluated everywhere, even where the true value is out of range, to check that
        // overflow and underflow are handled correctly:
        detail::ulps_class c = detail::classify_result<CoarseReal>(y_hi_acc, [&]() { return g(coarse_abscissas_[i]); }, w, &message);
        classes[i] = c;
        if (c != detail::ulps_class::ok)
        {
            if (first_message.empty() && !message.empty() && (c == detail::ulps_class::failed || c == detail::ulps_class::spurious_overflow || c == detail::ulps_class::mishandled))
            {
                first_message = message;
            }
            continue;
        }
        PreciseReal y_lo_acc = static_cast<PreciseReal>(w);
        PreciseReal absy = abs(y_hi_acc);
        PreciseReal dist = static_cast<PreciseReal>(nextafter(static_cast<CoarseReal>(absy), (std::numeric_limits<CoarseReal>::max)()) - static_cast<CoarseReal>(absy));
        ulps[i] = static_cast<CoarseReal>((y_lo_acc - y_hi_acc)/dist);
    }
    class_list_.emplace_back(classes);
    first_message_.emplace_back(first_message);
    ulp_list_.emplace_back(ulps);
    colors_.emplace_back(color);
    return *this;
}




} BOOST_MATH_NAMESPACE_END // namespace boost::math::tools
#endif
