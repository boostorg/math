//  (C) Copyright Nick Thompson 2020.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#ifndef BOOST_MATH_TOOLS_ULP_PLOT_HPP
#define BOOST_MATH_TOOLS_ULP_PLOT_HPP
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
#include <stdexcept>
#include <boost/math/tools/condition_numbers.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/algorithm/string/predicate.hpp>


// Design of this function comes from:
// https://blogs.mathworks.com/cleve/2017/01/23/ulps-plots-reveal-math-function-accurary/

// The envelope is the condition number of function evaluation.

namespace boost::math::tools {

namespace detail {
template<class F1, class F2, class Real>
void write_gridlines(std::ostream& fs, int horizontal_lines, int vertical_lines,
                     F1 x_scale, F2 y_scale, Real min_x, Real max_x, Real min_y, Real max_y,
                     int graph_width, int graph_height, int margin_left, std::string const & font_color)
{
  // Make a grid:
  for (int i = 1; i <= horizontal_lines; ++i) {
      Real y_cord_dataspace = min_y +  ((max_y - min_y)*i)/horizontal_lines;
      auto y = y_scale(y_cord_dataspace);
      fs << "<line x1='0' y1='" << y << "' x2='" << graph_width
         << "' y2='" << y
         << "' stroke='gray' stroke-width='1' opacity='0.5' stroke-dasharray='4' />\n";

      fs << "<text x='" <<  -margin_left/4 + 5 << "' y='" << y - 3
         << "' font-family='times' font-size='10' fill='" << font_color << "' transform='rotate(-90 "
         << -margin_left/4 + 8 << " " << y + 5 << ")'>"
         << std::setprecision(4) << y_cord_dataspace << "</text>\n";
   }

    for (int i = 1; i <= vertical_lines; ++i) {
        Real x_cord_dataspace = min_x +  ((max_x - min_x)*i)/vertical_lines;
        Real x = x_scale(x_cord_dataspace);
        fs << "<line x1='" << x << "' y1='0' x2='" << x
           << "' y2='" << graph_height
           << "' stroke='gray' stroke-width='1' opacity='0.5' stroke-dasharray='4' />\n";

        fs << "<text x='" <<  x - 10  << "' y='" << graph_height + 10
           << "' font-family='times' font-size='10' fill='" << font_color << "'>"
           << std::setprecision(4) << x_cord_dataspace << "</text>\n";
    }
}
}

template<class F, typename PreciseReal, typename CoarseReal>
class ulp_plot {
public:
    ulp_plot(F hi_acc_impl, CoarseReal a, CoarseReal b,
             bool perturb_abscissas = true, size_t samples = 10000, int random_seed = -1);

    void set_clip(int clip);

    void set_width(int width);

    void set_envelope_color(std::string const & color);

    void set_title(std::string const & title);

    void set_background_color(std::string const & background_color);

    void set_font_color(std::string const & font_color);

    void write_ulp_envelope(bool write_ulp);

    template<class G>
    void add_fn(G g, std::string const & color = "steelblue")
    {
        using std::abs;
        size_t samples = precise_abscissas_.size();
        std::vector<CoarseReal> ulps(samples);
        for (size_t i = 0; i < samples; ++i)
        {
            PreciseReal y_hi_acc = precise_ordinates_[i];
            PreciseReal y_lo_acc = g(coarse_abscissas_[i]);
            PreciseReal absy = abs(y_hi_acc);
            PreciseReal dist = nextafter(static_cast<CoarseReal>(absy), std::numeric_limits<CoarseReal>::max()) - static_cast<CoarseReal>(absy);
            ulps[i] = static_cast<CoarseReal>((y_lo_acc - y_hi_acc)/dist);
        }
        ulp_list_.emplace_back(ulps);
        colors_.emplace_back(color);
        return;
    }

    void set_horizontal_lines(int horizontal_lines);

    void set_vertical_lines(int vertical_lines);

    void write(std::string const & filename) const;

    friend std::ostream& operator<<(std::ostream& fs, ulp_plot const & plot)
    {
        using std::abs;
        using std::floor;
        using std::isnan;
        if (plot.ulp_list_.size() == 0)
        {
            throw std::domain_error("No functions added for comparison.");
        }
        if (plot.width_ <= 1)
        {
            throw std::domain_error("Width = " + std::to_string(plot.width_) + ", which is too small.");
        }

        PreciseReal worst_ulp_distance = 0;
        PreciseReal min_y = std::numeric_limits<PreciseReal>::max();
        PreciseReal max_y = std::numeric_limits<PreciseReal>::lowest();
        for (auto const & ulp_vec : plot.ulp_list_)
        {
            for (auto const & ulp : ulp_vec)
            {
                if (abs(ulp) > worst_ulp_distance)
                {
                    worst_ulp_distance = abs(ulp);
                }
                if (ulp < min_y)
                {
                    min_y = ulp;
                }
                if (ulp > max_y)
                {
                    max_y = ulp;
                }
            }
        }

        if (plot.clip_ > 0)
        {
            if (max_y > plot.clip_)
            {
                max_y = plot.clip_;
            }
            if (min_y < -plot.clip_)
            {
                min_y = -plot.clip_;
            }
        }

        int height = floor(double(plot.width_)/1.61803);
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
            return ((x-plot.a_)/(plot.b_ - plot.a_))*static_cast<CoarseReal>(graph_width);
        };

        auto y_scale = [&](PreciseReal y)->PreciseReal
        {
            return ((max_y - y)/(max_y - min_y) )*static_cast<PreciseReal>(graph_height);
        };

        fs << "<?xml version=\"1.0\" encoding='UTF-8' ?>\n"
           << "<svg xmlns='http://www.w3.org/2000/svg' width='"
           << plot.width_ << "' height='"
           << height << "'>\n"
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
                                    static_cast<CoarseReal>(min_y), static_cast<CoarseReal>(max_y), graph_width, graph_height, margin_left, plot.font_color_);
        }
        else
        {
            std::vector<double> ys{-3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0};
            for (size_t i = 0; i < ys.size(); ++i)
            {
                if (min_y <= ys[i] && ys[i] <= max_y)
                {
                    PreciseReal y_cord_dataspace = ys[i];
                    PreciseReal y = y_scale(y_cord_dataspace);
                    fs << "<line x1='0' y1='" << y << "' x2='" << graph_width
                       << "' y2='" << y
                       << "' stroke='gray' stroke-width='1' opacity='0.5' stroke-dasharray='4' />\n";

                    fs << "<text x='" <<  -margin_left/2 << "' y='" << y - 3
                       << "' font-family='times' font-size='10' fill='" << plot.font_color_ << "' transform='rotate(-90 "
                       << -margin_left/2 + 11 << " " << y + 5 << ")'>"
                       <<  std::setprecision(4) << y_cord_dataspace << "</text>\n";
                }
            }
            for (int i = 1; i <= plot.vertical_lines_; ++i)
            {
                CoarseReal x_cord_dataspace = plot.a_ +  ((plot.b_ - plot.a_)*i)/plot.vertical_lines_;
                CoarseReal x = x_scale(x_cord_dataspace);
                fs << "<line x1='" << x << "' y1='0' x2='" << x
                   << "' y2='" << graph_height
                   << "' stroke='gray' stroke-width='1' opacity='0.5' stroke-dasharray='4' />\n";

                fs << "<text x='" <<  x - 10  << "' y='" << graph_height + 10
                   << "' font-family='times' font-size='10' fill='" << plot.font_color_ << "'>"
                   << std::setprecision(4) << x_cord_dataspace << "</text>\n";
            }
        }

        int color_idx = 0;
        for (auto const & ulp : plot.ulp_list_)
        {
            std::string color = plot.colors_[color_idx++];
            for (size_t j = 0; j < ulp.size(); ++j)
            {
                if (isnan(ulp[j]))
                {
                    continue;
                }
                if (plot.clip_ > 0 && abs(ulp[j]) > plot.clip_)
                {
                    continue;
                }
                CoarseReal x = x_scale(plot.coarse_abscissas_[j]);
                PreciseReal y = y_scale(ulp[j]);
                fs << "<circle cx='" << x << "' cy='" << y << "' r='1' fill='" << color << "'/>\n";
            }
        }

        if (plot.ulp_envelope_)
        {
            std::string close_path = "' stroke='"  + plot.envelope_color_ + "' stroke-width='1' fill='none'></path>\n";
            size_t jstart = 0;
            if (plot.clip_ > 0)
            {
                while (plot.cond_[jstart] > plot.clip_)
                {
                    ++jstart;
                    if (jstart >= plot.cond_.size())
                    {
                        goto done;
                    }
                }
            }
            size_t jmin = jstart;
        new_top_path:
            if (jmin >= plot.cond_.size())
            {
                goto start_bottom_paths;
            }
            fs << "<path d='M" << x_scale(plot.coarse_abscissas_[jmin]) << " " << y_scale(plot.cond_[jmin]);

            for (size_t j = jmin + 1; j < plot.coarse_abscissas_.size(); ++j)
            {
                bool bad = isnan(plot.cond_[j]) || (plot.clip_ > 0 && plot.cond_[j] > plot.clip_);
                if (bad)
                {
                    ++j;
                    while ( (j < plot.coarse_abscissas_.size() - 2) && bad)
                    {
                        bad = isnan(plot.cond_[j]) || (plot.clip_ > 0 && plot.cond_[j] > plot.clip_);
                        ++j;
                    }
                    jmin = j;
                    fs << close_path;
                    goto new_top_path;
                }

                CoarseReal t = x_scale(plot.coarse_abscissas_[j]);
                PreciseReal y = y_scale(plot.cond_[j]);
                fs << " L" << t << " " << y;
            }
            fs << close_path;
        start_bottom_paths:
            jmin = jstart;
        new_bottom_path:
            if (jmin >= plot.cond_.size())
            {
                goto done;
            }
            fs << "<path d='M" << x_scale(plot.coarse_abscissas_[jmin]) << " " << y_scale(-plot.cond_[jmin]);

            for (size_t j = jmin + 1; j < plot.coarse_abscissas_.size(); ++j)
            {
                bool bad = isnan(plot.cond_[j]) || (plot.clip_ > 0 && plot.cond_[j] > plot.clip_);
                if (bad)
                {
                    ++j;
                    while ( (j < plot.coarse_abscissas_.size() - 2) && bad)
                    {
                        bad = isnan(plot.cond_[j]) || (plot.clip_ > 0 && plot.cond_[j] > plot.clip_);
                        ++j;
                    }
                    jmin = j;
                    fs << close_path;
                    goto new_bottom_path;
                }
                CoarseReal t = x_scale(plot.coarse_abscissas_[j]);
                PreciseReal y = y_scale(-plot.cond_[j]);
                fs << " L" << t << " " << y;
            }
            fs << close_path;
        }
    done:
        fs << "</g>\n"
           << "</svg>\n";
        return fs;
    }

private:
    std::vector<PreciseReal> precise_abscissas_;
    std::vector<CoarseReal> coarse_abscissas_;
    std::vector<PreciseReal> precise_ordinates_;
    std::vector<PreciseReal> cond_;
    std::list<std::vector<CoarseReal>> ulp_list_;
    std::vector<std::string> colors_;
    CoarseReal a_;
    CoarseReal b_;
    int clip_;
    int width_;
    std::string envelope_color_;
    bool ulp_envelope_;
    int horizontal_lines_;
    int vertical_lines_;
    std::string title_;
    std::string background_color_;
    std::string font_color_;
};

template<class F, typename PreciseReal, typename CoarseReal>
void ulp_plot<F, PreciseReal, CoarseReal>::set_envelope_color(std::string const & color)
{
    envelope_color_ = color;
}

template<class F, typename PreciseReal, typename CoarseReal>
void ulp_plot<F, PreciseReal, CoarseReal>::set_clip(int clip)
{
    clip_ = clip;
}

template<class F, typename PreciseReal, typename CoarseReal>
void ulp_plot<F, PreciseReal, CoarseReal>::set_width(int width)
{
    width_ = width;
}

template<class F, typename PreciseReal, typename CoarseReal>
void ulp_plot<F, PreciseReal, CoarseReal>::set_horizontal_lines(int horizontal_lines)
{
    horizontal_lines_ = horizontal_lines;
}

template<class F, typename PreciseReal, typename CoarseReal>
void ulp_plot<F, PreciseReal, CoarseReal>::set_vertical_lines(int vertical_lines)
{
    vertical_lines_ = vertical_lines;
}

template<class F, typename PreciseReal, typename CoarseReal>
void ulp_plot<F, PreciseReal, CoarseReal>::set_title(std::string const & title)
{
    title_ = title;
}

template<class F, typename PreciseReal, typename CoarseReal>
void ulp_plot<F, PreciseReal, CoarseReal>::set_background_color(std::string const & background_color)
{
    background_color_ = background_color;
}

template<class F, typename PreciseReal, typename CoarseReal>
void ulp_plot<F, PreciseReal, CoarseReal>::set_font_color(std::string const & font_color)
{
    font_color_ = font_color;
}

template<class F, typename PreciseReal, typename CoarseReal>
void ulp_plot<F, PreciseReal, CoarseReal>::write_ulp_envelope(bool write_ulp_envelope)
{
    ulp_envelope_ = write_ulp_envelope;
}

template<class F, typename PreciseReal, typename CoarseReal>
void ulp_plot<F, PreciseReal, CoarseReal>::write(std::string const & filename) const
{
    if (!boost::algorithm::ends_with(filename, ".svg"))
    {
        throw std::logic_error("Only svg files are supported at this time.");
    }
    std::ofstream fs(filename);
    fs << *this;
    fs.close();
}


template<class F, typename PreciseReal, typename CoarseReal>
ulp_plot<F, PreciseReal, CoarseReal>::ulp_plot(F hi_acc_impl, CoarseReal a, CoarseReal b,
             bool perturb_abscissas, size_t samples, int random_seed)
{
    static_assert(sizeof(PreciseReal) >= sizeof(CoarseReal), "PreciseReal must have larger size than CoarseReal");
    if (samples < 10)
    {
        throw std::domain_error("Must have at least 10 samples, samples = " + std::to_string(samples));
    }
    if (b <= a)
    {
        throw std::domain_error("On interval [a,b], b > a is required.");
    }
    a_ = a;
    b_ = b;

    std::mt19937_64 gen;
    if (random_seed == -1)
    {
        std::random_device rd;
        gen.seed(rd());
    }
    // Boost's uniform_real_distribution can generate quad and multiprecision random numbers; std's cannot:
    boost::random::uniform_real_distribution<PreciseReal> dis(a, b);
    precise_abscissas_.resize(samples);
    coarse_abscissas_.resize(samples);

    if (perturb_abscissas)
    {
        for(size_t i = 0; i < samples; ++i)
        {
            precise_abscissas_[i] = dis(gen);
        }
        std::sort(precise_abscissas_.begin(), precise_abscissas_.end());
        for (size_t i = 0; i < samples; ++i)
        {
            coarse_abscissas_[i] = static_cast<CoarseReal>(precise_abscissas_[i]);
        }
    }
    else
    {
        for(size_t i = 0; i < samples; ++i)
        {
            coarse_abscissas_[i] = static_cast<CoarseReal>(dis(gen));
        }
        std::sort(coarse_abscissas_.begin(), coarse_abscissas_.end());
        for (size_t i = 0; i < samples; ++i)
        {
            precise_abscissas_[i] = coarse_abscissas_[i];
        }
    }

    precise_ordinates_.resize(samples);
    for (size_t i = 0; i < samples; ++i)
    {
        precise_ordinates_[i] = hi_acc_impl(precise_abscissas_[i]);
    }

    cond_.resize(samples, std::numeric_limits<PreciseReal>::quiet_NaN());
    for (size_t i = 0 ; i < samples; ++i)
    {
        PreciseReal y = precise_ordinates_[i];
        if (y != 0)
        {
            // Maybe cond_ is badly names; should it be half_cond_?
            cond_[i] = boost::math::tools::evaluation_condition_number(hi_acc_impl, precise_abscissas_[i])/2;
            // Half-ULP accuracy is the correctly rounded result, so make sure the envelop doesn't go below this:
            if (cond_[i] < 0.5)
            {
                cond_[i] = 0.5;
            }
        }
        // else leave it as nan.
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



} // namespace boost::math::tools
#endif
