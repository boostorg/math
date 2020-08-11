//  (C) Copyright Nick Thompson 2020.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#ifndef BOOST_MATH_TOOLS_QQ_PLOT_HPP
#define BOOST_MATH_TOOLS_QQ_PLOT_HPP
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
#include <filesystem>
#include <functional>



namespace boost::math::statistics {

/*namespace detail {
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
}*/

template<typename Real>
class qq_plot {
public:
    qq_plot(std::function<Real(Real)> f1, std::function<Real(Real)> f2);

    qq_plot& width(int width);

    qq_plot& title(std::string const & title);

    qq_plot& background_color(std::string const & background_color);

    qq_plot& font_color(std::string const & font_color);

    qq_plot& horizontal_lines(int horizontal_lines);

    qq_plot& vertical_lines(int vertical_lines);

    qq_plot& xlabel(std::string const & label);

    qq_plot& ylabel(std::string const & label);

    qq_plot& line_color(std::string const & line_color);

    void write(std::filesystem::path const & filename) const;

    friend std::ostream& operator<<(std::ostream& fs, qq_plot const & plot)
    {
        using std::abs;
        using std::floor;
        using std::isnan;
        if (plot.width_ <= 1)
        {
            throw std::domain_error("Width = " + std::to_string(plot.width_) + ", which is too small.");
        }

        int height = static_cast<int>(floor(double(plot.width_)/1.61803));
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
        /*auto x_scale = [&](Real x)->Real
        {
            return ((x-plot.a_)/(plot.b_ - plot.a_))*static_cast<Real>(graph_width);
        };

        auto y_scale = [&](Real y)->Real
        {
            return ((max_y - y)/(max_y - min_y) )*static_cast<Real>(graph_height);
        };*/

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
        /*fs << "<g transform='translate(" << margin_left << ", " << margin_top << ")'>\n";
            // y-axis:
        fs  << "<line x1='0' y1='0' x2='0' y2='" << graph_height
            << "' stroke='gray' stroke-width='1'/>\n";
        Real x_axis_loc = y_scale(static_cast<Real>(0));
        fs << "<line x1='0' y1='" << x_axis_loc
            << "' x2='" << graph_width << "' y2='" << x_axis_loc
            << "' stroke='gray' stroke-width='1'/>\n";

        fs << "</g>\n"
           << "</svg>\n";*/
        return fs;
    }

private:
    int width_;
    int horizontal_lines_;
    int vertical_lines_;
    std::string title_;
    std::string background_color_;
    std::string font_color_;
    std::string line_color_;
    std::string xlabel_;
    std::string ylabel_;
};


template<typename Real>
qq_plot<Real>& qq_plot<Real>::width(int width)
{
    width_ = width;
    return *this;
}

template<typename Real>
qq_plot<Real>& qq_plot<Real>::horizontal_lines(int horizontal_lines)
{
    horizontal_lines_ = horizontal_lines;
    return *this;
}

template<typename Real>
qq_plot<Real>& qq_plot<Real>::vertical_lines(int vertical_lines)
{
    vertical_lines_ = vertical_lines;
    return *this;
}

template<typename Real>
qq_plot<Real>& qq_plot<Real>::title(std::string const & title)
{
    title_ = title;
    return *this;
}

template<typename Real>
qq_plot<Real>& qq_plot<Real>::background_color(std::string const & background_color)
{
    background_color_ = background_color;
    return *this;
}

template<typename Real>
qq_plot<Real>& qq_plot<Real>::font_color(std::string const & font_color)
{
    font_color_ = font_color;
    return *this;
}

template<typename Real>
qq_plot<Real>& qq_plot<Real>::line_color(std::string const & line_color)
{
    line_color_ = line_color;
    return *this;
}


template<typename Real>
void qq_plot<Real>::write(std::filesystem::path const & filename) const
{
    if (filename.empty())
    {
        throw std::logic_error("Filename is empty.");
    }
    if (filename.extension() != ".svg")
    {
        throw std::logic_error("Only svg files are supported at this time.");
    }
    std::ofstream fs(filename);
    fs << *this;
    fs.close();
}


template<typename Real>
qq_plot<Real>::qq_plot(std::function<Real(Real)> f1, std::function<Real(Real)> f2)
{
    width_ = 1100;
    line_color_ = "chartreuse";
    horizontal_lines_ = 8;
    vertical_lines_ = 10;
    title_ = "";
    background_color_ = "black";
    font_color_ = "white";
    xlabel_ = "";
    ylabel_ = "";
}

} // namespace boost::math::statistics
#endif
