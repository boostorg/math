// Copyright Paul A. Bristow 2017
// Copyright John Z. Maddock 2017

// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or
//  copy at http ://www.boost.org/LICENSE_1_0.txt).

/*! \brief Graph showing use of Lambert W function.

\details

*/

#include <boost/math/special_functions/lambert_w.hpp>
using boost::math::lambert_w0;
using boost::math::lambert_wm1;

#include "c:\Users\Paul\Desktop\lambert_w0.hpp" // for jm_tests
using boost::math::jm_lambert_w0;

#include <boost/math/special_functions.hpp>
using boost::math::isfinite;
#include <boost/svg_plot/svg_2d_plot.hpp>
using namespace boost::svg;

#include <iostream>
// using std::cout;
// using std::endl;
#include <exception>
#include <stdexcept>
#include <string>
#include <array>
#include <vector>
#include <utility>
using std::pair;
#include <map>
using std::map;
#include <set>
using std::multiset;
#include <limits>
using std::numeric_limits;
#include <cmath> //


  /*!
  */
int main()
{
  try
  {
    std::cout << "Lambert W graph example." << std::endl;

//[lambert_w_graph_1
//] [/lambert_w_graph_1]

    std::map<const double, double> w0s;   // Lambert W0 branch values.
    std::map<const double, double> wm1s;   // Lambert W-1 branch values.
    std::map<const double, double> w0s_big;   // Lambert W0 branch values for large z and W.
    std::map<const double, double> wm1s_big;   // Lambert W-1 branch values for small z and large -W.

    std::cout.precision(std::numeric_limits<double>::max_digits10);

    int count = 0;
    for (double z = -0.36787944117144232159552377016146086744581113103176804; z < 2.8; z += 0.001)
    {
      double w0 = jm_lambert_w0(z);
      //double w0 = lambert_w0(z);
      w0s[z] = w0;
 //     std::cout << "z " << z << ", w = " << w0 << std::endl;
      count++;
    }
    std::cout << "points " << count << std::endl;

    count = 0;
    for (double z = -0.3678794411714423215955237701614608727; z  < -0.001; z += 0.001)
    {
      double wm1 = lambert_wm1(z);
      wm1s[z] = wm1;
      count++;
    }
    std::cout << "points " << count << std::endl;

    svg_2d_plot data_plot;
    svg_2d_plot data_plot2;

    data_plot.title("Lambert W function.")
      .x_size(400)
      .y_size(300)
      .legend_on(true)
      .x_label("z")
      .y_label("W")
      //.x_label_on(true)
      //.y_label_on(true)
      //.xy_values_on(false)
      .x_range(-1, 3.)
      .y_range(-4., +1.)
      .x_major_interval(1.)
      .y_major_interval(1.)
      .x_major_grid_on(true)
      .y_major_grid_on(true)
      //.x_values_on(true)
      //.y_values_on(true)
      .y_values_rotation(horizontal)
      //.plot_window_on(true)
      .x_values_precision(3)
      .y_values_precision(3)
      .coord_precision(4) // Needed to avoid stepping on curves.
      .copyright_holder("Paul A. Bristow")
      .copyright_date("2017")
      //.background_border_color(black);
      ;
    
    // bigger W
    for (double z = -0.3678794411714423215955237701614608727; z < 10000.; z += 50.)
    {
      double w0 = lambert_w0(z);
      w0s_big[z] = w0;
      count++;
    }
    std::cout << "points " << count << std::endl;

    count = 0;
    for (double z = -0.3678794411714423215955237701614608727; z < -0.001; z += 0.001)
    {
      double wm1 = lambert_wm1(z);
      wm1s_big[z] = wm1;
      count++;
    }
    std::cout << "points " << count << std::endl;

    data_plot.plot(w0s, "W0 branch").line_color(red).shape(none).line_on(true).bezier_on(false).line_width(1);
    data_plot.plot(wm1s, "W-1 branch").line_color(blue).shape(none).line_on(true).bezier_on(false).line_width(1);
    data_plot.write("./lambert_w_graph");

    data_plot2.title("Lambert W function for larger z.")
      .x_size(400)
      .y_size(300)
      .legend_on(true)
      .x_label("z")
      .y_label("W")
      //.x_label_on(true)
      //.y_label_on(true)
      //.xy_values_on(false)
      .x_range(-1, 10000.)
      .y_range(-1., +8.)
      .x_major_interval(2000.)
      .y_major_interval(1.)
      .x_major_grid_on(true)
      .y_major_grid_on(true)
      //.x_values_on(true)
      //.y_values_on(true)
      .y_values_rotation(horizontal)
      //.plot_window_on(true)
      .x_values_precision(3)
      .y_values_precision(3)
      .coord_precision(4) // Needed to avoid stepping on curves.
      .copyright_holder("Paul A. Bristow")
      .copyright_date("2017")
      //.background_border_color(black);
      ;

    data_plot2.plot(w0s_big, "W0 branch").line_color(red).shape(none).line_on(true).bezier_on(false).line_width(1);
    data_plot2.plot(wm1s_big, "W-1 branch").line_color(blue).shape(none).line_on(true).bezier_on(false).line_width(1);
    data_plot2.write("./lambert_w_graph_big_w");

    // bezier_on(true); // ???
  }
  catch (std::exception& ex)
  {
    std::cout << ex.what() << std::endl;
  }
}  // int main()

   /*

   //[lambert_w_graph_1_output

   //] [/lambert_w_graph_1_output]
   */
