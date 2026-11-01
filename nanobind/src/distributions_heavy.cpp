// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// Stable-law and other distributions with large tables or unusual constructors, the empirical cdf, and
// find_location and find_scale, over the location-scale families.

#include <boost/math/distributions/holtsmark.hpp>
#include <boost/math/distributions/landau.hpp>
#include <boost/math/distributions/mapairy.hpp>
#include <boost/math/distributions/saspoint5.hpp>
#include <boost/math/distributions/von_mises.hpp>
#include <boost/math/distributions/kolmogorov_smirnov.hpp>
#include <boost/math/distributions/hyperexponential.hpp>
#include <boost/math/distributions/cauchy.hpp>
#include <boost/math/distributions/empirical_cumulative_distribution_function.hpp>
#include <boost/math/distributions/extreme_value.hpp>
#include <boost/math/distributions/find_location.hpp>
#include <boost/math/distributions/find_scale.hpp>
#include <boost/math/distributions/laplace.hpp>
#include <boost/math/distributions/logistic.hpp>
#include <boost/math/distributions/normal.hpp>

#include <nanobind/stl/vector.h>

#include <cmath>
#include <type_traits>
#include <vector>

#include "bindings.hpp"
#include "distributions.hpp"

namespace bm = boost::math;
namespace nb = nanobind;
using nb::arg;
using bmpy::add_dist_methods;
using pol = bmpy::dist_policy;

namespace {

// f(std::type_identity<D>()) for the bound class `type` among the Ds.
template <class... D>
struct families {
   template <class F>
   static double call(char const* name, nb::handle type, F f) {
      double r = 0;
      if (!((type.is(nb::type<D>()) && (r = bmpy::call_scalar(name, [&] { return f(std::type_identity<D>()); }), true)) || ...))
         throw nb::type_error((std::string(name) + ": the distribution must be one of the location-scale classes "
                               "normal, cauchy, laplace, logistic, extreme_value, holtsmark, mapairy or saspoint5").c_str());
      return r;
   }
};
using location_scale = families<bm::normal_distribution<double, pol>, bm::cauchy_distribution<double, pol>, bm::laplace_distribution<double, pol>,
                                bm::logistic_distribution<double, pol>, bm::extreme_value_distribution<double, pol>,
                                bm::holtsmark_distribution<double, pol>, bm::mapairy_distribution<double, pol>, bm::saspoint5_distribution<double, pol>>;

using ecdf = bm::empirical_cumulative_distribution_function<std::vector<double>>;

}  // namespace

void bmpy::bind_distributions_heavy(nb::module_& m) {
   {
      using D = bm::holtsmark_distribution<double, pol>;
      nb::class_<D> c(m, "holtsmark_distribution",
                     "Holtsmark distribution, the stable law with alpha = 3/2 and beta = 0, shifted by location and scaled by scale." BMPY_DIST_URL("holtsmark_dist"));
      c.def(bmpy::init<double, double>(), arg("location") = 0.0, arg("scale") = 1.0);
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("location", location), BMPY_PARAM("scale", scale)});
   }
   {
      using D = bm::landau_distribution<double, pol>;
      nb::class_<D> c(m, "landau_distribution",
                     "Landau distribution, the stable law with alpha = 1 and beta = 1, of energy loss in a thin layer." BMPY_DIST_URL("landau_dist"));
      c.def(bmpy::init<double, double>(), arg("location") = 0.0, arg("scale") = 1.0);
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("location", location), BMPY_PARAM("scale", scale)});
   }
   {
      using D = bm::mapairy_distribution<double, pol>;
      nb::class_<D> c(m, "mapairy_distribution",
                     "Map-Airy distribution, the stable law with alpha = 3/2 and beta = 1, shifted by location and scaled by scale." BMPY_DIST_URL("mapairy_dist"));
      c.def(bmpy::init<double, double>(), arg("location") = 0.0, arg("scale") = 1.0);
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("location", location), BMPY_PARAM("scale", scale)});
   }
   {
      using D = bm::saspoint5_distribution<double, pol>;
      nb::class_<D> c(m, "saspoint5_distribution",
                     "S-alpha-S point5 distribution, the stable law with alpha = 1/2 and beta = 0, shifted by location and scaled by scale." BMPY_DIST_URL("saspoint5_dist"));
      c.def(bmpy::init<double, double>(), arg("location") = 0.0, arg("scale") = 1.0);
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("location", location), BMPY_PARAM("scale", scale)});
   }
   {
      using D = bm::von_mises_distribution<double, pol>;
      nb::class_<D> c(m, "von_mises_distribution",
                     "Von Mises distribution on the circle: pdf exp(concentration cos(x - mean)) / (2 pi I0(concentration))." BMPY_DIST_URL("von_mises_dist"));
      c.def(bmpy::init<double, double>(), arg("mean") = 0.0, arg("concentration") = 1.0);
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("mean", mean), BMPY_PARAM("concentration", concentration)});
   }
   {
      using D = bm::kolmogorov_smirnov_distribution<double, pol>;
      nb::class_<D> c(m, "kolmogorov_smirnov_distribution",
                     "Kolmogorov-Smirnov distribution, in its limiting form, of D_n = sup_x |G(x) - G_n(x)|, G_n the empirical cdf of n observations." BMPY_DIST_URL("kolmogorov_smirnov_dist"));
      c.def(bmpy::init<double>(), arg("n"));
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("n", number_of_observations)});
   }
   {
      using D = bm::hyperexponential_distribution<double, pol>;
      nb::class_<D> c(m, "hyperexponential_distribution",
                     "Hyperexponential distribution: the mixture sum_i probabilities[i] rates[i] exp(-rates[i] x); probabilities default to equal." BMPY_DIST_URL("hyperexponential_dist"));
      c.def(bmpy::init<const std::vector<double>&>(), arg("rates"));
      c.def(bmpy::init<const std::vector<double>&, const std::vector<double>&>(), arg("probabilities"), arg("rates"));
      add_dist_methods(c);
      bmpy::add_dist_params(c, {BMPY_PARAM("probabilities", probabilities), BMPY_PARAM("rates", rates)});
   }
   // After the classes they take.
   m.def("find_location", [](nb::handle distribution, double z, double p, double scale) {
      return location_scale::call("find_location", distribution, [&](auto d) { return bm::find_location<typename decltype(d)::type>(z, p, scale, pol()); });
   }, arg(), arg(), arg(), arg(), nb::sig("def find_location(distribution, z, p, scale, /)"),
      "The location for which the distribution class, with the given scale, has cdf(z) = p: z - scale quantile(p) of the "
      "standard distribution." BMPY_URL("dist_ref/dist_algorithms"));
   m.def("find_scale", [](nb::handle distribution, double z, double p, double location) {
      return location_scale::call("find_scale", distribution, [&](auto d) { return bm::find_scale<typename decltype(d)::type>(z, p, location, pol()); });
   }, arg(), arg(), arg(), arg(), nb::sig("def find_scale(distribution, z, p, location, /)"),
      "The scale for which the distribution class, with the given location, has cdf(z) = p: (z - location) / quantile(p) "
      "of the standard distribution." BMPY_URL("dist_ref/dist_algorithms"));

   nb::class_<ecdf>(m, "empirical_cumulative_distribution_function",
                    "The empirical cdf of a sample: F(x) = (the number of samples <= x) / n. The constructor copies the data, "
                    "a 1D array, and sorts it unless sorted is True." BMPY_DIST_URL("empirical_cdf"))
      .def("__init__", [](ecdf* self, nb::handle data, bool sorted) {
         nb::object a = bmpy::numpy().attr("asarray")(data, arg("dtype") = bmpy::numpy().attr("float64"));
         if (nb::cast<int>(a.attr("ndim")) != 1)
            throw nb::value_error("empirical_cumulative_distribution_function: data must be one-dimensional");
         auto view = nb::cast<nb::ndarray<const double, nb::ndim<1>, nb::device::cpu>>(a);
         std::vector<double> v(view.shape(0));
         for (std::size_t i = 0; i < v.size(); ++i) {
            v[i] = view(i);
            if (std::isnan(v[i]))
               throw nb::value_error("empirical_cumulative_distribution_function: data contains NaN");
         }
         new (self) ecdf(std::move(v), sorted);
      }, arg(), nb::kw_only(), arg("sorted") = false, nb::sig("def __init__(self, data, /, *, sorted=False)"))
      .def("__call__", [](ecdf const& F, nb::handle x, int threads) {
         auto f = [&F]<class T>(double x) { return std::isnan(x) ? x : F(x); };
         return bmpy::call_array<false, bmpy::Real>("empirical_cumulative_distribution_function", f, {x}, {"x"}, threads);
      }, arg(), nb::kw_only(), arg("threads") = 1, nb::sig(bmpy::vec_signature("__call__", std::array<const char*, 1>{"x"}, true).c_str()),
         "F(x), elementwise.");
}
