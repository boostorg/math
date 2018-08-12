

#include <iostream>
#include <iomanip>
#include <random>
#include <boost/random.hpp>
#include "../stopwatch.hpp"
#include <boost/math/policies/error_handling.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/quadrature/naive_monte_carlo.hpp>



void test_host_pi(double error_goal)
{
   using boost::math::quadrature::naive_monte_carlo;


   watch w;
   auto g = [](std::vector<double> const & x)->double
   {
      double r = x[0] * x[0] + x[1] * x[1];
      if (r <= 1)
      {
         return 1;
      }
      return 0;
   };

   std::cout << "Regular host code:\n";

   std::vector<std::pair<double, double>> bounds2{ { -1.0, 1.0 },{ -1., 1. } };

   for (unsigned threads = 1; threads < 9; ++threads)
   {
      w.reset();
      naive_monte_carlo<double, decltype(g), boost::random::taus88> mc(g, bounds2, error_goal,
         /*singular =*/ false,/* threads = */ threads, /* seed = */ 128402);
      auto task = mc.integrate();
      double val = task.get();
      double elapsed = w.elapsed();
      boost::uintmax_t points = mc.calls();
      double err = fabs(val - boost::math::constants::pi<double>());
      std::cout << std::right << std::setw(15) << threads 
         << std::right << std::setw(15) << "-" 
         << std::right << std::setw(15) << points 
         << std::right << std::setw(15) << std::fixed << elapsed 
         << std::right << std::setw(15) << std::scientific << error_goal 
         << std::right << std::setw(15) << std::scientific << mc.variance()
         << std::right << std::setw(15) << std::scientific << mc.current_error_estimate()
         << std::right << std::setw(15) << std::scientific << err << std::endl;
   }
}

void test_host_hypersphere_10()
{
   using boost::math::quadrature::naive_monte_carlo;


   watch w;
   auto g = [](std::vector<double> const & x)->double
   {
      double location = 0;
      for (unsigned i = 0; i < 10; ++i)
         location += x[i] * x[i];
      return location <= 1 ? 1 : 0;
   };

   std::cout << "Regular host code:\n";

   double hypersphere10 = std::pow(boost::math::constants::pi<double>(), 5) / boost::math::tgamma(6.0);

   std::vector<std::pair<double, double>> bounds2;
   std::pair<double, double> point = { -1.0, 1.0 };
   for (unsigned i = 0; i < 10; ++i)
      bounds2.push_back(point);
   double error_goal = 0.02;
   for (unsigned threads = 1; threads < 9; ++threads)
   {
      w.reset();
      naive_monte_carlo<double, decltype(g), boost::random::taus88> mc(g, bounds2, error_goal,
         /*singular =*/ false,/* threads = */ threads, /* seed = */ 128402);
      auto task = mc.integrate();
      double val = task.get();
      double elapsed = w.elapsed();
      boost::uintmax_t points = mc.calls();
      double err = fabs(val - hypersphere10);
      std::cout << std::right << std::setw(15) << threads
         << std::right << std::setw(15) << "-"
         << std::right << std::setw(15) << points
         << std::right << std::setw(15) << std::fixed << elapsed
         << std::right << std::setw(15) << std::scientific << error_goal
         << std::right << std::setw(15) << std::scientific << mc.variance()
         << std::right << std::setw(15) << std::scientific << mc.current_error_estimate()
         << std::right << std::setw(15) << std::scientific << err << std::endl;
   }
}
