//  Copyright John Maddock 2016.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_MATH_OVERFLOW_ERROR_POLICY ignore_error

#ifdef __CUDACC__
#define BOOST_PP_VARIADICS 0
#endif

#include <iostream>
#include <iomanip>
#include <numeric>
#include "../stopwatch.hpp"
#include <boost/math/policies/error_handling.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/relative_difference.hpp>
#include <boost/math/special_functions/pow.hpp>
#include <boost/math/quadrature/cuda_naive_monte_carlo.hpp>


typedef double float_type;

 template <class Real>
 struct pi_calculator
 {
    __host__ __device__ Real operator()(const Real* abscissa)const
    {
       return abscissa[0] * abscissa[0] + abscissa[1] * abscissa[1] <= 1 ? 1 : 0;
    }
 };

 template <class Real, unsigned N>
 struct hypersphere
 {
    __host__ __device__ Real operator()(const Real* abscissa)const
    {
       Real location = 0;
       for (unsigned i = 0; i < N; ++i)
          location += abscissa[i] * abscissa[i];
       return location <= 1 ? 1 : 0;
    }
 };

 void test_host_pi(double);
 void test_host_hypersphere_10();

/**
 * Host main routine
 */
int main(void)
{

   using boost::math::quadrature::cuda_naive_monte_carlo;
   cudaDeviceProp deviceProp;
   cudaGetDeviceProperties(&deviceProp, 0);

   boost::uintmax_t max_threads = deviceProp.maxThreadsPerMultiProcessor * deviceProp.multiProcessorCount;

   std::cout << "maxThreadsPerBlock = " << deviceProp.maxThreadsPerBlock << std::endl;
   std::cout << "maxThreadsPerMultiProcessor = " << deviceProp.maxThreadsPerMultiProcessor << std::endl;
   std::cout << "multiProcessorCount = " << deviceProp.multiProcessorCount << std::endl << std::endl;
   std::cout << "Total max threads = " << max_threads << std::endl;

   std::cout << "Testing Pi calculation for circle formula.\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
   std::cout << std::right << std::setw(15);
   std::cout << "Threads" << std::right << std::setw(15) << "Init #/thread" 
      << std::right << std::setw(15) << "Total Points" 
      << std::right << std::setw(15) << "Time" 
      << std::right << std::setw(15) << "Error Goal" 
      << std::right << std::setw(15) << "Variance"
      << std::right << std::setw(15) << "Error Est."
      << std::right << std::setw(15) << "Error Actual" << std::endl;

   //
   // Do something to initialize the CUDA device:
   //
   do {
      std::vector<std::pair<double, double> > bounds = { { -1, 1 },{ -1, 1 } };
      cuda_naive_monte_carlo<double, pi_calculator<double>, thrust::random::taus88> init(pi_calculator<double>(), bounds);
      init.integrate(0.01);
   } while (0);


   for (double error_goal = 1e-3; error_goal > 1e-5; error_goal /= 2)
   {
      for (boost::uintmax_t calls_per_thread = 128; calls_per_thread < 2048; calls_per_thread *= 2)
      {
         try {
            std::vector<std::pair<double, double> > bounds = { {-1, 1}, {-1, 1} };

            watch w;
            cuda_naive_monte_carlo<double, pi_calculator<double>, thrust::random::taus88> integrator(pi_calculator<double>(), bounds);
            double val = integrator.integrate(error_goal, calls_per_thread);
            double elapsed = w.elapsed();
            double err = fabs(val - boost::math::constants::pi<double>());
            std::cout << std::right << std::setw(15) << std::fixed << max_threads 
               << std::right << std::setw(15) << calls_per_thread 
               << std::right << std::setw(15) << integrator.calls() 
               << std::right << std::setw(15) << std::fixed << elapsed 
               << std::right << std::setw(15) << std::scientific << error_goal
               << std::right << std::setw(15) << std::scientific << integrator.variance()
               << std::right << std::setw(15) << std::scientific << integrator.current_error_estimate()
               << std::right << std::setw(15) << std::scientific << err << std::endl;
         }
         catch (const std::exception& e)
         {
            std::cout << "Found exception: " << e.what() << std::endl;
         }
      }
   }

   for (double error_goal = 1e-3; error_goal > 1e-5; error_goal /= 2)
   {
      test_host_pi(error_goal);
   }

   std::cout << "Testing Hypersphere volume.\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
   std::cout << std::right << std::setw(15);
   std::cout << "Threads" << std::right << std::setw(15) << "Init #/thread"
      << std::right << std::setw(15) << "Total Points"
      << std::right << std::setw(15) << "Time"
      << std::right << std::setw(15) << "Error Goal"
      << std::right << std::setw(15) << "Variance"
      << std::right << std::setw(15) << "Error Est."
      << std::right << std::setw(15) << "Error Actual" << std::endl;

   double hypersphere10 = std::pow(boost::math::constants::pi<double>(), 5) / boost::math::tgamma(6.0);

   // initialized CUDA device code:
   {
      std::vector<std::pair<double, double> > bounds;
      std::pair<double, double> point = { -1.0, 1.0 };
      for (unsigned i = 0; i < 10; ++i)
         bounds.push_back(point);

      cuda_naive_monte_carlo<double, hypersphere<double, 10>, thrust::random::taus88> integrator(hypersphere<double, 10>(), bounds);
      double val = integrator.integrate(1e-2);
   }

   for (boost::uintmax_t calls_per_thread = 64; calls_per_thread < 4086; calls_per_thread *= 2)
   {
      try {
         double error_goal = 0.02;

         std::vector<std::pair<double, double> > bounds;
         std::pair<double, double> point = { -1.0, 1.0 };
         for (unsigned i = 0; i < 10; ++i)
            bounds.push_back(point);

         watch w;
         cuda_naive_monte_carlo<double, hypersphere<double, 10>, thrust::random::taus88> integrator(hypersphere<double, 10>(), bounds);
         double val = integrator.integrate(error_goal, calls_per_thread);
         double elapsed = w.elapsed();
         double err = fabs(val - hypersphere10);
         std::cout << std::right << std::setw(15) << std::fixed << max_threads
            << std::right << std::setw(15) << calls_per_thread
            << std::right << std::setw(15) << integrator.calls()
            << std::right << std::setw(15) << std::fixed << elapsed
            << std::right << std::setw(15) << std::scientific << error_goal
            << std::right << std::setw(15) << std::scientific << integrator.variance()
            << std::right << std::setw(15) << std::scientific << integrator.current_error_estimate()
            << std::right << std::setw(15) << std::scientific << err << std::endl;
      }
      catch (const std::exception& e)
      {
         std::cout << "Found exception: " << e.what() << std::endl;
      }
   }
   
   test_host_hypersphere_10();

   /* Example code from docs */

   {

      // Define a function to integrate:
      auto g = [] __device__ (const double* x)
      {
         constexpr const double pi = boost::math::constants::pi<double>();
         constexpr const double A = 1.0 / (pi * pi * pi);
         return A / (1.0 - cos(x[0])*cos(x[1])*cos(x[2]));
      };
      std::vector<std::pair<double, double>> bounds{ { 0, boost::math::constants::pi<double>() },{ 0, boost::math::constants::pi<double>() },{ 0, boost::math::constants::pi<double>() } };
      double error_goal = 0.001;
      cuda_naive_monte_carlo<double, decltype(g)> mc(g, bounds);

      double result = mc.integrate(error_goal);

      std::cout << "Integration result is: " << result << std::endl;


   }

   return 0;
}

