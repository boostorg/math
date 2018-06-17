//  Copyright John Maddock 2016.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_MATH_OVERFLOW_ERROR_POLICY ignore_error

#include <iostream>
#include <iomanip>
#include <numeric>
#include "../stopwatch.hpp"
#include <boost/math/policies/error_handling.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/relative_difference.hpp>

// For the CUDA runtime routines (prefixed with "cuda_")
#include <cuda_runtime.h>
#include <thrust/tabulate.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/random.h>

typedef double float_type;

/**
 * CUDA Kernel Device code
 *
 */

 template <class Real, class F, class Gen>
 struct cuda_naive_monte_carlo_func
 {
    F f;
    const boost::uint64_t seed;
    const Real* p_start_locations;
    const Real*  p_scales;
    Real*  p_storage;
    unsigned n_calls, n_dimentions;
    cuda_naive_monte_carlo_func(F f_, boost::uint64_t seed_, const Real* const&  starts, const Real* const& scales, unsigned calls, Real* const& storage, unsigned dims) 
      : f(f_), seed(seed_), p_start_locations(starts), p_scales(scales), n_calls(calls), p_storage(storage), n_dimentions(dims) {}
    __host__  __device__ Real operator()(unsigned id)
    {
      Gen gen(seed * (id + 1));
      Real* storage_base = p_storage + id * n_dimentions;
      thrust::uniform_real_distribution<Real> dist(0, 1);
      Real sum(0);
      Real c(0);
      for(unsigned i = 0; i < n_calls; ++i)
      {
         for(unsigned j = 0; j < n_dimentions; ++j)
            storage_base[j] = p_start_locations[j] + p_scales[j] * dist(gen);
         Real y = f(storage_base) - c;
         Real t = sum + y;
         c = (t - sum) - y;
         sum = t;
      }
      return sum / n_calls;
    }
 };

 template <class Real, class F, class Gen>
 struct cuda_naive_monte_carlo
 {
   cuda_naive_monte_carlo(const F& integrand,
                      std::vector<std::pair<Real, Real>> const & bounds,
                      uint64_t seed = 0) : m_gen(seed), m_volume(1), m_sigma(0), m_sigma_squares(0)
   {
      auto it = bounds.begin();
      while(it != bounds.end())
      {
         m_start_locations.push_back(it->first);
         m_scale_factors.push_back(it->second - it->first);
         m_volume *= (it->second - it->first);
         ++it;
      }
   }
   /*
   static void throw_on_cuda_error(cudaError_t code, const char *file, int line)
   {
      if (code != cudaSuccess)
      {
         std::stringstream ss;
         ss << file << "(" << line << ")";
         std::string file_and_line;
         ss >> file_and_line;
         throw thrust::system_error(code, thrust::cuda_category(), file_and_line);
      }
   }
   */
   template <class T>
   static T* to_pointer(T* p) { return p; }
   template <class T>
   static T* to_pointer(thrust::device_ptr<T> p) { return p.get(); }
   
   Real integrate(boost::uintmax_t n_calls, boost::uintmax_t calls_per_thread = 1024)
   {
      boost::uintmax_t threads;

      threads = (n_calls % calls_per_thread ? 1 : 0) + n_calls / calls_per_thread;

      thrust::device_vector<Real> starts = m_start_locations;
      thrust::device_vector<Real> scales = m_scale_factors;
      thrust::device_vector<Real> storage(threads * m_start_locations.size());
      cuda_naive_monte_carlo_func<Real, F, Gen> func(m_f, m_gen(), to_pointer(starts.data()), to_pointer(scales.data()), calls_per_thread, to_pointer(storage.data()), scales.size());

      thrust::device_vector<Real> sums(threads);
      thrust::tabulate(sums.begin(), sums.end(), func);
      //thrust::host_vector<Real> hv(sums);
      Real sum = thrust::reduce(sums.begin(), sums.end());
      return m_volume * sum / threads;
   }
   /*
   Real integrate_host(boost::uintmax_t n_calls, boost::uintmax_t calls_per_thread = 1024)
   {
      boost::uintmax_t threads;
      size_t estimated_need, free, total;

      threads = (n_calls % calls_per_thread ? 1 : 0) + n_calls / calls_per_thread;

      thrust::host_vector<Real> starts = m_start_locations;
      thrust::host_vector<Real> scales = m_scale_factors;
      thrust::host_vector<Real> storage(threads * m_start_locations.size());
      cuda_naive_monte_carlo_func<Real, F, Gen> func(m_f, m_gen(), to_pointer(starts.data()), to_pointer(scales.data()), calls_per_thread, to_pointer(storage.data()), scales.size());

      thrust::host_vector<Real> sums(threads);
      std::cout << "Running tabulate:" << std::endl;
      thrust::tabulate(sums.begin(), sums.end(), func);
      std::cout << "Copying vector:" << std::endl;
      thrust::host_vector<Real> host_sums = sums;
      std::cout << "Running reduce:" << std::endl;
      Real sum = thrust::reduce(sums.begin(), sums.end());
      return m_volume * sum / (threads);
   }*/

 private:
   Gen m_gen;
   Real m_volume, m_sigma, m_sigma_squares;
   F m_f;
   thrust::host_vector<Real> m_start_locations;
   thrust::host_vector<Real> m_scale_factors;
 };

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

 void test_host_pi();
 void test_host_hypersphere_10();

/**
 * Host main routine
 */
int main(void)
{
   cudaDeviceProp deviceProp;
   cudaGetDeviceProperties(&deviceProp, 0);

   boost::uintmax_t max_threads = deviceProp.maxThreadsPerMultiProcessor * deviceProp.multiProcessorCount;

   std::cout << "maxThreadsPerBlock = " << deviceProp.maxThreadsPerBlock << std::endl;
   std::cout << "maxThreadsPerMultiProcessor = " << deviceProp.maxThreadsPerMultiProcessor << std::endl;
   std::cout << "multiProcessorCount = " << deviceProp.multiProcessorCount << std::endl << std::endl;
   std::cout << "Total max threads = " << max_threads << std::endl;

   std::cout << "Testing Pi calculation for circle formula.\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
   std::cout << std::right << std::setw(15);
   std::cout << "Threads" << std::right << std::setw(15) << "Points" << std::right << std::setw(15) << "Time" << std::right << std::setw(15)
      << "Value" << std::right << std::setw(15) << "Error" << std::endl;

   for (boost::uintmax_t calls_per_thread = 512; calls_per_thread < 100000; calls_per_thread *= 2)
   {
      for (boost::uintmax_t threads = 1024; threads < max_threads * 2; threads *= 2)
      {
         try {
            std::vector<std::pair<double, double> > bounds = { {-1, 1}, {-1, 1} };

            watch w;
            cuda_naive_monte_carlo<double, pi_calculator<double>, thrust::random::taus88> integrator(pi_calculator<double>(), bounds, 128402);
            double val = integrator.integrate(threads * calls_per_thread, calls_per_thread);
            double elapsed = w.elapsed();
            double err = boost::math::relative_difference(val, boost::math::constants::pi<double>());
            std::cout << std::right << std::setw(15) << std::fixed << threads << std::right << std::setw(15)
               << threads * calls_per_thread << std::right << std::setw(15)
               << elapsed << std::right << std::setw(15)
               << val << std::right << std::setw(15) << std::scientific << err << std::endl;
         }
         catch (const std::exception& e)
         {
            std::cout << "Found exception: " << e.what() << std::endl;
         }
      }
   }

   test_host_pi();

   std::cout << "Testing Hypersphere volume.\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
   std::cout << std::right << std::setw(15);
   std::cout << "Threads" << std::right << std::setw(15) << "Points" << std::right << std::setw(15) << "Time" << std::right << std::setw(15)
      << "Value" << std::right << std::setw(15) << "Error" << std::endl;

   double hypersphere10 = std::pow(boost::math::constants::pi<double>(), 5) / boost::math::tgamma(6.0);

   for (boost::uintmax_t calls_per_thread = 512; calls_per_thread < 50000; calls_per_thread *= 2)
   {
      for (boost::uintmax_t threads = 1024; threads < max_threads * 2; threads *= 2)
      {
         try {
            std::vector<std::pair<double, double> > bounds;
            std::pair<double, double> point = { -1.0, 1.0 };
            for (unsigned i = 0; i < 10; ++i)
               bounds.push_back(point);

            watch w;
            cuda_naive_monte_carlo<double, hypersphere<double, 10>, thrust::random::taus88> integrator(hypersphere<double, 10>(), bounds, 128402);
            double val = integrator.integrate(threads * calls_per_thread, calls_per_thread);
            double elapsed = w.elapsed();
            double err = boost::math::relative_difference(val, hypersphere10);
            std::cout << std::right << std::setw(15) << std::fixed << threads << std::right << std::setw(15)
               << threads * calls_per_thread << std::right << std::setw(15)
               << elapsed << std::right << std::setw(15)
               << val << std::right << std::setw(15) << std::scientific << err << std::endl;
         }
         catch (const std::exception& e)
         {
            std::cout << "Found exception: " << e.what() << std::endl;
         }
      }
   }
   
   test_host_hypersphere_10();

   return 0;
}

