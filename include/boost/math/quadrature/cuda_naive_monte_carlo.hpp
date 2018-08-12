/*
 * Copyright Nick Thompson, 2018
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef BOOST_MATH_QUADRATURE_NAIVE_MONTE_CARLO_HPP
#define BOOST_MATH_QUADRATURE_NAIVE_MONTE_CARLO_HPP
#include <algorithm>
#include <vector>
#include <utility>
#include <random>
#include <boost/math/policies/error_handling.hpp>
#include <boost/cstdint.hpp>
 // For the CUDA runtime routines (prefixed with "cuda_")
#include <cuda_runtime.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/random.h>

namespace boost { namespace math { namespace quadrature { namespace detail{

   template <class T>
   struct pair_add
   {
      __host__  __device__  thrust::pair<T, T> operator()(const thrust::pair<T, T>& a, const thrust::pair<T, T>& b)const
      {
         return thrust::pair<T, T>(a.first + b.first, a.second + b.second);
      }
   };

   template <class F, class Gen, class Real>
   void __global__ cuda_naive_monte_carlo_device_proc(F f, Gen* seeds, thrust::pair<Real, Real>* sums, const Real* p_start_locations, const Real* p_scales, unsigned n_calls, Real* p_storage, unsigned n_dimentions, bool is_first)
   {
      int id = blockIdx.x * blockDim.x + threadIdx.x;

      Gen gen(seeds[id]);
      Real* storage_base = p_storage + id * n_dimentions;
      thrust::uniform_real_distribution<Real> dist(0, 1);
      Real sum(0);
      Real sigma_squared(0);
      Real c(0);
      for (unsigned i = 0; i < n_calls; ++i)
      {
         for (unsigned j = 0; j < n_dimentions; ++j)
            storage_base[j] = p_start_locations[j] + p_scales[j] * dist(gen);
         Real fv = (f)(storage_base);
         Real y = fv - c;
         Real t = sum + y;
         c = (t - sum) - y;
         sum = t;
         sigma_squared += fv * fv;
      }
      seeds[id] = gen;
      if (is_first)
         sums[id] = thrust::pair<Real, Real>(sum, sigma_squared);
      else
      {
         sums[id].first += sum;
         sums[id].second += sigma_squared;
      }
   }

   template <class F, class Gen, class Real>
   void __global__ cuda_naive_monte_carlo_fast_device_proc(F f, Gen* seeds, thrust::pair<Real, Real>* sums, const Real* p_start_locations, const Real* p_scales, unsigned n_calls, Real* p_storage, unsigned n_dimentions, bool is_first)
   {
      int id = blockIdx.x * blockDim.x + threadIdx.x;

      Gen gen(seeds[id]);
      Real* storage_base = p_storage + id * n_dimentions;
      thrust::uniform_real_distribution<Real> dist(0, 1);
      Real sum(0);
      Real sigma_squared(0);
      for (unsigned i = 0; i < n_calls; ++i)
      {
         for (unsigned j = 0; j < n_dimentions; ++j)
            storage_base[j] = p_start_locations[j] + p_scales[j] * dist(gen);
         Real fv = (f)(storage_base);
         sum += fv;
         sigma_squared += fv * fv;
      }
      seeds[id] = gen;
      if (is_first)
         sums[id] = thrust::pair<Real, Real>(sum, sigma_squared);
      else
      {
         sums[id].first += sum;
         sums[id].second += sigma_squared;
      }
   }

} // namespace detail

   template <class Real, class F, class ThreadGen = thrust::random::taus88, class MasterGen = std::random_device>
   struct cuda_naive_monte_carlo
   {
      cuda_naive_monte_carlo(const F& integrand,
         std::vector<std::pair<Real, Real>> const & bounds,
         const MasterGen& seed) : m_f(integrand), m_gen(seed), m_volume(1), m_sigma(0), m_sigma_squares(0), m_calls(0)
      {
         auto it = bounds.begin();
         while (it != bounds.end())
         {
            m_start_locations.push_back(it->first);
            m_scale_factors.push_back(it->second - it->first);
            m_volume *= (it->second - it->first);
            ++it;
         }
      }
      cuda_naive_monte_carlo(const F& integrand,
         std::vector<std::pair<Real, Real>> const & bounds) 
         : m_f(integrand), m_volume(1), m_sigma(0), m_sigma_squares(0), m_calls(0)
      {
         auto it = bounds.begin();
         while (it != bounds.end())
         {
            m_start_locations.push_back(it->first);
            m_scale_factors.push_back(it->second - it->first);
            m_volume *= (it->second - it->first);
            ++it;
         }
      }

      template <class T>
      static T* to_pointer(T* p) { return p; }
      template <class T>
      static T* to_pointer(thrust::device_ptr<T> p) { return p.get(); }

      Real integrate(Real error_request, boost::uintmax_t calls_per_thread = 1024, boost::uintmax_t max_calls_per_thread = 250000, bool is_compensated = true)
      {
         boost::uintmax_t threads;
         cudaDeviceProp deviceProp;
         cudaGetDeviceProperties(&deviceProp, 0);

         threads = deviceProp.maxThreadsPerMultiProcessor * deviceProp.multiProcessorCount;

         thrust::device_vector<Real> starts = m_start_locations;
         thrust::device_vector<Real> scales = m_scale_factors;
         thrust::device_vector<Real> storage(threads * m_start_locations.size());
         thrust::device_vector<thrust::pair<Real, Real> > sums(threads);
         thrust::device_vector<ThreadGen> seeds(threads);
         thrust::host_vector<ThreadGen> host_seeds(threads);
         typedef typename ThreadGen::result_type seed_type;
         std::uniform_int_distribution<seed_type> ui_dist((std::numeric_limits<seed_type>::min)(), (std::numeric_limits<seed_type>::max)());
         for (unsigned i = 0; i < host_seeds.size(); ++i)
            host_seeds[i].seed(ui_dist(m_gen));
         seeds = host_seeds;
         bool first_call = true;
         bool have_variance = false;
         Real sample_variance = 0;
         do {
            if (m_calls)
            {
               // If we've been called before, adjust calls_per_thread according
               // to the sample variance:
               Real err = sqrt(variance() / m_calls);
               boost::uintmax_t target_total = variance() > 0 ? (boost::uintmax_t)std::floor(1.05 * m_calls * err * err / (error_request * error_request)) : 2 * m_calls;
               calls_per_thread = 1 + target_total / threads;
               have_variance = true;
               sample_variance = variance();
            }
            do
            {
               if (!first_call)
               {
                  // Update how many calls per thread to try based on our variance estimate:
                  boost::uintmax_t target_total;
                  if (sample_variance > 0)
                     target_total = static_cast<boost::uintmax_t>(1.05 * m_calls * boost::math::pow<2>(sqrt(sample_variance / m_calls)) / (error_request * error_request) - m_calls);
                  else
                     target_total = 2 * m_calls;
                  calls_per_thread = 1 + target_total / threads;
               }
               if (calls_per_thread > max_calls_per_thread)
                  calls_per_thread = max_calls_per_thread;
               // std::cout << "Executing with calls_per_thread = " << calls_per_thread << std::endl;
               if(is_compensated)
                  detail::cuda_naive_monte_carlo_device_proc << <threads / 256, 256 >> > (m_f, to_pointer(seeds.data()), to_pointer(sums.data()), to_pointer(starts.data()), to_pointer(scales.data()), calls_per_thread, to_pointer(storage.data()), scales.size(), first_call);
               else
                  detail::cuda_naive_monte_carlo_fast_device_proc << <threads / 256, 256 >> > (m_f, to_pointer(seeds.data()), to_pointer(sums.data()), to_pointer(starts.data()), to_pointer(scales.data()), calls_per_thread, to_pointer(storage.data()), scales.size(), first_call);
               first_call = false;
               m_calls += threads * calls_per_thread;
               // If we haven't been called before then get an estimate of the sample 
               // variance based on the first sum, this removes the need to reduce
               // the whole thing which is relatively expensive:
               if (!have_variance)
               {
                  std::pair<Real, Real> first_result;
                  cudaMemcpy(&first_result, to_pointer(sums.data()), sizeof(first_result), cudaMemcpyDeviceToHost);
                  // It's tempting to use calls_per_thread as the sample size for the sums: but we may
                  // have been through this loop more than once and calls_per_thread is just the "extra"
                  // calls from the last run through.  So divide the total number of calls by the number
                  // of threads to get the actual sample size for each sub-sum:
                  sample_variance = m_volume * m_volume * (first_result.second - first_result.first * first_result.first / (m_calls / threads)) / ((m_calls / threads) - 1);
               }
               // std::cout << "Estimated error = " << sqrt(sample_variance / m_calls) << std::endl;
            } while ((sample_variance < 0) || (sqrt(sample_variance / m_calls) > error_request));

            // Reduce the results:
            thrust::pair<Real, Real> sum = thrust::reduce(sums.begin(), sums.end(), thrust::pair<Real, Real>(0, 0), detail::pair_add<Real>());
            m_sigma += sum.first;
            m_sigma_squares += sum.second;
            first_call = true;
            //
            // Since we used an estimate for the variance in the code above, now that
            // we have the true variance after a full reduction, we had better double check
            // the termination condition:
            //
            // std::cout << "Estimated error after reduce = " << sqrt(variance() / m_calls) << std::endl;
         } while ((variance() < 0) || (sqrt(variance() / m_calls) > error_request));

         return m_volume * m_sigma / m_calls;
      }

      Real variance()const
      {
         return m_volume * m_volume * (m_sigma_squares - m_sigma * m_sigma / m_calls) / (m_calls - 1);
      }
      Real current_error_estimate() const
      {
         using std::sqrt;
         return sqrt(variance() / m_calls);
      }
      uint64_t calls() const
      {
         return m_calls;
      }

   private:
      MasterGen m_gen;
      Real m_volume, m_sigma, m_sigma_squares;
      boost::uintmax_t m_calls;
      F m_f;
      thrust::host_vector<Real> m_start_locations;
      thrust::host_vector<Real> m_scale_factors;
   };

}}}
#endif
